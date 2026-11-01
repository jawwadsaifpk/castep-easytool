# -*- coding: utf-8 -*-
"""
Adaptive cut-off energy and k-point convergence for CASTEP.

Runs CASTEP job by job, extending the cut-off energy and then the k-point mesh until
the total energy (and, for optimisation tasks, the stress) changes by less than a
threshold for ``require_n_conv`` consecutive steps. Resumable: finished jobs are
detected and skipped.

The k-point mesh is lattice-aware (see kmesh.py), so it works for any cell shape, and
the report gives both the converged MP grid and an equivalent KPOINTS_MP_SPACING.

Based on the run_adaptive_conv.py workflow; the .conv config format stays compatible.
"""
import math
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from . import cellio, kmesh

# Keys and defaults. Types are inferred from the defaults; ``None`` defaults are "auto".
DEFAULTS = {
    "castep_command": "castep <seedname>",
    "seedname": "",
    "n_atoms": 0,             # 0 = count from the .cell file
    "threshold_mev": 1.0,     # meV/atom
    "require_n_conv": 2,
    "ecut_start": 400,
    "ecut_step": 100,
    "ecut_max": 1200,
    "kspacing_fixed": 0.07,   # 1/A: moderate fixed mesh during the cut-off scan
    "kpn_fixed": None,        # explicit 'n1 n2 n3' overrides kspacing_fixed; None = auto
    "kpn_start": 1,           # N = points along the densest axis (others scale with the cell)
    "kpn_step": 1,
    "kpn_max": 0,             # 0/auto = until kspacing_min
    "kspacing_min": 0.02,     # 1/A: finest spacing considered when kpn_max is auto
    "stress": "auto",         # auto|true|false: add a stress criterion (optimisation tasks)
    "threshold_stress_gpa": 0.1,
    "final_check": True,      # verify the result one cut-off step higher
    "reuse": True,            # chain REUSE .check files along the cut-off scan
}
_AUTO_WORDS = ("auto", "none", "")
_PSEUDO_SUFFIXES = (".usp", ".recpot", ".uspcc", ".upf", ".uspso", ".pspot")
COLOR_E, COLOR_RATE, COLOR_S = "#2E75B6", "#C00000", "#2E8B57"


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def _to_bool(s):
    v = s.strip().lower()
    if v in ("true", "yes", "1", "on"):
        return True
    if v in ("false", "no", "0", "off"):
        return False
    raise ValueError("expected true/false")


def default_config(seed=""):
    cfg = {k: (list(v) if isinstance(v, list) else v) for k, v in DEFAULTS.items()}
    cfg["seedname"] = seed
    return cfg


def parse_config(path):
    """Read a .conv file. Returns (cfg, warnings). Raises ValueError on bad values."""
    cfg = default_config()
    warns = []
    for n, raw in enumerate(Path(path).read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        key, val = line.split(":", 1)
        key, val = key.strip().lower(), val.strip()
        if key == "kpn_ratio":
            warns.append(f"line {n}: 'kpn_ratio' is obsolete (the k-mesh now follows the cell shape automatically)")
            continue
        if key not in cfg:
            warns.append(f"line {n}: unknown setting '{key}' ignored")
            continue
        default = DEFAULTS[key]
        try:
            if key == "kpn_fixed":
                if val.lower() in _AUTO_WORDS:
                    cfg[key] = None
                else:
                    nums = [int(x) for x in val.split()]
                    if len(nums) != 3 or min(nums) < 1:
                        raise ValueError("needs three positive integers or 'auto'")
                    cfg[key] = nums
            elif key == "kpn_max":
                cfg[key] = 0 if val.lower() in _AUTO_WORDS else int(float(val))
            elif key == "stress":
                v = val.lower()
                cfg[key] = "auto" if v == "auto" else ("true" if _to_bool(v) else "false")
            elif isinstance(default, bool):
                cfg[key] = _to_bool(val)
            elif isinstance(default, int):
                cfg[key] = int(float(val))
            elif isinstance(default, float):
                cfg[key] = float(val)
            else:
                cfg[key] = val
        except ValueError as e:
            raise ValueError(f"line {n}: bad value for '{key}': {val!r} ({e})")
    if not cfg["seedname"]:
        raise ValueError("'seedname' is not set in the config file")
    if cfg["ecut_step"] <= 0 or cfg["kpn_step"] <= 0:
        raise ValueError("ecut_step and kpn_step must be positive")
    if cfg["ecut_start"] > cfg["ecut_max"]:
        raise ValueError("ecut_start is larger than ecut_max")
    if cfg["kpn_start"] < 1 or (cfg["kpn_max"] and cfg["kpn_start"] > cfg["kpn_max"]):
        raise ValueError("kpn_start must be at least 1 and not larger than kpn_max")
    if cfg["require_n_conv"] < 1:
        raise ValueError("require_n_conv must be at least 1")
    if cfg["kspacing_fixed"] <= 0 or cfg["kspacing_min"] <= 0:
        raise ValueError("kspacing_fixed and kspacing_min must be positive")
    return cfg, warns


def render_conv(cfg, notes=None):
    """Text of a .conv file for ``cfg``. ``notes`` maps a key to an extra comment line."""
    notes = notes or {}
    seed = cfg["seedname"]

    def line(key, value, comment=None):
        out = [f"# {c}" for c in (comment.split("\n") if comment else [])]
        if key in notes:
            out.append(f"# {notes[key]}")
        out.append(f"{key:<22}: {value}")
        return "\n".join(out)

    kfix = "auto" if cfg["kpn_fixed"] is None else " ".join(map(str, cfg["kpn_fixed"]))
    kmax = "auto" if not cfg["kpn_max"] else cfg["kpn_max"]
    parts = [
        "# ============================================================\n"
        f"# castep-easytool adaptive convergence config - {seed}\n"
        f"# Run:  castep-easytool autoconv {seed}.conv\n"
        f"# Needs {seed}.cell and {seed}.param next to this file.\n"
        "# ============================================================",
        line("castep_command", cfg["castep_command"],
             "CASTEP command; <seedname> is replaced by each job name.\n"
             "Serial: castep <seedname>   Parallel: mpirun -np 24 castep.mpi <seedname>"),
        line("seedname", seed),
        line("n_atoms", cfg["n_atoms"], "Atoms in the cell, for meV/atom normalisation (0 = count from the .cell)."),
        "# ---- Convergence criteria ----",
        line("threshold_mev", cfg["threshold_mev"], "Energy tolerance (meV/atom)."),
        line("require_n_conv", cfg["require_n_conv"], "Consecutive steps that must all be below the tolerance."),
        line("stress", cfg["stress"],
             "Also require the stress to converge: auto (on for geometry optimisation with a free cell) | true | false."),
        line("threshold_stress_gpa", cfg["threshold_stress_gpa"], "Largest allowed change of any stress component (GPa)."),
        line("final_check", str(cfg["final_check"]).lower(),
             "Verify the result with one extra job one cut-off step higher."),
        "# ---- Phase 1: cut-off scan ----",
        line("ecut_start", cfg["ecut_start"], "Cut-off energy range (eV)."),
        line("ecut_step", cfg["ecut_step"]),
        line("ecut_max", cfg["ecut_max"], "Safety cap: the scan stops here even if not converged."),
        line("kspacing_fixed", cfg["kspacing_fixed"],
             "k-point spacing (1/Angstrom, CASTEP convention) used during the cut-off scan."),
        line("kpn_fixed", kfix, "Explicit 'n1 n2 n3' grid for the cut-off scan (overrides kspacing_fixed), or auto."),
        line("reuse", str(cfg["reuse"]).lower(), "Chain REUSE .check files along the cut-off scan (faster)."),
        "# ---- Phase 2: k-point scan ----",
        line("kpn_start", cfg["kpn_start"],
             "N = k-points along the axis with the longest reciprocal vector; the other axes scale with\n"
             "the cell shape (N x N x N for cubic/fcc/bcc cells, e.g. N x N x N/2 for a tall cell)."),
        line("kpn_step", cfg["kpn_step"]),
        line("kpn_max", kmax, "Last N, or auto = until the spacing reaches kspacing_min."),
        line("kspacing_min", cfg["kspacing_min"], "Finest spacing (1/Angstrom) considered when kpn_max is auto."),
    ]
    return "\n\n".join(parts) + "\n"


# ---------------------------------------------------------------------------
# CASTEP file helpers
# ---------------------------------------------------------------------------

_NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"


def parse_energy(path):
    """Last 'Final energy, E' in a .castep file (falls back to the free energy), or None."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for pat in (rf"Final energy, E\s*=\s*({_NUM})\s*eV", rf"Final free energy \(E-TS\)\s*=\s*({_NUM})\s*eV"):
        m = re.findall(pat, text)
        if m:
            return float(m[-1].replace("D", "E"))
    return None


def parse_stress(path):
    """
    Last stress tensor of a .castep file as [xx, yy, zz, xy, xz, yz] in GPa, or None.
    Handles both 'Symmetrised Stress Tensor' and 'Stress Tensor' boxes.
    """
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    heads = list(re.finditer(r"\*+\s*(?:Symmetrised\s+)?Stress Tensor\s*\*+", text))
    if not heads:
        return None
    block = text[heads[-1].end(): heads[-1].end() + 1500]
    rows = re.findall(rf"^\s*\*\s*[xyz]\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*\*", block, re.MULTILINE)
    if len(rows) < 3:
        return None
    m = np.array([[float(x) for x in r] for r in rows[:3]])
    return [float(m[0, 0]), float(m[1, 1]), float(m[2, 2]), float(m[0, 1]), float(m[0, 2]), float(m[1, 2])]


def pressure_from_stress(s):
    return None if s is None else -(s[0] + s[1] + s[2]) / 3.0


def is_complete(path):
    """True if CASTEP finished (timing summary written) and produced an energy."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return bool(re.search(r"^\s*Total time\s*=", text, re.MULTILINE)) and "Final energy" in text


def count_atoms_in_cell(cell_text):
    """Number of atoms in the cell text (0 if it cannot be read)."""
    try:
        return len(cellio.read_cell(cell_text)["species"])
    except (ValueError, KeyError):
        return 0


def write_cell(template, dest, kgrid):
    """Write the .cell with exactly one k-point specification: KPOINTS_MP_GRID."""
    out = re.sub(r"^\s*kpoints?_mp_(?:grid|spacing|offset)\b.*$\n?", "", template,
                 flags=re.IGNORECASE | re.MULTILINE)
    out = re.sub(r"%BLOCK\s+KPOINTS?_LIST.*?%ENDBLOCK\s+KPOINTS?_LIST\s*\n?", "", out,
                 flags=re.IGNORECASE | re.DOTALL)
    out = out.rstrip() + f"\n\nKPOINTS_MP_GRID : {kgrid[0]} {kgrid[1]} {kgrid[2]}\n"
    Path(dest).write_text(out, encoding="utf-8")


def _set_keyword(text, key, value):
    """Replace the whole line of ``key`` (case-insensitive) or append it."""
    pat = re.compile(rf"^\s*{key}\s*[:=]?.*$", re.IGNORECASE | re.MULTILINE)
    line = f"{key} : {value}"
    if pat.search(text):
        return pat.sub(line, text)
    return text.rstrip() + f"\n{line}\n"


def write_param(template, dest, ecut, reuse=None, stress=False):
    """
    Write the .param for one job: SinglePoint, finite-basis correction off (it would
    hide the cut-off error), cut-off set.
    """
    out = re.sub(r"^\s*basis_precision\b.*\n?", "", template, flags=re.IGNORECASE | re.MULTILINE)
    out = _set_keyword(out, "task", "SinglePoint")
    out = _set_keyword(out, "cut_off_energy", f"{float(ecut):.1f} eV")
    out = _set_keyword(out, "finite_basis_corr", "0")
    if stress:
        out = _set_keyword(out, "calculate_stress", "true")
    if reuse:
        out = _set_keyword(out, "reuse", reuse)
    else:
        out = re.sub(r"^\s*reuse\s*[:=].*\n?", "", out, flags=re.IGNORECASE | re.MULTILINE)
    Path(dest).write_text(out, encoding="utf-8")


def detect_xc(param_text):
    m = re.search(r"^\s*xc_functional\s*[:=]\s*(\S+)", param_text, re.IGNORECASE | re.MULTILINE)
    return m.group(1).upper() if m else ""


def resolve_stress(cfg, param_text, cell_info):
    """Whether the stress criterion is active, and a short reason."""
    mode = cfg["stress"]
    if mode == "true":
        return True, "requested"
    if mode == "false":
        return False, "disabled"
    task = re.search(r"^\s*task\s*[:=]\s*(\S+)", param_text, re.IGNORECASE | re.MULTILINE)
    task = task.group(1).lower() if task else ""
    if cellio.is_fixed_cell(cell_info):
        return False, "cell is fixed (FIX_ALL_CELL)"
    if task.startswith("geom"):
        return True, "geometry optimisation with a free cell"
    return False, f"task is {task or 'unspecified'}, not a geometry optimisation"


# ---------------------------------------------------------------------------
# Running + convergence
# ---------------------------------------------------------------------------

def split_command(command, jobname):
    return shlex.split(command.replace("<seedname>", jobname), posix=(sys.platform != "win32"))


def run_castep(jobname, jobdir, command, extra_flag=None, quiet=False, timeout=None):
    """Run one CASTEP job in ``jobdir``; True on exit code 0."""
    token = f"{extra_flag} {jobname}" if extra_flag else jobname
    cmd = split_command(command, token)
    if not quiet:
        print(f"  -> Running: {' '.join(cmd)}", flush=True)
    t0 = time.time()
    try:
        ret = subprocess.call(cmd, cwd=str(jobdir), timeout=timeout,
                              stdout=subprocess.DEVNULL if quiet else None,
                              stderr=subprocess.DEVNULL if quiet else None)
    except FileNotFoundError:
        print(f"    X Command not found: '{cmd[0]}'. Check castep_command (use --cmd).")
        return False
    except subprocess.TimeoutExpired:
        print(f"    X Timed out after {timeout} s")
        return False
    if not quiet:
        print(f"    {'OK' if ret == 0 else 'X Exit code ' + str(ret)} after {time.time() - t0:.1f} s")
    return ret == 0


def check_convergence(energies, stresses, threshold_mev, n_atoms, require_n=2, thr_stress=None):
    """
    (converged, dE, dS): dE[i] = |E[i+1]-E[i]| in meV/atom; dS[i] = largest change of any
    stress component in GPa (empty if stresses are not used). Converged when the last
    ``require_n`` steps all satisfy dE < threshold (and dS < thr_stress when given).
    """
    if len(energies) < 2:
        return False, [], []
    dE = [abs(energies[i + 1] - energies[i]) * 1000 / n_atoms for i in range(len(energies) - 1)]
    dS = []
    if thr_stress is not None and len(stresses) == len(energies):
        dS = [float(np.max(np.abs(np.array(stresses[i + 1]) - np.array(stresses[i]))))
              for i in range(len(stresses) - 1)]
    if len(dE) < require_n:
        return False, dE, dS
    ok = all(d < threshold_mev for d in dE[-require_n:])
    if thr_stress is not None and dS:
        ok = ok and all(d < thr_stress for d in dS[-require_n:])
    return ok, dE, dS


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_plot(ecut_data, kpn_data, n_atoms, cfg, use_stress, out_paths, seed="", xc="", suffix=""):
    """Rows: cut-off scan and k-point scan; columns: E - E_ref, |dE|, and |d sigma| if stress is used."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    thr, thr_s = cfg["threshold_mev"], cfg["threshold_stress_gpa"]
    ncol = 3 if use_stress else 2
    with plt.rc_context({"font.size": 12, "axes.labelsize": 13, "axes.titlesize": 14,
                         "legend.fontsize": 10, "axes.linewidth": 1.2,
                         "xtick.direction": "in", "ytick.direction": "in"}):
        fig, axes = plt.subplots(2, ncol, figsize=(6 * ncol, 8), squeeze=False)

        def rate(ax, xs, ys, color, thr_val, unit, title):
            ax.semilogy(xs, np.maximum(ys, 1e-12), "s-", color=color, lw=2, ms=7)
            ax.axhline(thr_val, color="gray", ls="--", lw=1.2, label=f"{thr_val:g} {unit} threshold")
            ax.set(ylabel=f"|Δ| ({unit})", title=title)
            ax.legend()
            ax.grid(True, alpha=0.3)

        def pending(row, text):
            for ax in axes[row]:
                ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes,
                        fontsize=13, color="gray")
                ax.set_axis_off()

        d = ecut_data
        if len(d["E"]) >= 2:
            xs, Es = np.array(d["x"]), np.array(d["E"])
            ax = axes[0][0]
            ax.plot(xs, Es - Es[-1], "o-", color=COLOR_E, lw=2, ms=7)
            ax.axhline(0, color="gray", lw=0.8, ls="--")
            ax.set(xlabel="Cut-off energy (eV)", ylabel="E − E$_{ref}$ (eV)", title="Cut-off energy convergence")
            ax.grid(True, alpha=0.3)
            rate(axes[0][1], xs[1:], np.abs(np.diff(Es)) * 1000 / n_atoms, COLOR_RATE, thr, "meV/atom", "Energy change")
            axes[0][1].set_xlabel("Cut-off energy (eV)")
            if use_stress and all(s is not None for s in d["S"]):
                S = np.array(d["S"])
                rate(axes[0][2], xs[1:], np.max(np.abs(np.diff(S, axis=0)), axis=1), COLOR_S, thr_s, "GPa", "Stress change")
                axes[0][2].set_xlabel("Cut-off energy (eV)")
        else:
            pending(0, "Cut-off data pending…")

        d = kpn_data
        if len(d["E"]) >= 2:
            xk, Ek = np.arange(len(d["labels"])), np.array(d["E"])
            tick = [f"{l}\nΔ={s:.3f}" for l, s in zip(d["labels"], d["spacing"])]
            rot = 30 if len(xk) > 5 else 0
            ax = axes[1][0]
            ax.plot(xk, Ek - Ek[-1], "o-", color=COLOR_E, lw=2, ms=7)
            ax.axhline(0, color="gray", lw=0.8, ls="--")
            ax.set_xticks(xk)
            ax.set_xticklabels(tick, rotation=rot, fontsize=9)
            ax.set(xlabel="k-point grid (Δ = spacing, 1/Å)", ylabel="E − E$_{ref}$ (eV)", title="k-point convergence")
            ax.grid(True, alpha=0.3)
            rate(axes[1][1], xk[1:], np.abs(np.diff(Ek)) * 1000 / n_atoms, COLOR_RATE, thr, "meV/atom", "Energy change")
            axes[1][1].set_xticks(xk[1:])
            axes[1][1].set_xticklabels(tick[1:], rotation=rot, fontsize=9)
            axes[1][1].set_xlabel("k-point grid")
            if use_stress and all(s is not None for s in d["S"]):
                S = np.array(d["S"])
                rate(axes[1][2], xk[1:], np.max(np.abs(np.diff(S, axis=0)), axis=1), COLOR_S, thr_s, "GPa", "Stress change")
                axes[1][2].set_xticks(xk[1:])
                axes[1][2].set_xticklabels(tick[1:], rotation=rot, fontsize=9)
                axes[1][2].set_xlabel("k-point grid")
        else:
            pending(1, "k-point data pending…")

        xc_txt = f" ({xc})" if xc else ""
        fig.suptitle(f"{seed} convergence tests{xc_txt}{suffix}\n{datetime.now():%Y-%m-%d %H:%M}",
                     fontsize=13, y=1.0)
        fig.tight_layout()
        for p in out_paths:
            fig.savefig(str(p), dpi=300 if str(p).endswith("png") else None, bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

def _fmt_grid(g):
    return f"{g[0]}x{g[1]}x{g[2]}"


def _ecut_values(cfg):
    return list(range(cfg["ecut_start"], cfg["ecut_max"] + 1, cfg["ecut_step"]))


def _load_templates(cfg, conv_dir):
    seed = cfg["seedname"]
    cell_path, param_path = Path(conv_dir) / f"{seed}.cell", Path(conv_dir) / f"{seed}.param"
    for p in (cell_path, param_path):
        if not p.exists():
            raise FileNotFoundError(f"'{p.name}' not found in {conv_dir} (the seedname must match; "
                                    f"`castep-easytool autoconv <structure.cif>` can create all inputs)")
    cell_text = cell_path.read_text(encoding="utf-8", errors="replace")
    param_text = param_path.read_text(encoding="utf-8", errors="replace")
    try:
        info = cellio.read_cell(cell_text)
    except ValueError as e:
        raise ValueError(f"cannot read the lattice/positions of '{cell_path.name}': {e}")
    return cell_text, param_text, info


def plan(cfg, conv_dir, phase="both"):
    """Everything needed to describe/run the scans, without running anything."""
    cell_text, param_text, info = _load_templates(cfg, conv_dir)
    lat = info["lattice"]
    n_atoms = cfg["n_atoms"] or len(info["species"])
    use_stress, why = resolve_stress(cfg, param_text, info)
    fixed = cfg["kpn_fixed"] or kmesh.grid_from_spacing(lat, cfg["kspacing_fixed"])
    steps = kmesh.ladder(lat, cfg["kpn_start"], cfg["kpn_step"], cfg["kpn_max"], cfg["kspacing_min"])
    return {"cell_text": cell_text, "param_text": param_text, "info": info, "lattice": lat,
            "n_atoms": n_atoms, "use_stress": use_stress, "stress_reason": why,
            "fixed_grid": fixed, "k_steps": steps, "ecuts": _ecut_values(cfg), "phase": phase}


def plan_lines(cfg, pl, conv_dir):
    out = [f"Seed            : {cfg['seedname']}  (inputs in {conv_dir})",
           f"Atoms           : {pl['n_atoms']}",
           "Reciprocal cell : |b|/2pi = " + ", ".join(f"{x:.4f}" for x in kmesh.recip_lengths(pl['lattice'])) + " 1/A",
           f"Criteria        : dE < {cfg['threshold_mev']:g} meV/atom"
           + (f" and d(stress) < {cfg['threshold_stress_gpa']:g} GPa" if pl["use_stress"] else "")
           + f", {cfg['require_n_conv']} consecutive steps",
           f"Stress criterion: {'ON' if pl['use_stress'] else 'off'} ({pl['stress_reason']})",
           f"CASTEP command  : {cfg['castep_command']}"]
    ev = pl["ecuts"]
    if pl["phase"] in ("both", "ecut"):
        out.append(f"Cut-off scan    : {ev[0]}..{ev[-1]} eV, step {cfg['ecut_step']}, fixed k-grid "
                   f"{_fmt_grid(pl['fixed_grid'])}, up to {len(ev)} jobs"
                   + ("  (REUSE chained)" if cfg["reuse"] else ""))
    if pl["phase"] in ("both", "kpoints"):
        ks = pl["k_steps"]
        out.append(f"k-point scan    : {_fmt_grid(ks[0][1])} (spacing {ks[0][2]:.3f}) .. "
                   f"{_fmt_grid(ks[-1][1])} (spacing {ks[-1][2]:.3f} 1/A), up to {len(ks)} jobs")
    if pl["phase"] == "both":
        out.append(f"Worst case      : {len(ev) + len(pl['k_steps']) + (1 if cfg['final_check'] else 0)} CASTEP jobs")
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_adaptive(cfg, conv_dir, workdir=None, phase="both", ecut_fixed=None, plot_only=False):
    """
    Execute the convergence workflow. Returns a dict with the recommended parameters,
    the report text and the problems found. Raises FileNotFoundError/ValueError for
    setup problems.
    """
    conv_dir = Path(conv_dir).resolve()
    seed = cfg["seedname"]
    workdir = Path(workdir) if workdir else conv_dir / "serial"
    plots_dir = conv_dir / "plots"
    pl = plan(cfg, conv_dir, phase)
    lat, n_atoms, use_stress = pl["lattice"], pl["n_atoms"], pl["use_stress"]
    cell_tmpl, param_tmpl = pl["cell_text"], pl["param_text"]
    xc = detect_xc(param_tmpl)
    thr, req, thr_s = cfg["threshold_mev"], cfg["require_n_conv"], cfg["threshold_stress_gpa"]
    cmd = cfg["castep_command"]

    workdir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(exist_ok=True)
    for f in conv_dir.iterdir():  # pseudopotential files referenced by the .cell
        if f.suffix.lower() in _PSEUDO_SUFFIXES and not (workdir / f.name).exists():
            shutil.copy2(f, workdir / f.name)
            print(f"  Copied pseudopotential: {f.name}")

    latest = [plots_dir / f"{seed}_convergence_latest.png"]
    ecut_data = {"x": [], "E": [], "S": []}
    kpn_data = {"x": [], "E": [], "S": [], "labels": [], "spacing": []}
    problems = []
    state = {"failed": False, "stress_missing": False}
    stress_on = {"v": use_stress}

    def replot(suffix):
        make_plot(ecut_data, kpn_data, n_atoms, cfg, stress_on["v"], latest, seed, xc, suffix)

    def job(jobname, kgrid, ecut, reuse):
        """Run (or reuse) one job; returns (energy|None, stress|None, ok)."""
        castep_file = workdir / f"{jobname}.castep"
        if is_complete(castep_file):
            print(f"  [SKIP ] {jobname}  (already complete)")
        elif plot_only:
            print(f"  [SKIP ] {jobname}  (--plot-only)")
            return None, None, True
        else:
            write_cell(cell_tmpl, workdir / f"{jobname}.cell", kgrid)
            write_param(param_tmpl, workdir / f"{jobname}.param", ecut, reuse=reuse, stress=stress_on["v"])
            if not run_castep(jobname, workdir, cmd):
                problems.append(f"Job failed: {jobname} (see {jobname}.castep and any *.err file in {workdir})")
                print(f"  [ABORT] Job failed: {jobname}. Stopping.")
                state["failed"] = True
                return None, None, False
        E = parse_energy(castep_file)
        S = parse_stress(castep_file) if stress_on["v"] else None
        if E is not None and stress_on["v"] and S is None and not state["stress_missing"]:
            state["stress_missing"] = True
            stress_on["v"] = False
            problems.append("No stress tensor found in the CASTEP output; the stress criterion was dropped "
                            "(energy only).")
            print("  ! No stress tensor in the output: continuing with the energy criterion only.")
        return E, S, True

    def add(data, E, S, **extra):
        data["E"].append(E)
        data["S"].append(S)
        for k, v in extra.items():
            data[k].append(v)

    def stress_list(data):
        return data["S"] if stress_on["v"] and all(s is not None for s in data["S"]) else []

    def status_line(label, E, dE, dS, conv):
        de = f"{dE[-1]:.4f}" if dE else "---"
        ds = f"  dσ={dS[-1]:.3f} GPa" if dS else ""
        print(f"  {label}  E={E:.6f} eV  dE={de} meV/atom{ds}  {'CONVERGED' if conv else 'x'}")

    # -- Phase 1: cut-off --------------------------------------------------
    ecut_conv, ecut_ok = ecut_fixed, ecut_fixed is not None
    fixed = pl["fixed_grid"]
    if phase in ("both", "ecut"):
        print("\n" + "=" * 66)
        print(f"PHASE 1 - cut-off convergence  (fixed k-grid: {_fmt_grid(fixed)})")
        print(f"Criteria: dE < {thr:g} meV/atom" + (f", dσ < {thr_s:g} GPa" if stress_on['v'] else "")
              + f" for {req} consecutive steps")
        print("=" * 66)
        prev = None
        ecut_conv, ecut_ok = None, False
        for ecut in pl["ecuts"]:
            name = f"{seed}_cut_{float(ecut):.1f}_kpn_{_fmt_grid(fixed)}"
            reuse = f"{prev}.check" if (cfg["reuse"] and prev and (workdir / f"{prev}.check").exists()) else None
            E, S, ok = job(name, fixed, ecut, reuse)
            if E is not None:
                add(ecut_data, E, S, x=float(ecut))
                prev = name
            conv, dE, dS = check_convergence(ecut_data["E"], stress_list(ecut_data), thr, n_atoms, req,
                                             thr_s if stress_on["v"] else None)
            if E is not None:
                status_line(f"Ecut={ecut:6.1f} eV", E, dE, dS, conv)
            replot(" [Phase 1 - in progress]")
            if not ok:
                break
            if conv:
                ecut_conv, ecut_ok = ecut, True
                print(f"\n  Cut-off CONVERGED at {ecut_conv} eV")
                break
        if state["failed"]:
            ecut_conv = int(ecut_data["x"][-1]) if ecut_data["x"] else None
        elif not ecut_ok:
            ecut_conv = int(ecut_data["x"][-1]) if ecut_data["x"] else cfg["ecut_start"]
            print(f"\n  Cut-off NOT converged within {cfg['ecut_max']} eV (last tested: {ecut_conv} eV). "
                  "Increase ecut_max in the .conv file to extend.")
            problems.append(f"Cut-off not converged up to {ecut_conv} eV")

    # -- Phase 2: k-points -------------------------------------------------
    kconv_grid, kconv_n, kok = None, None, False
    if phase in ("both", "kpoints") and not state["failed"]:
        if ecut_conv is None:
            raise ValueError("a k-point-only run needs the cut-off: use --ecut")
        print("\n" + "=" * 66)
        print(f"PHASE 2 - k-point convergence  (fixed cut-off: {ecut_conv:g} eV)")
        print(f"Criteria: dE < {thr:g} meV/atom" + (f", dσ < {thr_s:g} GPa" if stress_on['v'] else "")
              + f" for {req} consecutive steps")
        print("=" * 66)
        for n, grid, spc in pl["k_steps"]:
            name = f"{seed}_cut_{float(ecut_conv):.1f}_kpn_{_fmt_grid(grid)}"
            E, S, ok = job(name, grid, ecut_conv, None)  # different meshes: no REUSE
            if E is not None:
                add(kpn_data, E, S, x=n, labels=_fmt_grid(grid).replace("x", "×"), spacing=spc)
            conv, dE, dS = check_convergence(kpn_data["E"], stress_list(kpn_data), thr, n_atoms, req,
                                             thr_s if stress_on["v"] else None)
            if E is not None:
                status_line(f"k-grid={_fmt_grid(grid):>9} (Δ={spc:.3f})", E, dE, dS, conv)
            replot(" [Phase 2 - in progress]")
            if not ok:
                break
            if conv:
                kconv_grid, kconv_n, kok = grid, n, True
                print(f"\n  k-mesh CONVERGED at {_fmt_grid(grid)}")
                break
        if not kok and not state["failed"]:
            if kpn_data["x"]:
                kconv_n = kpn_data["x"][-1]
                kconv_grid = next(g for nn, g, _ in pl["k_steps"] if nn == kconv_n)
            last = pl["k_steps"][-1]
            print(f"\n  k-mesh NOT converged up to {_fmt_grid(last[1])} (spacing {last[2]:.3f} 1/A). "
                  "Increase kpn_max (or lower kspacing_min) to extend.")
            problems.append(f"k-points not converged up to {_fmt_grid(kconv_grid or last[1])}")
        elif state["failed"] and kpn_data["x"]:
            kconv_n = kpn_data["x"][-1]
            kconv_grid = next(g for nn, g, _ in pl["k_steps"] if nn == kconv_n)

    # -- Final check one cut-off step higher --------------------------------
    verified = None
    if cfg["final_check"] and phase == "both" and ecut_ok and kok and not state["failed"]:
        e_hi = int(ecut_conv + cfg["ecut_step"])
        print("\n" + "=" * 66)
        print(f"FINAL CHECK - cut-off {e_hi} eV at the converged k-mesh {_fmt_grid(kconv_grid)}")
        print("=" * 66)
        name_ref = f"{seed}_cut_{float(ecut_conv):.1f}_kpn_{_fmt_grid(kconv_grid)}"
        name_hi = f"{seed}_cut_{float(e_hi):.1f}_kpn_{_fmt_grid(kconv_grid)}"
        f_ref = workdir / f"{name_ref}.castep"
        E_ref = parse_energy(f_ref)
        S_ref = parse_stress(f_ref) if stress_on["v"] else None
        E_hi, S_hi, ok = job(name_hi, kconv_grid, e_hi, None)
        if E_hi is not None and E_ref is not None:
            d_e = abs(E_hi - E_ref) * 1000 / n_atoms
            d_s = (float(np.max(np.abs(np.array(S_hi) - np.array(S_ref))))
                   if (stress_on["v"] and S_hi is not None and S_ref is not None) else None)
            verified = d_e < thr and (d_s is None or d_s < thr_s)
            print(f"  dE = {d_e:.4f} meV/atom" + (f", dσ = {d_s:.3f} GPa" if d_s is not None else "")
                  + f"  ->  {'VERIFIED' if verified else 'STILL CHANGING'}")
            if not verified:
                problems.append(f"Final check failed: raising the cut-off to {e_hi} eV at the converged k-mesh "
                                f"still changes the result (dE = {d_e:.3f} meV/atom). Re-run with a larger ecut_start/"
                                "ecut_step or a tighter require_n_conv.")

    # -- Outputs -------------------------------------------------------------
    rec = []
    if ecut_conv is not None and ecut_data["x"]:
        rec.append(f"  cut_off_energy      : {ecut_conv:g} eV" + ("" if ecut_ok else "   (NOT converged - last tested)"))
    k_out = None
    if kconv_grid is not None:
        sp = kmesh.spacing_for_grid(lat, kconv_grid)
        sep = kmesh.per_axis_spacing(lat, kconv_grid)
        k_out = {"grid": kconv_grid, "spacing": sp, "separation": sep}
        flag = "" if kok else "   (NOT converged - last tested)"
        rec += [f"  KPOINTS_MP_GRID     : {kconv_grid[0]} {kconv_grid[1]} {kconv_grid[2]}{flag}",
                f"  KPOINTS_MP_SPACING  : {sp:g} 1/ang   (CASTEP convention, no 2pi; reproduces the grid on this cell)",
                "  k-point separation  : " + "  ".join(f"{x:.4f}" for x in sep) + " 1/A along b1 b2 b3"]
    suffix = ""
    if rec:
        parts = ([f"Ecut={ecut_conv:g} eV"] if ecut_conv is not None and ecut_data["x"] else [])
        if kconv_grid is not None:
            parts.append(f"k-grid={_fmt_grid(kconv_grid)} (Δ={k_out['spacing']:.3f} 1/Å)")
        suffix = "\n" + ", ".join(parts)
    final = [plots_dir / f"{seed}_convergence_final.png", plots_dir / f"{seed}_convergence_final.pdf"]
    make_plot(ecut_data, kpn_data, n_atoms, cfg, stress_on["v"], final + latest, seed, xc, suffix)

    crit = f"dE < {thr:g} meV/atom" + (f" and d(stress) < {thr_s:g} GPa" if stress_on["v"] else "")
    rep = [f"Adaptive convergence report - {seed}",
           f"Generated : {datetime.now():%Y-%m-%d %H:%M:%S}",
           f"Criteria  : {crit}, {req} consecutive steps | {n_atoms} atoms", "=" * 66]
    if ecut_data["x"]:
        rep.append(f"\nPhase 1 - cut-off (fixed k-grid: {_fmt_grid(fixed)})")
        for i, (x, E) in enumerate(zip(ecut_data["x"], ecut_data["E"])):
            d = f"{abs(E - ecut_data['E'][i - 1]) * 1000 / n_atoms:.4f}" if i else "---"
            rep.append(f"  Ecut={x:7.1f} eV  E={E:>18.9f} eV  dE={d:>10} meV/atom")
    if kpn_data["x"]:
        rep.append(f"\nPhase 2 - k-points (fixed cut-off: {ecut_conv:g} eV)")
        for i, (lbl, E, spc) in enumerate(zip(kpn_data["labels"], kpn_data["E"], kpn_data["spacing"])):
            d = f"{abs(E - kpn_data['E'][i - 1]) * 1000 / n_atoms:.4f}" if i else "---"
            rep.append(f"  k-grid={lbl:>9}  spacing={spc:6.4f}  E={E:>18.9f} eV  dE={d:>10} meV/atom")
    if verified is not None:
        rep.append(f"\nFinal check (cut-off + {cfg['ecut_step']} eV at the converged k-mesh): "
                   + ("verified" if verified else "STILL CHANGING"))
    rep += ["", "=" * 66, "RECOMMENDED PRODUCTION PARAMETERS"] + (rec or ["  (no data yet)"]) + ["=" * 66]
    if problems:
        rep += ["", "Notes:"] + [f"  - {p}" for p in problems]
    report_text = "\n".join(rep)
    report_path = conv_dir / f"{seed}_convergence_report.txt"
    report_path.write_text(report_text + "\n", encoding="utf-8")
    return {"ecut": ecut_conv, "ecut_converged": ecut_ok, "kgrid": kconv_grid, "k": k_out,
            "kpn_converged": kok, "verified": verified, "report": report_text,
            "report_path": report_path, "plot": final[0], "problems": problems}
