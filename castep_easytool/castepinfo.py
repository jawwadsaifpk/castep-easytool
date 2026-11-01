# -*- coding: utf-8 -*-
"""
Run summary and SCF-convergence analysis for CASTEP ``.castep`` output files.
Pure Python (re + numpy + matplotlib).
"""
import re

import numpy as np

_NUM = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'


def _f(s):
    return float(s.replace("D", "E").replace("d", "e"))


def _last(pattern, text, flags=re.MULTILINE):
    m = list(re.finditer(pattern, text, flags))
    return m[-1] if m else None


def _first(pattern, text, flags=re.MULTILINE):
    return re.search(pattern, text, flags)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def parse_summary(text):
    """Collect the key inputs and results of a CASTEP run into a dict (None = not found)."""
    # Parameters come from the latest run header if the file holds several runs.
    start = text.rfind("Run started")
    seg = text[start:] if start >= 0 else text

    def val(pattern, src=seg, cast=str, group=1):
        m = _last(pattern, src)
        if not m:
            return None
        try:
            return cast(m.group(group).strip())
        except ValueError:
            return None

    d = {}
    d["run_started"] = val(r'Run started:\s*(.+)')
    d["task"] = val(r'^\s*type of calculation\s*:\s*(.+)$')
    d["functional"] = val(r'^\s*using functional\s*:\s*(.+)$')
    d["dispersion"] = val(r'DFT\+D: Semi-empirical dispersion correction\s*:\s*(.+)$')
    d["cutoff_eV"] = val(rf'plane wave basis set cut-off\s*:\s*({_NUM})\s*eV', cast=_f)
    spin = _last(r'treating system as (non-spin-polarized|spin-polarized)', seg)
    d["spin_polarized"] = (spin.group(1) == "spin-polarized") if spin else None
    d["n_electrons"] = val(rf'number of\s+electrons\s*:\s*({_NUM})', cast=_f)
    d["n_bands"] = val(r'number of bands\s*:\s*(\d+)', cast=int)
    d["smearing"] = val(r'smearing scheme\s*:\s*(.+)$')
    d["smearing_width_eV"] = val(rf'smearing width\s*:\s*({_NUM})\s*eV', cast=_f)
    d["scf_tol_eV_atom"] = val(rf'total energy / atom convergence tol\.\s*:\s*({_NUM})', cast=_f)
    d["max_scf_cycles"] = val(r'max\. number of SCF cycles\s*:\s*(\d+)', cast=int)
    d["optimizer"] = val(r'optimization method\s*:\s*(.+)$')
    d["n_ions"] = val(r'Total number of ions in cell\s*=\s*(\d+)', cast=int)
    d["n_species"] = val(r'Total number of species in cell\s*=\s*(\d+)', cast=int)
    d["space_group"] = val(r'Space group of crystal\s*=\s*(.+)$')
    mp = _last(r'MP grid size for SCF calculation is\s+(\d+)\s+(\d+)\s+(\d+)', seg)
    d["kgrid"] = tuple(int(x) for x in mp.groups()) if mp else None
    d["n_kpoints"] = val(r'Number of kpoints used\s*=\s*(\d+)', cast=int)

    vols = re.findall(r'Current cell volume\s*=\s*([0-9.]+)\s*A\*\*3', text)
    d["volume_initial"] = float(vols[0]) if vols else None
    d["volume_final"] = float(vols[-1]) if vols else None

    d["final_energy_eV"] = val(rf'Final energy, E\s*=\s*({_NUM})\s*eV', text, _f)
    d["final_free_energy_eV"] = val(rf'Final free energy \(E-TS\)\s*=\s*({_NUM})\s*eV', text, _f)
    d["final_enthalpy_eV"] = val(rf'Final Enthalpy\s*=\s*({_NUM})\s*eV', text, _f)
    d["pressure_GPa"] = val(rf'\*\s*Pressure:\s*({_NUM})', text, _f)
    d["total_time_s"] = val(rf'^Total time\s*=\s*({_NUM})\s*s', text, _f)
    d["peak_memory_kB"] = val(r'Peak Memory Use\s*=\s*(\d+)\s*kB', text, int)
    d["n_scf_loops"] = len(re.findall(r'^\s*Initial\s+' + _NUM + r'.*<-- SCF', text, re.MULTILINE))
    d["geomopt_converged"] = (
        True if re.search(r'Geometry optimi[sz]ation completed successfully', text, re.I) else
        (False if re.search(r'finished iteration', text) else None))
    d["finished"] = d["total_time_s"] is not None

    # Warnings: keep unique texts, first occurrences first.
    warns = []
    for m in re.finditer(r'^\s*(?:\S+\s+)?[Ww]arning[: ].*$', text, re.MULTILINE):
        w = m.group(0).strip()
        if w not in warns:
            warns.append(w)
    d["warnings"] = warns
    return d


def summary_lines(d, source=""):
    def fmt(label, value, unit=""):
        return None if value is None else f"{label:<22}: {value}{(' ' + unit) if unit else ''}"

    spin = None if d["spin_polarized"] is None else ("spin-polarised" if d["spin_polarized"] else "non-spin-polarised")
    kgrid = None
    if d["kgrid"]:
        kgrid = "x".join(map(str, d["kgrid"])) + (f"  ({d['n_kpoints']} irreducible k-points)" if d["n_kpoints"] else "")
    elif d["n_kpoints"]:
        kgrid = f"{d['n_kpoints']} k-points"
    smear = None
    if d["smearing"]:
        smear = d["smearing"] + (f", {d['smearing_width_eV']:g} eV" if d["smearing_width_eV"] is not None else "")
    vol = None
    if d["volume_initial"] is not None:
        vol = f"{d['volume_initial']:.3f}"
        if d["volume_final"] is not None and abs(d["volume_final"] - d["volume_initial"]) > 1e-6:
            vol += f" -> {d['volume_final']:.3f}  ({100 * (d['volume_final'] / d['volume_initial'] - 1):+.2f} %)"
    t = d["total_time_s"]
    ttxt = None if t is None else (f"{t:.1f} s" + (f"  ({t / 3600:.2f} h)" if t > 600 else ""))

    rows = [
        "Calculation", fmt("  Task", d["task"]), fmt("  Optimiser", d["optimizer"]),
        fmt("  Started", d["run_started"]),
        "Inputs", fmt("  XC functional", d["functional"]), fmt("  Dispersion", d["dispersion"]),
        fmt("  Cut-off energy", d["cutoff_eV"], "eV"), fmt("  k-point grid", kgrid),
        fmt("  Spin", spin), fmt("  Smearing", smear),
        fmt("  SCF tolerance", d["scf_tol_eV_atom"], "eV/atom"),
        "System", fmt("  Ions / species", None if d["n_ions"] is None else f"{d['n_ions']} / {d['n_species']}"),
        fmt("  Electrons / bands", None if d["n_electrons"] is None else f"{d['n_electrons']:g} / {d['n_bands']}"),
        fmt("  Space group", d["space_group"]), fmt("  Cell volume (A^3)", vol),
        "Results", fmt("  Final enthalpy", None if d["final_enthalpy_eV"] is None else f"{d['final_enthalpy_eV']:.6f}", "eV"),
        fmt("  Final energy E", None if d["final_energy_eV"] is None else f"{d['final_energy_eV']:.6f}", "eV"),
        fmt("  Free energy E-TS", None if d["final_free_energy_eV"] is None else f"{d['final_free_energy_eV']:.6f}", "eV"),
        fmt("  Energy per ion", None if (d["final_free_energy_eV"] is None or not d["n_ions"])
            else f"{d['final_free_energy_eV'] / d['n_ions']:.6f}", "eV"),
        fmt("  Pressure", d["pressure_GPa"], "GPa"),
        fmt("  SCF loops", d["n_scf_loops"] or None),
        "Run", fmt("  Total time", ttxt),
        fmt("  Peak memory", None if d["peak_memory_kB"] is None else f"{d['peak_memory_kB'] / 1e6:.2f}", "GB"),
    ]
    if d["geomopt_converged"] is not None:
        rows.append(fmt("  Geometry optimisation", "converged" if d["geomopt_converged"] else "NOT converged"))
    rows.append(fmt("  Status", "finished" if d["finished"] else "unfinished (no timing summary: still running, killed or crashed)"))
    out = [f"CASTEP run summary" + (f" - {source}" if source else ""), "=" * 60]
    section_pending = None
    for r in rows:
        if r is None:
            continue
        if not r.startswith(" "):
            section_pending = r
            continue
        if section_pending:
            out += ["", section_pending + ":"]
            section_pending = None
        out.append(r)
    if d["warnings"]:
        out += ["", f"Warnings ({len(d['warnings'])} unique):"] + [f"  {w}" for w in d["warnings"][:10]]
        if len(d["warnings"]) > 10:
            out.append(f"  ... and {len(d['warnings']) - 10} more")
    else:
        out += ["", "Warnings: none"]
    return out


# ---------------------------------------------------------------------------
# SCF convergence
# ---------------------------------------------------------------------------

def parse_scf(text):
    """
    Return a list of SCF loops. Each loop is a list of dicts with ``cycle``
    (0 = 'Initial'), ``energy`` (eV), ``gain`` (eV/atom, None for Initial) and ``time`` (s).
    """
    loops, cur = [], None
    for line in text.splitlines():
        if "<-- SCF" not in line:
            continue
        tok = line.split("<--")[0].split()
        if len(tok) < 4:
            continue
        if tok[0] == "Initial":
            cur = []
            loops.append(cur)
            cycle = 0
        elif tok[0].isdigit():
            cycle = int(tok[0])
            if cur is None:
                cur = []
                loops.append(cur)
        else:
            continue
        try:
            energy = _f(tok[1])
            time = _f(tok[-1])
            gain = _f(tok[-2]) if cycle > 0 and len(tok) >= 5 else None
        except ValueError:
            continue
        cur.append({"cycle": cycle, "energy": energy, "gain": gain, "time": time})
    loops = [lp for lp in loops if lp]
    if not loops:
        raise ValueError("no SCF cycles found - is this a CASTEP .castep output file?")
    return loops


def scf_lines(loops, tol=None, max_cycles=None, source=""):
    n = [len(lp) - 1 for lp in loops]  # cycles excluding 'Initial'
    last = loops[-1]
    out = []
    if source:
        out.append(f"File            : {source}")
    out += [f"SCF loops       : {len(loops)}",
            f"Cycles per loop : min {min(n)}, max {max(n)}, mean {np.mean(n):.1f}",
            f"Last loop       : {n[-1]} cycles, final energy {last[-1]['energy']:.8f} eV"]
    if tol is not None:
        out.append(f"Tolerance       : {tol:g} eV/atom (CASTEP requires it to be met over its convergence window)")
    if max_cycles:
        hit = [i + 1 for i, k in enumerate(n) if k >= max_cycles]
        if hit:
            out.append(f"WARNING         : loop(s) {hit} reached the maximum of {max_cycles} SCF cycles - "
                       "likely not converged (try a smaller mixing amplitude, more SCF cycles, or check the structure)")
        else:
            out.append(f"Status          : all loops converged within the {max_cycles}-cycle limit")
    slow = [(i + 1, k) for i, k in enumerate(n) if k > 2 * max(np.median(n), 5)]
    if slow:
        out.append("Slow loops      : " + ", ".join(f"#{i} ({k} cycles)" for i, k in slow[:8]))
    return out


def plot_scf(loops, loop_index=-1, tol=None, title=None):
    """Two panels: |energy gain| per cycle for one loop, and cycles needed by every loop."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    lp = loops[loop_index]
    idx = (loop_index % len(loops)) + 1
    x = np.array([r["cycle"] for r in lp if r["gain"] is not None])
    y = np.abs(np.array([r["gain"] for r in lp if r["gain"] is not None]))
    pos = y > 0
    a1.plot(x[pos], y[pos], "o-", color="tab:blue")
    if tol:
        a1.axhline(tol, color="k", ls="--", lw=1, label=f"tol = {tol:.1e} eV/atom")
        a1.legend(fontsize=8)
    a1.set_yscale("log")
    a1.set_xlabel("SCF cycle")
    a1.set_ylabel("|energy gain| (eV/atom)")
    a1.set_title(f"SCF loop #{idx} of {len(loops)}")
    a1.xaxis.set_major_locator(MaxNLocator(integer=True))
    a1.grid(alpha=0.3)

    n = [len(l) - 1 for l in loops]
    a2.bar(np.arange(1, len(n) + 1), n, color="tab:orange")
    a2.set_xlabel("SCF loop (one per ionic step / trial step)")
    a2.set_ylabel("SCF cycles")
    a2.set_title("Cycles per SCF loop")
    a2.xaxis.set_major_locator(MaxNLocator(integer=True))
    a2.grid(alpha=0.3, axis="y")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
