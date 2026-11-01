# -*- coding: utf-8 -*-
"""
Geometry-optimisation convergence analysis for CASTEP ``.castep`` output files.

Pure Python (re + numpy + matplotlib): no grep/awk, so it works on Windows, and it
is not tied to one optimiser (LBFGS, BFGS, FIRE, DELOC, TPSD ...).

Inspired by ``check_converge.py`` from https://github.com/tjz21/Castep_Tools
"""
import re
from pathlib import Path

import numpy as np

# "<METHOD>: finished iteration  N with enthalpy= -1.10666689E+004 eV"
_ITER_RE = re.compile(
    r'^\s*([A-Za-z]+):\s*finished iteration\s+(\d+)\s+with enthalpy=\s*([-+0-9.EeDd]+)\s*eV',
    re.MULTILINE,
)
# "| |F|max | 4.157448E-001 | 5.000000E-002 | eV/A | No  |"
_ROW_RE = re.compile(
    r'^\s*\|\s*(.+?)\s*\|\s*([-+0-9.EeDd]+)\s*\|\s*([-+0-9.EeDd]+)\s*\|\s*(\S+)\s*\|\s*(Yes|No)\s*\|',
    re.MULTILINE,
)
_VOL_RE = re.compile(r'Current cell volume\s*=\s*([0-9.]+)\s*A\*\*3')
_SUCCESS_RE = re.compile(r'Geometry optimi[sz]ation completed successfully', re.IGNORECASE)

# Display metadata for the criteria CASTEP prints, keyed by the name in the table.
CRITERIA = {
    "dE/ion": ("Energy change", r"$\Delta E$ per ion (eV)", "tab:blue"),
    "|F|max": ("Max force", r"$|F|_{\max}$ (eV/Å)", "tab:red"),
    "|dR|max": ("Max displacement", r"$|dR|_{\max}$ (Å)", "tab:green"),
    "Smax": ("Max stress", r"$S_{\max}$ (GPa)", "tab:purple"),
}


def _f(s):
    return float(s.replace("D", "E").replace("d", "e"))


def parse_geomopt(text):
    """
    Parse a geometry-optimisation run from the text of a ``.castep`` file.

    Returns a dict with: ``method``, ``converged`` (bool, CASTEP's own success
    message), ``steps`` (list of dicts: ``n``, ``enthalpy``, ``volume`` and
    ``criteria`` = {name: (value, tolerance, units, ok)}), and ``criteria_names``.
    Raises ValueError if no finished optimisation iteration is found.
    """
    matches = list(_ITER_RE.finditer(text))
    if not matches:
        raise ValueError("no 'finished iteration' lines found - is this a "
                         "GeometryOptimization .castep file?")

    steps = []
    method = matches[0].group(1)
    for k, m in enumerate(matches):
        n = int(m.group(2))
        # A repeated/lower iteration number means the job was restarted from
        # scratch inside the same file: keep only the latest attempt.
        if steps and n <= steps[-1]["n"]:
            steps = []
        end = matches[k + 1].start() if k + 1 < len(matches) else len(text)
        # The convergence table directly follows the "finished iteration" line.
        block = text[m.end():end]
        crit = {}
        for r in _ROW_RE.finditer(block):
            name = r.group(1)
            if name in crit:  # table of the next section, stop at first repeat
                break
            crit[name] = (_f(r.group(2)), _f(r.group(3)), r.group(4), r.group(5) == "Yes")
        # Accepted geometry of this iteration = last cell volume printed before
        # the "finished iteration" line (earlier ones are line-search trial steps).
        vols = _VOL_RE.findall(text[:m.start()])
        steps.append({
            "n": n,
            "enthalpy": _f(m.group(3)),
            "volume": float(vols[-1]) if vols else float("nan"),
            "criteria": crit,
        })
        method = m.group(1)

    names = []
    for s in steps:
        for name in s["criteria"]:
            if name not in names:
                names.append(name)
    return {
        "method": method,
        "converged": bool(_SUCCESS_RE.search(text)),
        "steps": steps,
        "criteria_names": names,
    }


def summary_lines(data, source=""):
    """Plain-text summary: overall status plus the criteria at the last iteration."""
    steps = data["steps"]
    last = steps[-1]
    lines = []
    if source:
        lines.append(f"File            : {source}")
    lines += [
        f"Optimiser       : {data['method']}",
        f"Iterations      : {last['n']}",
        f"Final enthalpy  : {last['enthalpy']:.6f} eV",
    ]
    if not np.isnan(last["volume"]):
        lines.append(f"Final volume    : {last['volume']:.4f} A^3")
    all_ok = bool(last["criteria"]) and all(c[3] for c in last["criteria"].values())
    if data["converged"]:
        status = "CONVERGED (CASTEP reported successful completion)"
    elif all_ok:
        status = "all criteria met at the last iteration, but no completion message (run cut off?)"
    else:
        status = "NOT CONVERGED"
    lines += [f"Status          : {status}", "",
              f"{'Criterion':<10} {'Final value':>14} {'Tolerance':>14}  {'Units':<5} Converged"]
    failing = []
    for name in data["criteria_names"]:
        if name not in last["criteria"]:
            continue
        v, tol, units, ok = last["criteria"][name]
        lines.append(f"{name:<10} {v:>14.6e} {tol:>14.6e}  {units:<5} {'Yes' if ok else 'No'}")
        if not ok:
            failing.append(f"{name} ({v:.3g} vs tol {tol:.3g}, x{v / tol:.1f})" if tol else name)
    if failing and not data["converged"]:
        lines += ["", "Not yet converged: " + ", ".join(failing)]
    return lines


def plot_convergence(data, title=None):
    """Return a matplotlib Figure: enthalpy, volume (if it changes) and each criterion."""
    import matplotlib
    import matplotlib.pyplot as plt

    steps = data["steps"]
    n = np.array([s["n"] for s in steps])
    enth = np.array([s["enthalpy"] for s in steps])
    vol = np.array([s["volume"] for s in steps])

    panels = [("Enthalpy", "Enthalpy (eV)", "tab:orange", n, enth, None, None)]
    if np.all(np.isfinite(vol)) and np.ptp(vol) > 1e-6 * abs(vol[0]):
        panels.append(("Cell volume", r"Volume (Å$^3$)", "tab:brown", n, vol, None, None))

    for name in data["criteria_names"]:
        ttl, ylabel, color = CRITERIA.get(name, (name, name, "tab:gray"))
        pts = [(s["n"], s["criteria"][name]) for s in steps if name in s["criteria"]]
        x = np.array([p[0] for p in pts])
        y = np.array([p[1][0] for p in pts])
        ok = np.array([p[1][3] for p in pts])
        tol = pts[-1][1][1]
        panels.append((ttl, ylabel, color, x, y, tol, ok))

    ncols = 2
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(11, 3.3 * nrows + 0.8), squeeze=False)
    axs = axs.ravel()

    for ax, (ttl, ylabel, color, x, y, tol, ok) in zip(axs, panels):
        if tol is None:
            ax.plot(x, y, "o-", color=color)
            ax.ticklabel_format(useOffset=False, style="plain", axis="y")
            ax.annotate(f"Final = {y[-1]:.4f}", xy=(0.97, 0.95), xycoords="axes fraction",
                        ha="right", va="top", fontsize=9)
        else:
            # Log axis: CASTEP prints exactly 0 for dE/ion and dR at iteration 0.
            pos = y > 0
            ax.plot(x[pos], y[pos], "-", color=color, lw=1, zorder=1)
            ax.scatter(x[pos & ok], y[pos & ok], s=36, color="tab:green", zorder=2,
                       label="criterion met")
            ax.scatter(x[pos & ~ok], y[pos & ~ok], s=36, facecolors="white",
                       edgecolors=color, zorder=2, label="not met")
            ax.axhline(tol, color="k", ls="--", lw=1, label=f"tol = {tol:.2e}")
            ax.set_yscale("log")
            ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
            state = "converged" if ok[-1] else "NOT converged"
            ax.annotate(f"Final = {y[-1]:.2e}  ({state})", xy=(0.03, 0.05),
                        xycoords="axes fraction", ha="left", va="bottom", fontsize=9,
                        color="darkgreen" if ok[-1] else "firebrick")
        ax.set_title(ttl, fontsize=11)
        ax.set_xlabel("Optimisation step")
        ax.set_ylabel(ylabel)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.grid(alpha=0.3)
    for ax in axs[len(panels):]:
        ax.set_visible(False)

    status = "converged" if data["converged"] else "NOT converged"
    fig.suptitle(f"{title or 'CASTEP geometry optimisation'} - {data['method']}, "
                 f"{steps[-1]['n']} steps, {status}",
                 fontsize=13, color="darkgreen" if data["converged"] else "firebrick")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig
