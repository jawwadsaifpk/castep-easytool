# -*- coding: utf-8 -*-
"""
Self-contained reader, band-gap analysis and plotter for CASTEP band structures.

This module replaces the previous dependency on ``sumo``. It only needs
``numpy`` and ``matplotlib``.

Inputs
------
* ``<seed>.bands`` -- CASTEP eigenvalue file (lattice, k-points, eigenvalues,
  electron counts, Fermi energy).
* ``<seed>.cell``  -- optional source of high-symmetry k-point labels, read from
  the ``[BS|SPECTRAL]_KPOINT[S]_[PATH|LIST]`` block (labels are taken from the
  trailing ``! label`` comment on each line, as written by
  ``castep-easytool bands``).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903

# Matching tolerance (fractional coordinates) between .bands k-points and
# vertices of the path block in the .cell file.
_KPT_TOL = 1e-4
# Gaps smaller than this (eV) are reported as metallic / zero-gap.
_GAP_TOL = 1e-4


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class BandStructure:
    lattice: np.ndarray                       # (3, 3) Angstrom, rows = vectors
    kpoints: np.ndarray                       # (nk, 3) fractional
    eigenvalues: Dict[int, np.ndarray]        # spin -> (nbands, nk) in eV
    efermi: List[float]                       # eV, one per spin channel (or one)
    n_electrons: List[float]                  # one per spin channel (or one)
    labels: List[Optional[str]] = field(default_factory=list)  # per k-point
    # (i, j): the path jumps from k-point i to k-point j; k-points strictly
    # between them only exist because CASTEP joins consecutive vertices.
    breaks: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def n_spins(self) -> int:
        return len(self.eigenvalues)

    @property
    def is_spin_polarized(self) -> bool:
        return self.n_spins == 2

    @property
    def valid(self) -> np.ndarray:
        """Boolean mask of k-points that lie on a real part of the path."""
        mask = np.ones(len(self.kpoints), dtype=bool)
        for i, j in self.breaks:
            mask[i + 1:j] = False
        return mask

    @property
    def reciprocal(self) -> np.ndarray:
        """Reciprocal lattice (rows, 1/Angstrom, including the 2*pi factor)."""
        return 2.0 * np.pi * np.linalg.inv(self.lattice).T


@dataclass
class ChannelGap:
    """Band-edge analysis of one spin channel."""
    metallic: bool
    gap: Optional[float]      # eV; None if it cannot be determined
    direct: bool
    vbm: Optional[float]      # eV
    cbm: Optional[float]      # eV
    vbm_k: Optional[int]      # k-point index
    cbm_k: Optional[int]
    n_occ: Optional[int]      # number of fully occupied bands (None for metals)


# ---------------------------------------------------------------------------
# .bands reader
# ---------------------------------------------------------------------------

def _floats(tokens) -> List[float]:
    return [float(t.replace("D", "E").replace("d", "e")) for t in tokens]


def read_bands(bands_file) -> BandStructure:
    """Parse a CASTEP ``.bands`` file. Energies are converted from Hartree to eV."""
    lines = [ln.strip() for ln in Path(bands_file).read_text(encoding="utf-8", errors="replace").splitlines()]
    lines = [ln for ln in lines if ln]

    n_k = n_spins = None
    n_el: List[float] = []
    n_bands: List[int] = []
    efermi: List[float] = []
    lattice = None
    body_start = None

    for i, ln in enumerate(lines):
        low = ln.lower()
        if low.startswith("number of k-points"):
            n_k = int(ln.split()[-1])
        elif low.startswith("number of spin components"):
            n_spins = int(ln.split()[-1])
        elif low.startswith("number of electrons"):
            n_el = _floats(ln.split()[3:])
        elif low.startswith("number of eigenvalues"):
            n_bands = [int(x) for x in ln.split()[3:]]
        elif low.startswith("fermi energ"):
            if "atomic" not in low:
                raise ValueError("Band data are not in atomic units; unsupported .bands file.")
            efermi = _floats(ln.split(")")[-1].split())
        elif low.startswith("unit cell vectors"):
            lattice = np.array([_floats(lines[i + 1 + j].split()[:3]) for j in range(3)])
            body_start = i + 4
            break

    if n_k is None or n_spins is None or lattice is None or body_start is None or not n_el or not n_bands or not efermi:
        raise ValueError("Could not read the header of the .bands file (unexpected format).")
    if n_spins not in (1, 2):
        raise ValueError(f"Unsupported number of spin components: {n_spins}")
    if len(n_bands) < n_spins:
        n_bands = n_bands * n_spins

    # Body: "K-point idx kx ky kz weight", then "Spin component s" + eigenvalues
    kpts: List[Tuple[int, np.ndarray]] = []
    eig: Dict[int, List[List[float]]] = {s: [] for s in range(n_spins)}
    current: Optional[List[float]] = None
    spin = 0
    for ln in lines[body_start:]:
        low = ln.lower()
        if low.startswith("k-point"):
            vals = _floats(ln.split()[1:])
            kpts.append((int(vals[0]), np.array(vals[1:4])))
            current = None
        elif low.startswith("spin component"):
            spin = int(ln.split()[-1]) - 1
            current = []
            eig[spin].append(current)
        elif current is not None:
            current.append(_floats(ln.split())[0])

    if len(kpts) != n_k:
        raise ValueError(f"Expected {n_k} k-points but found {len(kpts)}.")
    for s in range(n_spins):
        if len(eig[s]) != n_k or any(len(e) != n_bands[s] for e in eig[s]):
            raise ValueError("Eigenvalue block sizes do not match the header.")

    order = np.argsort([k[0] for k in kpts], kind="stable")
    kfrac = np.array([kpts[i][1] for i in order])
    eigenvalues = {
        s: np.array(eig[s])[order].T * HARTREE_TO_EV for s in range(n_spins)
    }

    return BandStructure(
        lattice=lattice * BOHR_TO_ANGSTROM,
        kpoints=kfrac,
        eigenvalues=eigenvalues,
        efermi=[e * HARTREE_TO_EV for e in efermi],
        n_electrons=n_el,
        labels=[None] * n_k,
    )


# ---------------------------------------------------------------------------
# .cell reader (k-point labels only)
# ---------------------------------------------------------------------------

_BLOCK_RE = re.compile(r"^%block\s+(bs|spectral)_kpoints?_(path|list)\b", re.IGNORECASE)
_ENDBLOCK_RE = re.compile(r"^%endblock\s+(bs|spectral)_kpoints?_(path|list)\b", re.IGNORECASE)
_UNIT_RE = re.compile(r"^(bs|spectral)_kpoints?_(path|list)_unit\b", re.IGNORECASE)


def read_cell_kpath(cell_file) -> Tuple[List[Tuple[np.ndarray, Optional[str], bool]], bool]:
    """
    Read the band-path vertices from a CASTEP ``.cell`` file.

    Returns ``([(frac_coords, label_or_None, break_before), ...], non_fractional_units)``.
    The label is the first word of the trailing ``! comment`` (or ``#``/``;``).
    A line containing ``BREAK`` (written by ``castep-easytool bands`` where the
    seekpath path is discontinuous) marks the next vertex as ``break_before``. If both a PATH and a LIST block
    exist, the PATH block is used.
    """
    blocks: Dict[str, List[Tuple[np.ndarray, Optional[str], bool]]] = {}
    non_fractional = False
    current = None
    pending_break = False
    for raw in Path(cell_file).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        m = _BLOCK_RE.match(line)
        if m:
            current = m.group(2).lower()
            blocks.setdefault(current, [])
            continue
        if _ENDBLOCK_RE.match(line):
            current = None
            continue
        if current is None:
            if _UNIT_RE.match(line):
                non_fractional = True
            continue
        if "break" in line.lower():
            pending_break = True
            continue
        if line[0] in "#;!":
            continue
        parts = re.split(r"[#;!]", line, maxsplit=1)
        toks = parts[0].split()
        if len(toks) < 3:
            continue
        try:
            coords = np.array(_floats(toks[:3]))
        except ValueError:
            continue
        label = None
        if len(parts) > 1:
            words = parts[1].split()
            if words:
                label = words[0].strip(",.:;")
        blocks[current].append((coords, label or None, pending_break and bool(blocks[current])))
        pending_break = False

    for key in ("path", "list"):
        if blocks.get(key):
            return blocks[key], non_fractional
    return [], non_fractional


def assign_labels(bs: BandStructure, vertices) -> int:
    """
    Walk along the k-points in order and attach path-vertex labels to the
    k-points that coincide with successive vertices. Returns the number of
    vertices that were matched.
    """
    labels: List[Optional[str]] = [None] * len(bs.kpoints)
    breaks: List[Tuple[int, int]] = []
    ptr = 0
    prev_idx = None
    for i, k in enumerate(bs.kpoints):
        if ptr >= len(vertices):
            break
        coords, label, brk = vertices[ptr]
        if np.allclose(k, coords, atol=_KPT_TOL):
            labels[i] = label
            if brk and prev_idx is not None:
                breaks.append((prev_idx, i))
            prev_idx = i
            ptr += 1
    bs.labels = labels
    bs.breaks = breaks
    return ptr


def kpoint_name(bs: BandStructure, idx: int) -> str:
    """Printable name (unicode Greek) for a k-point, or its coordinates."""
    label = bs.labels[idx] if idx < len(bs.labels) else None
    if label:
        return _plain_label(label)
    k = bs.kpoints[idx]
    return f"({k[0]:.3f},{k[1]:.3f},{k[2]:.3f})"


_GREEK = {
    "G": "Gamma", "GAMMA": "Gamma", "DELTA": "Delta", "LAMBDA": "Lambda",
    "SIGMA": "Sigma", "PI": "Pi", "OMEGA": "Omega", "THETA": "Theta",
}
_GREEK_UNICODE = {
    "Gamma": "Γ", "Delta": "Δ", "Lambda": "Λ", "Sigma": "Σ",
    "Pi": "Π", "Omega": "Ω", "Theta": "Θ",
}


def _split_label(label: str) -> Tuple[str, str]:
    base, _, sub = label.partition("_")
    return base, sub


def _plain_label(label: str) -> str:
    base, sub = _split_label(label)
    name = _GREEK.get(base.upper())
    base = _GREEK_UNICODE[name] if name else base
    return f"{base}_{sub}" if sub else base


def _mathtext_label(label: str) -> str:
    if "|" in label:
        return "|".join(_mathtext_label(p) for p in label.split("|"))
    base, sub = _split_label(label)
    name = _GREEK.get(base.upper())
    core = rf"\{name}" if name else rf"\mathrm{{{base}}}"
    return f"${core}_{{{sub}}}$" if sub else f"${core}$"


# ---------------------------------------------------------------------------
# Band-gap analysis
# ---------------------------------------------------------------------------

def _channel_gap(bs: BandStructure, spin: int) -> ChannelGap:
    E = bs.eigenvalues[spin]                    # (nbands, nk)
    valid = bs.valid
    n_b = E.shape[0]
    n_el = bs.n_electrons[spin if spin < len(bs.n_electrons) else 0]
    occ_f = n_el if bs.is_spin_polarized else n_el / 2.0
    n_occ = int(round(occ_f))

    if abs(occ_f - n_occ) > 1e-3:
        # Fractional occupation per channel -> must be a metal (or magnetic metal).
        return ChannelGap(True, 0.0, False, None, None, None, None, None)

    if not (0 < n_occ < n_b):
        # No empty bands in the file (or no occupied bands): fall back to the
        # Fermi level to decide, and report "undetermined" if there is no gap.
        ef = bs.efermi[spin if spin < len(bs.efermi) else 0]
        Ev = E[:, valid]
        crossing = np.any((Ev.min(axis=1) < ef - 1e-5) & (Ev.max(axis=1) > ef + 1e-5))
        if crossing:
            return ChannelGap(True, 0.0, False, None, None, None, None, None)
        n_occ = int(np.sum(Ev.max(axis=1) <= ef))
        if not (0 < n_occ < n_b):
            return ChannelGap(False, None, False, None, None, None, None, n_occ)

    vk = int(np.argmax(np.where(valid, E[n_occ - 1], -np.inf)))
    ck = int(np.argmin(np.where(valid, E[n_occ], np.inf)))
    vbm, cbm = float(E[n_occ - 1, vk]), float(E[n_occ, ck])
    gap = cbm - vbm
    if gap <= _GAP_TOL:
        return ChannelGap(True, 0.0, False, vbm, cbm, vk, ck, n_occ)
    direct = bool(np.allclose(bs.kpoints[vk], bs.kpoints[ck], atol=_KPT_TOL))
    return ChannelGap(False, gap, direct, vbm, cbm, vk, ck, n_occ)


def analyse_gaps(bs: BandStructure) -> Tuple[Dict[int, ChannelGap], ChannelGap]:
    """Return (per-spin results, overall result)."""
    per_spin = {s: _channel_gap(bs, s) for s in bs.eigenvalues}
    if len(per_spin) == 1:
        return per_spin, per_spin[0]

    chans = list(per_spin.values())
    if any(c.metallic for c in chans):
        return per_spin, ChannelGap(True, 0.0, False, None, None, None, None, None)
    if any(c.gap is None for c in chans):
        return per_spin, ChannelGap(False, None, False, None, None, None, None, None)
    v = max(chans, key=lambda c: c.vbm)
    c = min(chans, key=lambda c: c.cbm)
    gap = c.cbm - v.vbm
    if gap <= _GAP_TOL:
        return per_spin, ChannelGap(True, 0.0, False, v.vbm, c.cbm, v.vbm_k, c.cbm_k, None)
    direct = bool(np.allclose(bs.kpoints[v.vbm_k], bs.kpoints[c.cbm_k], atol=_KPT_TOL))
    return per_spin, ChannelGap(False, gap, direct, v.vbm, c.cbm, v.vbm_k, c.cbm_k, None)


def bandgap_report(bs: BandStructure, system: str) -> List[str]:
    """Text lines for the ``_bandgap.txt`` summary (same layout as before)."""
    per_spin, overall = analyse_gaps(bs)

    lines = [
        "Band Gap Summary",
        "=" * 40,
        f"System          : {system}",
        f"Spin-polarised  : {'Yes' if bs.is_spin_polarized else 'No'}",
        "",
    ]

    def row(label: str, r: ChannelGap) -> str:
        if r.metallic:
            return f"  {label:<22} 0.0000 eV   {'Metallic':<10}   -"
        if r.gap is None:
            return f"  {label:<22} n/a          {'Unknown':<10}   (no empty bands in file)"
        dtype = "Direct" if r.direct else "Indirect"
        trans = f"{kpoint_name(bs, r.vbm_k)} → {kpoint_name(bs, r.cbm_k)}"
        return f"  {label:<22} {r.gap:.4f} eV   {dtype:<10}   {trans}"

    if bs.is_spin_polarized:
        lines += [
            row("Spin-up band gap:", per_spin[0]),
            row("Spin-down band gap:", per_spin[1]),
            "",
            row("Overall band gap:", overall),
        ]
    else:
        lines.append(row("Band gap:", overall))
    return lines


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _path_geometry(bs: BandStructure):
    """x-coordinates along the path, branch index ranges and tick positions."""
    cart = bs.kpoints @ bs.reciprocal
    dx = np.linalg.norm(np.diff(cart, axis=0), axis=1)
    labelled = [lab is not None for lab in bs.labels]

    # Discontinuities: explicit BREAKs from the .cell file (k-points between
    # the two ends are dropped) and two adjacent labelled k-points.
    jumps = list(bs.breaks)
    jumps += [(i, i + 1) for i in range(len(dx))
              if labelled[i] and labelled[i + 1] and (i, i + 1) not in jumps]
    jumps.sort()
    for i, j in jumps:
        dx[i:j] = 0.0
    x = np.concatenate([[0.0], np.cumsum(dx)])

    branches = []
    start = 0
    for i, j in jumps:
        branches.append((start, i + 1))
        start = j
    branches.append((start, len(x)))

    valid = bs.valid
    ticks: List[Tuple[float, str]] = []
    for i, lab in enumerate(bs.labels):
        if lab is None or not valid[i]:
            continue
        if ticks and abs(ticks[-1][0] - x[i]) < 1e-9:
            if ticks[-1][1] != lab:
                ticks[-1] = (x[i], f"{ticks[-1][1]}|{lab}")
        else:
            ticks.append((x[i], lab))
    return x, branches, ticks


def plot_band_structure(bs: BandStructure, emin: float = -6.0, emax: float = 6.0):
    """
    Build a matplotlib Figure of the band structure.

    Energies are plotted relative to the VBM for semiconductors/insulators and
    relative to the Fermi energy for metals. Valence bands are blue and
    conduction bands orange; metals are drawn entirely in blue. For
    spin-polarised calculations spin-up is solid orange and spin-down dashed
    blue.
    """
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    per_spin, overall = analyse_gaps(bs)
    if overall.metallic or overall.vbm is None:
        zero = bs.efermi[0]
        zero_name = "Fermi level"
    else:
        zero = overall.vbm
        zero_name = "VBM"

    x, branches, ticks = _path_geometry(bs)

    fig = Figure(figsize=(5.5, 5.0))
    ax = fig.add_subplot(111)

    blue, orange = "#3953a4", "#f5a31a"
    for spin, E in bs.eigenvalues.items():
        res = per_spin[spin]
        for b in range(E.shape[0]):
            if bs.is_spin_polarized:
                colour, ls = (orange, "-") if spin == 0 else (blue, "--")
            elif res.metallic or res.n_occ is None:
                colour, ls = blue, "-"
            else:
                colour, ls = (blue if b < res.n_occ else orange), "-"
            for lo, hi in branches:
                ax.plot(x[lo:hi], E[b, lo:hi] - zero, color=colour, ls=ls, lw=1.6)

    for pos, _ in ticks:
        ax.axvline(pos, color="black", lw=0.8)
    ax.axhline(0.0, color="0.5", lw=0.6, ls=":")
    ax.tick_params(direction="in", top=False, right=True, width=1.0)
    for sp in ax.spines.values():
        sp.set_linewidth(1.0)

    ax.set_xticks([p for p, _ in ticks])
    ax.set_xticklabels([_mathtext_label(l) for _, l in ticks])
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(emin, emax)
    ax.set_ylabel("Energy (eV)")
    if bs.is_spin_polarized:
        ax.legend(
            handles=[
                Line2D([], [], color=orange, ls="-", label="Spin up"),
                Line2D([], [], color=blue, ls="--", label="Spin down"),
            ],
            loc="upper right", frameon=False,
        )
    fig.tight_layout()
    return fig
