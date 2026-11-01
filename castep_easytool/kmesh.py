# -*- coding: utf-8 -*-
"""
Lattice-aware Monkhorst-Pack grids for k-point convergence scans.

Convention (same as CASTEP's KPOINTS_MP_SPACING, in 1/Angstrom *without* the 2*pi):
    r_i = |b_i| / (2*pi) = length of the i-th row of inv(A).T
    n_i = ceil(r_i / spacing)

The scan is a ladder indexed by N, the number of k-points along the axis with the
longest reciprocal vector. The other axes scale in proportion and are rounded up:
    n_i = ceil(N * r_i / r_max),   spacing_N = r_max / N
For a cubic or fcc/bcc primitive cell this is exactly N x N x N; for hexagonal,
tetragonal, orthorhombic or low-symmetry cells the grid follows the cell shape.
"""
import math

import numpy as np

_EPS = 1e-9


def recip_lengths(lattice):
    """r_i = |b_i|/(2 pi) in 1/Angstrom for a 3x3 lattice (rows = lattice vectors, Angstrom)."""
    a = np.asarray(lattice, dtype=float)
    if a.shape != (3, 3) or abs(np.linalg.det(a)) < 1e-12:
        raise ValueError("invalid lattice (not a non-singular 3x3 matrix)")
    return np.linalg.norm(np.linalg.inv(a).T, axis=1)


def grid_from_spacing(lattice, spacing):
    """MP grid for a spacing in 1/Angstrom (CASTEP convention)."""
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    r = recip_lengths(lattice)
    return [max(1, int(math.ceil(x / spacing - _EPS))) for x in r]


def ladder_grid(lattice, n):
    """Grid with ``n`` points along the longest reciprocal axis, others in proportion."""
    r = recip_lengths(lattice)
    return [max(1, int(math.ceil(n * x / r.max() - _EPS))) for x in r]


def auto_n_max(lattice, spacing_min=0.02, cap=20, n_start=1):
    """Ladder size N at which the spacing reaches ``spacing_min`` (1/Angstrom), capped."""
    r = recip_lengths(lattice)
    return int(min(cap, max(n_start, math.ceil(r.max() / spacing_min - _EPS))))


def spacing_for_grid(lattice, grid):
    """
    A KPOINTS_MP_SPACING value that reproduces ``grid`` on this lattice under the
    ceil convention, chosen in the middle of the valid interval so that rounding in
    CASTEP cannot flip a grid dimension. Returns a float (1/Angstrom).
    """
    r = recip_lengths(lattice)
    n = np.asarray(grid, dtype=float)
    lo = float(np.max(r / n))
    upper = [x / (k - 1) for x, k in zip(r, n) if k > 1]
    hi = min(upper) if upper else math.inf
    if hi <= lo * (1 + 1e-6):
        raise ValueError(f"grid {list(grid)} cannot be produced by any spacing on this lattice")
    delta = min(0.02, (hi / lo - 1) / 2) if math.isfinite(hi) else 0.02
    target = lo * (1 + delta)
    for digits in (4, 5, 6, 8):
        v = round(target, digits)
        if lo * (1 + 1e-7) < v < hi * (1 - 1e-7) and list(grid_from_spacing(lattice, v)) == list(grid):
            return v
    return target


def per_axis_spacing(lattice, grid):
    """Actual k-point separation r_i / n_i along each reciprocal axis (1/Angstrom)."""
    return [float(x / n) for x, n in zip(recip_lengths(lattice), grid)]


def ladder(lattice, n_start=1, n_step=1, n_max=0, spacing_min=0.02):
    """
    List of (N, grid, spacing_N) from n_start to n_max (0 = automatic from
    ``spacing_min``). Consecutive duplicate grids are impossible by construction.
    """
    if n_max <= 0:
        n_max = auto_n_max(lattice, spacing_min, n_start=n_start)
    r_max = float(recip_lengths(lattice).max())
    return [(n, ladder_grid(lattice, n), r_max / n) for n in range(n_start, n_max + 1, n_step)]


def fixed_grid(lattice, spacing):
    """Grid for the cut-off scan (moderate fixed spacing)."""
    return grid_from_spacing(lattice, spacing)
