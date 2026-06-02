"""Faddeev-Skyrme energy and soliton diagnostics (numpy reference).

    E = E2 + c4 * E4,   E2 = sum_a |grad n^a|^2,   E4 = sum_{i<j} F_ij^2 ,

with E2 the FORWARD-difference Dirichlet energy and E4 the faithful area-form
quartic (see `area_form`).  Forward differences are used (not central) because
central differences have a Nyquist/checkerboard null space -- their derivative
of a +/-1 checkerboard is identically zero, so a central-difference energy does
not penalise checkerboard noise and a minimiser pumps in free sub-grid noise
that corrupts the charge and opens an unwinding channel.

Under x -> lambda x, E2 ~ lambda and E4 ~ 1/lambda (Derrick), so a soliton sits
at finite size lambda* = sqrt(c4 Ehat4 / E2) where E2 = E4 at the minimum
(virial).  `soliton_size` returns the E2/E4 split so that balance can be checked.
"""
from __future__ import annotations

import numpy as np

from .area_form import plaquette_area_form


def _fwd_grads(field, dx):
    """Forward (nearest-neighbour) difference gradient, periodic.  No Nyquist
    null space (unlike central differences)."""
    gx = (np.roll(field, -1, 0) - field) / dx
    gy = (np.roll(field, -1, 1) - field) / dx
    gz = (np.roll(field, -1, 2) - field) / dx
    return gx, gy, gz


def energy_density(n1, n2, n3, dx, c4):
    """Per-cell energy density e2 + c4 * e4 (its sum * dx**3 == faddeev_energy)."""
    g = [_fwd_grads(c, dx) for c in (n1, n2, n3)]
    e2 = sum(gx**2 + gy**2 + gz**2 for (gx, gy, gz) in g)
    n = np.stack([n1, n2, n3])
    e4 = np.zeros_like(n1)
    for (i, j) in ((0, 1), (1, 2), (0, 2)):
        e4 = e4 + (plaquette_area_form(n, i, j) / dx**2) ** 2
    return e2 + c4 * e4


def faddeev_energy(n1, n2, n3, dx, c4):
    """Total Faddeev energy E = E2 + c4 * E4 (forward-diff E2, area-form E4)."""
    return float(np.sum(energy_density(n1, n2, n3, dx, c4)) * dx**3)


def soliton_size(n1, n2, n3, grid, c4):
    """Energy-weighted RMS radius R_E about the energy centroid, plus the E2/E4
    split.  R_E scales linearly with the soliton; E2/E4 -> 1 at the Derrick
    minimum (virial), a convergence check.  Returns {"R_E", "E2", "E4"}."""
    dx = grid.dx
    X, Y, Z = (np.asarray(c) for c in grid.coords())
    e = energy_density(n1, n2, n3, dx, c4)
    g = [_fwd_grads(c, dx) for c in (n1, n2, n3)]
    E2 = float(np.sum(sum(gx**2 + gy**2 + gz**2 for (gx, gy, gz) in g)) * dx**3)
    n = np.stack([n1, n2, n3])
    e4 = sum((plaquette_area_form(n, i, j) / dx**2) ** 2
             for (i, j) in ((0, 1), (1, 2), (0, 2)))
    E4 = float(np.sum(e4) * dx**3) * c4
    w = e / e.sum()
    cx, cy, cz = (w * X).sum(), (w * Y).sum(), (w * Z).sum()
    r2 = (X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2
    return {"R_E": float(np.sqrt((w * r2).sum())), "E2": E2, "E4": E4}
