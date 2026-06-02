"""Geometric Whitehead Hopf charge Q_H of a lattice S^2 field.

The Hopf charge Q_H in pi_3(S^2) = Z is the linking number of two generic
preimages of the map n: S^3 -> S^2.  Computed here as the Whitehead integral:
form the magnetic 2-form B = (F_23, F_31, F_12) from the faithful area form
(see `area_form`), solve curl A = B spectrally in Coulomb gauge, and

    Q_H = (1 / 16 pi^2) integral A . B .

Using the area-form F (rather than the naive same-index product) is what makes
Q_H robust on a relaxed lattice field -- the same lesson that stabilises the
hopfion.

This is the FIELD-theoretic Hopf charge.  It is distinct from the algebraic
carrier-knot invariant `nwt_substrate.topology.hopf_charge(p, m) = p * m`, which
labels a (p, m) torus knot; the two coincide in value for the field that
realises that knot but are computed from entirely different data.
"""
from __future__ import annotations

import numpy as np

from .area_form import plaquette_area_form


def whitehead_hopf_charge(n1, n2, n3, grid):
    """Whitehead Hopf charge Q_H = (1/16 pi^2) int A.B of the unit field n.

    grid is a `BoxGrid` (supplies the spectral wavenumbers and dx).  Returns a
    float that is an integer for a smooth, well-resolved hopfion."""
    KX, KY, KZ, K2 = grid.k_vectors()
    dx = grid.dx
    n = np.stack([n1, n2, n3])
    Bx = plaquette_area_form(n, 1, 2) / dx**2
    By = plaquette_area_form(n, 2, 0) / dx**2
    Bz = plaquette_area_form(n, 0, 1) / dx**2
    Bxh, Byh, Bzh = np.fft.fftn(Bx), np.fft.fftn(By), np.fft.fftn(Bz)
    inv = np.zeros_like(K2)
    np.divide(1.0, K2, out=inv, where=K2 > 0)       # skip the k=0 mode
    Ax = np.real(np.fft.ifftn(-1j * (KY * Bzh - KZ * Byh) * inv))
    Ay = np.real(np.fft.ifftn(-1j * (KZ * Bxh - KX * Bzh) * inv))
    Az = np.real(np.fft.ifftn(-1j * (KX * Byh - KY * Bxh) * inv))
    H = np.sum(Ax * Bx + Ay * By + Az * Bz) * dx**3
    return float(H / (16.0 * np.pi**2))
