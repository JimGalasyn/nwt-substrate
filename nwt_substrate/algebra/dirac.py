"""
Lorentzian Dirac matrices via Wick rotation of Cl(0,7).

Picking 4 of the 7 imaginary octonion directions and Wick-rotating one gives
Lorentz Cl(1,3) with signature (+,-,-,-).  Default convention:

    gamma^0 = i L_4
    gamma^i = L_i  for i = 1, 2, 3

so {gamma^mu, gamma^nu} = 2 eta^{mu nu} I_8 with eta = diag(+,-,-,-).

The remaining 3 imaginary directions (L_5, L_6, L_7) generate an
internal SU(2) that commutes with all gamma^mu.

Verified in tests/test_dirac.py to machine precision.
"""

from __future__ import annotations

import numpy as np

from .clifford import left_mult_matrix


# Lorentzian metric (mostly-minus convention)
ETA_LORENTZ = np.diag([1.0, -1.0, -1.0, -1.0])


def make_lorentz_gammas(T: np.ndarray, lorentz_indices=(1, 2, 3, 4)):
    """
    Build Lorentzian Dirac matrices [gamma^0, gamma^1, gamma^2, gamma^3].

    Parameters
    ----------
    T : ndarray
        Octonion multiplication table (from algebra.octonions.make_octonion_table).
    lorentz_indices : tuple of 4 ints
        Octonion directions to use for Lorentz.  Default (1,2,3,4) means:
        spatial = (e_1, e_2, e_3), time = e_4 (Wick rotated).

        Must be 4 distinct integers in {1, ..., 7} (= N_CL07_IMAGINARY).

    Returns
    -------
    list of 4 complex 8x8 ndarrays
        [gamma^0, gamma^1, gamma^2, gamma^3] satisfying
        {gamma^mu, gamma^nu} = 2 eta^{mu nu} I_8.

    Substrate identities asserted at runtime:
        len(lorentz_indices) == N_LORENTZ_FROM_CL07         (= 4)
        each index in {1, ..., N_CL07_IMAGINARY}             (= 7)
        each gamma has shape (DIM_OCTONION, DIM_OCTONION)   (= 8x8)

    The remaining 3 imaginary directions (e.g. e_5, e_6, e_7) generate
    an internal SU(2) of dimension DIM_INTERNAL_SU2 = 3 that commutes
    with all γ^μ.
    """
    from ..isa.constants import (
        N_LORENTZ_FROM_CL07,
        N_CL07_IMAGINARY,
        DIM_OCTONION,
    )
    if len(lorentz_indices) != N_LORENTZ_FROM_CL07:
        raise ValueError(
            f"lorentz_indices must have length {N_LORENTZ_FROM_CL07} "
            f"(= N_LORENTZ_FROM_CL07 = N_CL07_IMAGINARY - "
            f"DIM_INTERNAL_SU2); got {len(lorentz_indices)}"
        )
    for idx in lorentz_indices:
        if not (1 <= int(idx) <= N_CL07_IMAGINARY):
            raise ValueError(
                f"lorentz_indices values must be in [1, {N_CL07_IMAGINARY}] "
                f"(= N_CL07_IMAGINARY); got {idx}"
            )

    spatial = [int(idx) for idx in lorentz_indices[:3]]
    time = int(lorentz_indices[3])
    L_spatial = [left_mult_matrix(idx, T).astype(complex) for idx in spatial]
    L_time = left_mult_matrix(time, T).astype(complex)
    gamma_0 = 1j * L_time
    gammas = [gamma_0] + L_spatial

    # Cross-shim shape identity: every gamma is DIM_OCTONION × DIM_OCTONION
    for k, g in enumerate(gammas):
        if g.shape != (DIM_OCTONION, DIM_OCTONION):
            raise AssertionError(
                f"Substrate identity violated: gamma^{k} has shape "
                f"{g.shape}, expected ({DIM_OCTONION}, {DIM_OCTONION}) "
                f"(= DIM_OCTONION × DIM_OCTONION)"
            )
    return gammas


def slash_p(p_mu: np.ndarray, gammas) -> np.ndarray:
    """
    Build slash{p} = gamma^mu p_mu using mostly-minus metric.

    Parameters
    ----------
    p_mu : 4-vector (E, p_x, p_y, p_z)  -- contravariant components
    gammas : list of 4 Dirac matrices from make_lorentz_gammas

    Returns
    -------
    8x8 complex ndarray
        slash{p} satisfying (slash{p})^2 = p^2 I_8.
    """
    sp = gammas[0] * p_mu[0]
    for i in range(1, 4):
        sp = sp - gammas[i] * p_mu[i]
    return sp


def make_gamma5(gammas) -> np.ndarray:
    """
    gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3.

    Anticommutes with all gamma^mu, squares to I_8, has zero trace.
    """
    return 1j * gammas[0] @ gammas[1] @ gammas[2] @ gammas[3]


def slash_eps(eps_mu: np.ndarray, gammas) -> np.ndarray:
    """
    slash{epsilon} = gamma^mu epsilon_mu (for photon polarization or any 4-vec).
    """
    sl = gammas[0] * eps_mu[0]
    for i in range(1, 4):
        sl = sl - gammas[i] * eps_mu[i]
    return sl


def lorentz_dot(p, q) -> float | complex:
    """p . q = p^0 q^0 - p^i q^i (mostly-minus metric)."""
    return p[0] * q[0] - sum(p[i] * q[i] for i in range(1, 4))
