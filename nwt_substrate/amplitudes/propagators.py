"""
Feynman propagators for substrate-algebraic amplitude calculations.

Conventions (mostly-minus metric eta = diag(+,-,-,-)):

  Fermion (Dirac):
      S_F(P) = (slash{P} + m I_8) / (P^2 - m^2)

  Scalar:
      D_F(P) = 1 / (P^2 - m^2)

  Photon (Feynman gauge):
      D_F^{mu nu}(P) = -eta^{mu nu} / P^2

  Massive vector (unitary gauge):
      D_F^{mu nu}(P) = -(eta^{mu nu} - P^mu P^nu / m^2) / (P^2 - m^2)
"""

from __future__ import annotations

import numpy as np

from ..algebra.dirac import slash_p, ETA_LORENTZ


def _p_dot_p(p_mu: np.ndarray) -> complex:
    """Mostly-minus inner product P.P = P^0 P^0 - P^i P^i."""
    return complex(p_mu[0] ** 2 - p_mu[1] ** 2 - p_mu[2] ** 2 - p_mu[3] ** 2)


def fermion_propagator(p_mu: np.ndarray, gammas, m: float) -> np.ndarray:
    """
    Dirac fermion Feynman propagator:  S_F(p) = (slash{p} + m I_8) / (p^2 - m^2).

    Returns 8x8 complex matrix.  No iε prescription (good for tree level).
    """
    sp = slash_p(p_mu, gammas)
    I = np.eye(8, dtype=complex)
    p_sq = _p_dot_p(p_mu)
    denom = p_sq - m ** 2
    if abs(denom) < 1e-15:
        raise ValueError(f"On-shell propagator pole: P^2 = m^2 = {p_sq}")
    return (sp + m * I) / denom


def scalar_propagator(p_mu: np.ndarray, m: float) -> complex:
    """Scalar Feynman propagator:  D_F(p) = 1 / (p^2 - m^2)."""
    p_sq = _p_dot_p(p_mu)
    denom = p_sq - m ** 2
    if abs(denom) < 1e-15:
        raise ValueError(f"On-shell propagator pole: P^2 = m^2 = {p_sq}")
    return 1.0 / denom


def photon_propagator(p_mu: np.ndarray) -> np.ndarray:
    """
    Photon Feynman-gauge propagator:  D_F^{mu nu}(p) = -eta^{mu nu} / p^2.

    Returns 4x4 complex matrix (Lorentz indices).
    """
    p_sq = _p_dot_p(p_mu)
    if abs(p_sq) < 1e-15:
        raise ValueError("On-shell photon (p^2=0); use external polarization not propagator.")
    return -ETA_LORENTZ.astype(complex) / p_sq


def massive_vector_propagator(p_mu: np.ndarray, m: float) -> np.ndarray:
    """
    Massive vector boson (unitary gauge):
        D_F^{mu nu}(p) = -(eta^{mu nu} - p^mu p^nu / m^2) / (p^2 - m^2).
    """
    p_sq = _p_dot_p(p_mu)
    denom = p_sq - m ** 2
    if abs(denom) < 1e-15:
        raise ValueError(f"On-shell massive vector: P^2 = m^2 = {p_sq}")
    eta = ETA_LORENTZ.astype(complex)
    pp = np.outer(p_mu, p_mu)
    return -(eta - pp / m ** 2) / denom
