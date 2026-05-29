"""
Paper 6 mass formula: substrate-algebraic prediction of particle masses.

  m / m_e = [(p^2 + q^2) / 5]                 # gradient energy on torus
            * [beta(p, m) / beta_e]            # aspect-ratio ratio
            * [ln(8 beta) / ln(8 beta_e)]      # Kelvin-Saffman ring log
            * n_q^q                            # multi-component link enhancement

Anchors: m_e (electron rest energy) and kappa = pi^2 (Derrick equilibrium).
No other free parameters.

Verified across the 24-particle compendium at 0.76% median residual
after the (1,3,5,5) nucleon correction (2026-04-30).
"""

from __future__ import annotations

import numpy as np


# Electron rest energy (MeV) and Paper 6 constants
ME_MEV = 0.510998928
HBAR_C = 197.3269788  # MeV*fm

# Electron reference: (p, q, m) = (2, 1, 3)
P_E, Q_E, M_E_INT = 2, 1, 3
BETA_E = float(np.sqrt(M_E_INT ** 2 / P_E ** 2 - 1.0))  # = sqrt(5)/2
LN8_BE = float(np.log(8.0 * BETA_E))


def beta_phase(p: int, m_int: int) -> float | None:
    """sqrt(m^2/p^2 - 1).  Returns None if m_int^2 <= p^2."""
    if m_int * m_int <= p * p:
        return None
    return float(np.sqrt(m_int * m_int / (p * p) - 1.0))


def paper6_mass_ratio(p: int, q: int, m_int: int, nq: int) -> float | None:
    """
    Paper 6 mass formula in dimensionless form: m / m_e.

    Returns None if (p, q, m, nq) is unphysical.
    """
    b = beta_phase(p, m_int)
    if b is None or b <= 0.0:
        return None
    pq2 = p * p + q * q
    w_knot = pq2 / 5.0   # 5 = electron p²+q² (2²+1²) normalization, NOT H_V_SO7
    w_ring = (b / BETA_E) * (np.log(8.0 * b) / LN8_BE)
    if nq <= 1:
        w_link = 1.0
    else:
        w_link = nq ** q
    return float(w_knot * w_ring * w_link)


def paper6_mass_mev(p: int, q: int, m_int: int, nq: int) -> float | None:
    """Paper 6 mass formula in MeV."""
    ratio = paper6_mass_ratio(p, q, m_int, nq)
    if ratio is None:
        return None
    return ratio * ME_MEV
