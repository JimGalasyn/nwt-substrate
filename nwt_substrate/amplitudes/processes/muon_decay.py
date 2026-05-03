"""
Muon decay  mu- -> e- nu_mu nu_e_bar  via 4-fermion V-A interaction.

Effective theory amplitude (Fermi):

    iM = -(i G_F/sqrt(2)) [ubar(p_nu_mu) gamma^mu (1-gamma^5) u(p_mu)]
                          [ubar(p_e)     gamma_mu (1-gamma^5) v(p_nu_e_bar)]

Spin-averaged |M|^2 in the massless final-state limit (Peskin eq. 11.66):

    |M|^2_avg = 64 G_F^2 (p_mu . p_nu_e_bar) (p_e . p_nu_mu)

Reduced over the Dalitz plot in the muon rest frame with x = 2 E_e/m_mu,
y = 2 E_nu_e/m_mu, z = 2 - x - y:

    p_mu . p_nu_e_bar = m_mu^2 y / 2
    p_e . p_nu_mu     = m_mu^2 (1-y) / 2

so

    |M|^2_avg(x, y) = 16 G_F^2 m_mu^4 * y * (1-y)

The rate then integrates over the Dalitz region x + y >= 1 with x, y in [0, 1].

Standard textbook result:

    Gamma = G_F^2 m_mu^5 / (192 pi^3)   (massless final states)

Tau_mu = 1/Gamma ~ 2.197 us.

This module provides:
  - the Michel differential spectrum (1D)
  - the Dalitz density (2D)
  - both integrated to Gamma, with checks against the closed form
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad, dblquad

from ..cross_sections import G_F


LATEX_EXPRESSION = (
    r"$i\mathcal{M} = \frac{-i\,G_F}{\sqrt{2}}\,"
    r"\bar{u}_{\nu_\mu}\gamma^{\mu}(1-\gamma^5)u_\mu\,"
    r"\bar{u}_e\gamma_\mu(1-\gamma^5)v_{\bar{\nu}_e}$"
)


# ---------------------------------------------------------------------------
# Michel differential spectrum (1D in x_e)
# ---------------------------------------------------------------------------

def dGamma_dx_michel(x: float, m_mu: float, GF: float = G_F) -> float:
    """
    Michel spectrum  dGamma/dx_e = (G_F^2 m_mu^5)/(96 pi^3) * x^2 (3 - 2x).

    Massless final-state limit; ignores polarization (already spin-averaged
    over muon).  Integrates to G_F^2 m_mu^5 / (192 pi^3) over x in [0, 1].
    """
    if x < 0 or x > 1:
        return 0.0
    return (GF ** 2 * m_mu ** 5) / (96.0 * np.pi ** 3) * x ** 2 * (3.0 - 2.0 * x)


def gamma_via_michel(m_mu: float = 0.10566, GF: float = G_F) -> float:
    """Total muon decay rate by 1D integration of the Michel spectrum."""
    integral, _ = quad(lambda x: dGamma_dx_michel(x, m_mu, GF), 0.0, 1.0)
    return integral


# ---------------------------------------------------------------------------
# Full 2D Dalitz density
# ---------------------------------------------------------------------------

def dalitz_density_VA(x: float, y: float, m_mu: float,
                      GF: float = G_F) -> float:
    """
    Differential decay rate  d^2 Gamma / (dx dy)  on the muon Dalitz plot,
    with x = 2 E_e/m_mu, y = 2 E_nu_e/m_mu.

    Derivation (V-A with (1 - gamma^5) un-normalized convention):
        |M|^2_avg = 32 G_F^2 m_mu^4 * y (1 - y)            (massless final)
        dPS_3    = (m_mu^2 / (32 (2 pi)^3)) dx dy          (rest frame)
        d Gamma / (dx dy) = (1/(2 m_mu)) * |M|^2 * dPS_3
                          = (G_F^2 m_mu^5) / (2 (2 pi)^3) * y (1 - y)
                          = (G_F^2 m_mu^5) / (16 pi^3)   * y (1 - y)

    Dalitz region: x + y >= 1 with x, y in [0, 1].  Outside, returns 0.

    Integrates over the Dalitz triangle to G_F^2 m_mu^5 / (192 pi^3).
    """
    if x < 0 or x > 1 or y < 0 or y > 1 or (x + y) < 1.0:
        return 0.0
    return (GF ** 2 * m_mu ** 5) / (16.0 * np.pi ** 3) * y * (1.0 - y)


def gamma_via_dalitz(m_mu: float = 0.10566, GF: float = G_F) -> float:
    """
    Total muon decay rate by 2D integration over the Dalitz plot.

    Integration region: x in [0, 1], y in [max(0, 1 - x), 1].
    """
    def integrand(y, x):
        return dalitz_density_VA(x, y, m_mu, GF)

    # dblquad signature: dblquad(f, a, b, gfun, hfun)  integrates over x in [a,b],
    # then y in [gfun(x), hfun(x)].
    integral, _ = dblquad(integrand,
                          0.0, 1.0,                    # x range
                          lambda x: max(0.0, 1.0 - x),  # y_lower
                          lambda x: 1.0)                # y_upper
    return integral


def closed_form(m_mu: float = 0.10566, GF: float = G_F) -> float:
    """Closed-form textbook  Gamma = G_F^2 m_mu^5 / (192 pi^3)."""
    return (GF ** 2) * (m_mu ** 5) / (192.0 * np.pi ** 3)


# ---------------------------------------------------------------------------
# Reporting helper
# ---------------------------------------------------------------------------

def lifetime_seconds(Gamma_GeV: float) -> float:
    """Convert Gamma in GeV to lifetime in seconds (1 GeV = 1.519e24 / s)."""
    return 1.0 / (Gamma_GeV * 1.519267e24)
