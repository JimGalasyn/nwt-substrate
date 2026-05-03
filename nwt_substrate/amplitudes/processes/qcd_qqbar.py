"""
QCD scattering:  q(p1) + qbar(p2)  ->  q'(p3) + qbar'(p4)   (different flavor)

Single tree-level diagram: s-channel gluon exchange.

Substrate-algebra expression:

    iM = (-i g_s^2 / s) sum_a (T^a)_{i' i} (T^a)_{j j'}
                     * [vbar(p2) gamma^mu u(p1)] [ubar(p3) gamma_mu v(p4)]

Color sum / average for spin-averaged |M|^2:

    sum over all color indices and gluon color a:
        Sum_{i,i',j,j',a} |T^a_{i'i} T^a_{jj'}|^2
        = Sum_a [Tr(T^a T^a)] x [Tr(T^a T^a)]
        = Sum_a (1/2 * 1/2) = 8 * 1/4 = 2

    averaged over initial colors (1/N_c)^2 = 1/9:
        color factor = (1/9) * 2 = 2/9

QED-analog kinematic factor: same as e+e- -> mu+mu- (factorized two-trace).

Textbook QCD formula (massless limit, Peskin chapter 17 form):

    |M|^2_avg = (2/9) * g_s^4 * (1 + cos^2(theta))
              = (4/9) * g_s^4 * (t^2 + u^2)/s^2

For massive final-state quarks (still massless initial), replace
the kinematic factor with the e+e- -> mu+mu- form.

Substrate verification: structurally the same as e+e- -> mu+mu- (one
diagram, two factorized fermion-line traces, no interference issue),
so the (2E)^N normalization handles spin-counting correctly and the
color factor enters as a multiplicative 2/9.
"""

from __future__ import annotations

import numpy as np

from . import eemumu
from ...algebra.su3 import N_C


LATEX_EXPRESSION = (
    r"$i\mathcal{M} = \frac{-i\,g_s^2}{s}\,"
    r"(T^a)_{i'i}(T^a)_{jj'}\,"
    r"[\bar{v}_2\gamma^{\mu}u_1][\bar{u}_3\gamma_{\mu}v_4]$"
)


def color_factor_qqbar() -> float:
    """
    Color average + sum factor for q qbar -> q' qbar' s-channel:

        sum_a [Tr(T^a T^a)]^2 / (N_c)^2 = 8 * (1/4) / 9 = 2/9
    """
    return 2.0 / 9.0


def M_squared_avg(E_cm: float, theta: float, m_q: float, m_qprime: float,
                  gammas, alpha_s: float = 0.1179) -> tuple:
    """
    Spin- and color-averaged |M|^2 for q qbar -> q' qbar' (different flavor)
    via s-channel gluon, computed using substrate Dirac structure.

    Reuses the e+e- -> mu+mu- substrate amplitude and applies the QCD
    couplings + color factor:

        |M|^2_QCD = (2/9) * (g_s/e)^4 * |M|^2_QED_eemumu

    Args:
        E_cm     : CM energy (GeV)
        theta    : scattering angle in CM frame
        m_q      : initial-state quark mass (use small value for ~massless)
        m_qprime : outgoing quark mass
        gammas   : Dirac matrices (from algebra.dirac.make_lorentz_gammas)
        alpha_s  : strong coupling (default M_Z value)

    Returns:
        (|M|^2_avg, p1, p2, p3, p4)
    """
    # Use the e+e- -> mu+mu- substrate amplitude with electron charge
    # e_dummy = 1 so we can replace with g_s explicitly.
    e_dummy = 1.0
    M_sq_no_e, p1, p2, p3, p4 = eemumu.M_squared_avg(
        E_cm, theta, m_q, m_qprime, gammas, e=e_dummy
    )
    # M_sq_no_e includes a factor of e^4 = 1 (from e_dummy = 1).
    # Replace with g_s^4 = (4 pi alpha_s)^2 = 16 pi^2 alpha_s^2.
    g_s_sq = 4.0 * np.pi * alpha_s
    g_s_4 = g_s_sq * g_s_sq

    # Color factor: 2/9 (color avg + sum)
    color = color_factor_qqbar()

    return color * g_s_4 * M_sq_no_e, p1, p2, p3, p4


def textbook_M_squared(E_cm: float, theta: float, m_qprime: float,
                       alpha_s: float = 0.1179) -> float:
    """
    Textbook QCD |M|^2_avg for q qbar -> q' qbar' (different flavor),
    massless final-state limit:

        |M|^2_avg = (2/9) * g_s^4 * [(1 + 4 m_qprime^2/s) + (1 - 4 m_qprime^2/s) cos^2(theta)]

    For m_qprime = 0:
        |M|^2_avg = (2/9) * g_s^4 * (1 + cos^2(theta))
    """
    s = E_cm ** 2
    g_s_sq = 4.0 * np.pi * alpha_s
    g_s_4 = g_s_sq * g_s_sq
    r = 4.0 * (m_qprime ** 2) / s
    cos2 = float(np.cos(theta)) ** 2
    return color_factor_qqbar() * g_s_4 * ((1.0 + r) + (1.0 - r) * cos2)
