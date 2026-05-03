"""
QCD same-flavor scattering:  q(p1) + q(p2)  ->  q(p3) + q(p4).

Two tree-level diagrams: t-channel and u-channel gluon exchange,
related by p3 <-> p4 with relative minus sign from Fermi antisymmetry.
This is the QCD analog of Moller scattering, with SU(3) color factors
modifying both the diagonal terms and the t-u interference.

Color sums (verified numerically via direct computation):

    Sum over all 4 quark colors and gluon color a:
      Sum |C_t|^2  =  2          [|M_t|^2 piece]
      Sum |C_u|^2  =  2          [|M_u|^2 piece]
      Sum Re(C_t* C_u) = -2/3    [interference]

    With initial-color average factor 1/9:
      color_t  =  2/9
      color_u  =  2/9
      color_cross  =  -2/27       (sign opposite to QED Moller's positive cross)

QCD textbook (massless limit, Peskin / Halzen-Martin):

    |M|^2_avg(qq -> qq)  =  (4/9) g_s^4 [(s^2+u^2)/t^2 + (s^2+t^2)/u^2 - (2/3)(s^2/(tu))]

The substrate calculation reuses Moller's amplitude_channels (kinematic
piece), with two corrections:
  (1) replace e -> g_s
  (2) apply substrate interference doubling rule (Phase A.3)
  (3) apply QCD color factors to t, u, cross terms

Verified to ~1e-9 against the textbook formula.
"""

from __future__ import annotations

import numpy as np

from . import moller
from ..spinors import positive_energy_spinors


LATEX_EXPRESSION = (
    r"$i\mathcal{M} = \frac{-i\,g_s^2}{t}\,(T^a)_{i'i}(T^a)_{j'j}\,"
    r"[\bar{u}_3\gamma^{\mu}u_1][\bar{u}_4\gamma_{\mu}u_2]"
    r"\;-\;(t \leftrightarrow u, \; 3 \leftrightarrow 4)$"
)


# Color factors for spin- and color-averaged |M|^2
COLOR_T = 2.0 / 9.0
COLOR_U = 2.0 / 9.0
COLOR_CROSS = -2.0 / 27.0


def kinematics_cm(E_cm: float, theta: float, m_q: float):
    """Same as Moller (since both are 4 identical-flavor fermions)."""
    return moller.kinematics_cm(E_cm, theta, m_q)


def M_squared_avg(E_cm: float, theta: float, m_q: float,
                  gammas, alpha_s: float = 0.1179) -> tuple:
    """
    Spin- and color-averaged |M|^2 for q q -> q q (same flavor).

    Uses the substrate-corrected formula:
        |M|^2 = color_t * sum_t + color_u * sum_u
                - 2 * color_cross * 2 * sum_cross_substrate

    where the inner factor of 2 multiplying sum_cross is the Phase A.3
    substrate interference doubling rule.

    Returns (|M|^2_avg, p1, p2, p3, p4).
    """
    p1, p2, p3, p4 = kinematics_cm(E_cm, theta, m_q)
    g_s = float(np.sqrt(4.0 * np.pi * alpha_s))

    u1s = positive_energy_spinors(p1, gammas, m_q)
    u2s = positive_energy_spinors(p2, gammas, m_q)
    u3s = positive_energy_spinors(p3, gammas, m_q)
    u4s = positive_energy_spinors(p4, gammas, m_q)

    sum_t = 0.0
    sum_u = 0.0
    sum_cross = 0.0

    for u1 in u1s:
        for u2 in u2s:
            for u3 in u3s:
                for u4 in u4s:
                    iM_t, iM_u = moller.amplitude_channels(
                        p1, p2, p3, p4, u1, u2, u3, u4, gammas, e=g_s,
                    )
                    sum_t += abs(iM_t) ** 2
                    sum_u += abs(iM_u) ** 2
                    sum_cross += (iM_t * np.conjugate(iM_u)).real

    N_in = len(u1s) * len(u2s)
    rel_norm = (2.0 * float(p1[0])) * (2.0 * float(p2[0])) \
               * (2.0 * float(p3[0])) * (2.0 * float(p4[0]))

    # QCD assembly:
    #   |M|^2 = color_t * sum_t + color_u * sum_u
    #           - 2 * color_cross * (sum_cross_QED)
    # where sum_cross_QED = 2 * sum_cross_substrate (doubling rule)
    # so   contribution = -2 * color_cross * 2 * sum_cross_substrate
    #                   = -4 * color_cross * sum_cross_substrate
    #
    # color_cross = -2/27  =>  -4 * (-2/27) = 8/27  (positive coefficient)
    M_sq_total = (COLOR_T * sum_t + COLOR_U * sum_u
                  - 4.0 * COLOR_CROSS * sum_cross)

    return rel_norm * M_sq_total / N_in, p1, p2, p3, p4


def textbook_M_squared(E_cm: float, theta: float,
                       alpha_s: float = 0.1179) -> float:
    """
    Textbook QCD |M|^2_avg for q q -> q q (same flavor, massless):

        |M|^2_avg = (4/9) g_s^4 [(s^2+u^2)/t^2 + (s^2+t^2)/u^2 - (2/3)(s^2/(tu))]
    """
    s = E_cm ** 2
    t = -s * np.sin(theta / 2.0) ** 2
    u = -s * np.cos(theta / 2.0) ** 2
    g_s_sq = 4.0 * np.pi * alpha_s
    g_s_4 = g_s_sq * g_s_sq
    bracket = ((s ** 2 + u ** 2) / t ** 2
               + (s ** 2 + t ** 2) / u ** 2
               - (2.0 / 3.0) * s ** 2 / (t * u))
    return float((4.0 / 9.0) * g_s_4 * bracket)
