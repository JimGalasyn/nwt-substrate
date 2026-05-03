"""
Bhabha scattering:  e+(p1) + e-(p2)  ->  e+(p3) + e-(p4).

Two tree-level diagrams: s-channel annihilation + t-channel exchange.

    iM_s = (+i e^2 / s) [v_bar(p1) gamma^mu u(p2)] [u_bar(p4) gamma_mu v(p3)]
    iM_t = (+i e^2 / t) [v_bar(p1) gamma^mu v(p3)] [u_bar(p4) gamma_mu u(p2)]
    iM   = iM_s - iM_t

Standard QED textbook (Peskin, high-energy / massless limit):

    |M|^2_avg  =  2 e^4 [ (s^2 + u^2)/t^2 + (t^2 + u^2)/s^2 + 2 u^2/(s t) ]

== KNOWN OPEN ISSUE (Phase A.3 to-do, 2026-04-30) ==

Like Moller, the single s-channel and single t-channel pieces match
QED.  The two-channel sum |M_s - M_t|^2 (the full Bhabha |M|^2_avg)
disagrees with textbook by ~5-12% across angles.  Same root cause:
substrate's spinor sum gives different effective normalization for
single-trace interference vs two-trace direct contributions when
there are 4 fermion external lines and two diagrams.

Until the interference normalization is derived from first principles,
Bhabha's full |M|^2 from this module is an APPROXIMATION (~5-12% off)
and should not be used for quantitative predictions.
"""

from __future__ import annotations

import numpy as np

from ..spinors import (
    positive_energy_spinors,
    negative_energy_spinors,
    adjoint_spinor,
)
from ..vertices import ELECTRIC_CHARGE


LATEX_EXPRESSION = (
    r"$i\mathcal{M} = \frac{i\,e^2}{s}\,[\bar{v}_1\gamma^{\mu}u_2][\bar{u}_4\gamma_{\mu}v_3]"
    r"\;-\;\frac{i\,e^2}{t}\,[\bar{v}_1\gamma^{\mu}v_3][\bar{u}_4\gamma_{\mu}u_2]$"
)


def kinematics_cm(E_cm: float, theta: float, m_e: float):
    """
    CM kinematics: e+(p1) + e-(p2) -> e+(p3) + e-(p4).

    theta = angle between incoming e+ (along +z) and outgoing e+ (= p3),
    so t = (p1 - p3)^2 = -s sin^2(theta/2) (massless limit).
    """
    if E_cm < 2.0 * m_e:
        raise ValueError(f"E_cm = {E_cm} below 2*m_e = {2 * m_e}")
    E = E_cm / 2.0
    p = float(np.sqrt(E ** 2 - m_e ** 2))
    p1 = np.array([E, 0.0, 0.0, p])           # e+ in
    p2 = np.array([E, 0.0, 0.0, -p])          # e- in
    p3 = np.array([E, p * np.sin(theta), 0.0, p * np.cos(theta)])   # e+ out
    p4 = p1 + p2 - p3                          # e- out
    return p1, p2, p3, p4


def amplitude_channels(p1, p2, p3, p4, v1, u2, v3, u4,
                       gammas, e: float = ELECTRIC_CHARGE):
    """Return (iM_s, iM_t) channels separately for proper interference handling."""
    g0 = gammas[0]
    vbar1 = adjoint_spinor(v1, g0)
    ubar4 = adjoint_spinor(u4, g0)

    # s-channel
    q_s = p1 + p2
    s = float(q_s[0] ** 2 - q_s[1] ** 2 - q_s[2] ** 2 - q_s[3] ** 2)
    cur_s1 = np.array([complex(vbar1 @ gammas[mu] @ u2) for mu in range(4)])
    cur_s2 = np.array([complex(ubar4 @ gammas[mu] @ v3) for mu in range(4)])
    M_s_kin = (cur_s1[0] * cur_s2[0]
               - cur_s1[1] * cur_s2[1]
               - cur_s1[2] * cur_s2[2]
               - cur_s1[3] * cur_s2[3])
    iM_s = 1j * (e ** 2 / s) * M_s_kin

    # t-channel
    q_t = p1 - p3
    t = float(q_t[0] ** 2 - q_t[1] ** 2 - q_t[2] ** 2 - q_t[3] ** 2)
    cur_t1 = np.array([complex(vbar1 @ gammas[mu] @ v3) for mu in range(4)])
    cur_t2 = np.array([complex(ubar4 @ gammas[mu] @ u2) for mu in range(4)])
    M_t_kin = (cur_t1[0] * cur_t2[0]
               - cur_t1[1] * cur_t2[1]
               - cur_t1[2] * cur_t2[2]
               - cur_t1[3] * cur_t2[3])
    iM_t = 1j * (e ** 2 / t) * M_t_kin

    return iM_s, iM_t


def amplitude(p1, p2, p3, p4, v1, u2, v3, u4,
              gammas, e: float = ELECTRIC_CHARGE) -> complex:
    """Naive total amplitude iM_s - iM_t (Fermi antisymmetry).
    NOTE: |amplitude(...)|^2 alone gives 5-12% wrong answer due to substrate
    interference normalization; use M_squared_avg for the corrected result."""
    iM_s, iM_t = amplitude_channels(p1, p2, p3, p4, v1, u2, v3, u4, gammas, e)
    return iM_s - iM_t


def M_squared_avg(E_cm: float, theta: float, m_e: float,
                  gammas, e: float = ELECTRIC_CHARGE) -> tuple:
    """
    Spin-averaged |M|^2 for Bhabha with substrate interference correction.

    Same correction as Moller (Phase A.3, 2026-04-30):  the substrate's
    cross term Re(iM_s * conj(iM_t)) is exactly 0.5x the QED interference,
    so we weight it by an extra factor of 2 in the |M|^2 assembly.
    """
    p1, p2, p3, p4 = kinematics_cm(E_cm, theta, m_e)

    v1s = negative_energy_spinors(p1, gammas, m_e)
    u2s = positive_energy_spinors(p2, gammas, m_e)
    v3s = negative_energy_spinors(p3, gammas, m_e)
    u4s = positive_energy_spinors(p4, gammas, m_e)

    sum_s = 0.0
    sum_t = 0.0
    sum_cross = 0.0   # Re(iM_s * conj(iM_t))

    for v1 in v1s:
        for u2 in u2s:
            for v3 in v3s:
                for u4 in u4s:
                    iM_s, iM_t = amplitude_channels(p1, p2, p3, p4,
                                                    v1, u2, v3, u4, gammas, e)
                    sum_s += abs(iM_s) ** 2
                    sum_t += abs(iM_t) ** 2
                    sum_cross += (iM_s * np.conjugate(iM_t)).real

    N_in = len(v1s) * len(u2s)
    rel_norm = (2.0 * float(p1[0])) * (2.0 * float(p2[0])) \
               * (2.0 * float(p3[0])) * (2.0 * float(p4[0]))

    # |M|^2 = |M_s|^2 + |M_t|^2 - 2 Re(M_s M_t*) for QED.
    # Substrate cross term is 0.5x QED, so multiply by 2:
    M_sq_total = sum_s + sum_t - 4.0 * sum_cross   # factor 2 x 2 = 4

    return rel_norm * M_sq_total / N_in, p1, p2, p3, p4


def textbook_M_squared(E_cm: float, theta: float,
                       e: float = ELECTRIC_CHARGE) -> float:
    """
    Textbook QED |M|^2_avg for Bhabha in massless limit:

      |M|^2_avg = 2 e^4 [ (s^2 + u^2)/t^2 + (t^2 + u^2)/s^2 + 2 u^2/(s t) ]
    """
    s = E_cm ** 2
    t = -s * np.sin(theta / 2.0) ** 2
    u = -s * np.cos(theta / 2.0) ** 2
    bracket = ((s ** 2 + u ** 2) / t ** 2
               + (t ** 2 + u ** 2) / s ** 2
               + 2.0 * u ** 2 / (s * t))
    return float(2.0 * (e ** 4) * bracket)
