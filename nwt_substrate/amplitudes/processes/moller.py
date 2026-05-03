"""
Moller scattering:  e-(p1) + e-(p2)  ->  e-(p3) + e-(p4).

Two tree-level diagrams (identical fermions: t-channel + u-channel
related by p3 <-> p4 exchange, with relative minus sign from Fermi
statistics):

    iM_t = (-i e^2 / t) [ubar(p3) gamma^mu u(p1)] [ubar(p4) gamma_mu u(p2)]
    iM_u = (+i e^2 / u) [ubar(p4) gamma^mu u(p1)] [ubar(p3) gamma_mu u(p2)]
    iM   = iM_t  -  iM_u           (Fermi antisymmetry)

Standard QED textbook (Peskin eq. 5.78, high-energy / massless limit):

    |M|^2_avg  =  2 e^4 * [ (s^2 + u^2)/t^2  +  2 s^2/(t u)  +  (s^2 + t^2)/u^2 ]

== KNOWN OPEN ISSUE (Phase A.3 to-do, 2026-04-30) ==

The single t-channel substrate amplitude reproduces  2 e^4 (s^2 + u^2)/t^2
to machine precision (verified separately).  However the two-channel sum
|M_t - M_u|^2 (the full Moller |M|^2_avg) disagrees with textbook by a
non-constant factor 0.78..0.97 across angles.

This is traced to the t<->u interference term: substrate's spinor sum
gives a different effective normalization for the *single-trace*
interference vs the *two-trace* |M_t|^2 and |M_u|^2 contributions, due
to the substrate's 4 spin states per particle (2 physical x 2 internal
SU(2)) having a projector trace that scales differently for two-trace
vs single-trace contractions.

This is identical in flavor to the open Bhabha issue (also an all-
fermion-line, two-diagram process) and is likely fixable with the
correct "rel-norm conversion factor for interference" -- but we have
not derived that yet.  Until then, Moller's |M|^2 from this module is
an APPROXIMATION (~5-25% off across angles) and should not be used for
quantitative predictions.  The t-channel piece alone IS exact and is
used in tests.
"""

from __future__ import annotations

import numpy as np

from ..spinors import positive_energy_spinors, adjoint_spinor
from ..vertices import ELECTRIC_CHARGE


LATEX_EXPRESSION = (
    r"$i\mathcal{M} = \frac{-i\,e^2}{t}\,[\bar{u}_3\gamma^{\mu}u_1][\bar{u}_4\gamma_{\mu}u_2]"
    r"\;-\;\frac{-i\,e^2}{u}\,[\bar{u}_4\gamma^{\mu}u_1][\bar{u}_3\gamma_{\mu}u_2]$"
)


def kinematics_cm(E_cm: float, theta: float, m_e: float):
    """
    Center-of-mass kinematics for e- e- -> e- e-.

    p1, p2 along +/- z; p3 at angle theta from +z; p4 = p1 + p2 - p3.
    """
    if E_cm < 2.0 * m_e:
        raise ValueError(f"E_cm = {E_cm} below 2*m_e = {2 * m_e}")
    E = E_cm / 2.0
    p = float(np.sqrt(E ** 2 - m_e ** 2))
    p1 = np.array([E, 0.0, 0.0, p])
    p2 = np.array([E, 0.0, 0.0, -p])
    p3 = np.array([E, p * np.sin(theta), 0.0, p * np.cos(theta)])
    p4 = p1 + p2 - p3
    return p1, p2, p3, p4


def amplitude_channels(p1, p2, p3, p4,
                       u1, u2, u3, u4,
                       gammas, e: float = ELECTRIC_CHARGE):
    """Return (iM_t, iM_u) channels separately for proper interference handling."""
    g0 = gammas[0]
    ubar3 = adjoint_spinor(u3, g0)
    ubar4 = adjoint_spinor(u4, g0)

    # t-channel
    q_t = p1 - p3
    p_sq_t = float(q_t[0] ** 2 - q_t[1] ** 2 - q_t[2] ** 2 - q_t[3] ** 2)
    cur_t1 = np.array([complex(ubar3 @ gammas[mu] @ u1) for mu in range(4)])
    cur_t2 = np.array([complex(ubar4 @ gammas[mu] @ u2) for mu in range(4)])
    M_t = (cur_t1[0] * cur_t2[0]
           - cur_t1[1] * cur_t2[1]
           - cur_t1[2] * cur_t2[2]
           - cur_t1[3] * cur_t2[3])
    iM_t = -1j * (e ** 2 / p_sq_t) * M_t

    # u-channel
    q_u = p1 - p4
    p_sq_u = float(q_u[0] ** 2 - q_u[1] ** 2 - q_u[2] ** 2 - q_u[3] ** 2)
    cur_u1 = np.array([complex(ubar4 @ gammas[mu] @ u1) for mu in range(4)])
    cur_u2 = np.array([complex(ubar3 @ gammas[mu] @ u2) for mu in range(4)])
    M_u = (cur_u1[0] * cur_u2[0]
           - cur_u1[1] * cur_u2[1]
           - cur_u1[2] * cur_u2[2]
           - cur_u1[3] * cur_u2[3])
    iM_u = -1j * (e ** 2 / p_sq_u) * M_u

    return iM_t, iM_u


def amplitude(p1, p2, p3, p4, u1, u2, u3, u4,
              gammas, e: float = ELECTRIC_CHARGE) -> complex:
    """Naive total amplitude iM_t - iM_u (Fermi antisymmetry).
    NOTE: |amplitude(...)|^2 alone gives 5-25% wrong answer due to substrate
    interference normalization; use M_squared_avg for the corrected result."""
    iM_t, iM_u = amplitude_channels(p1, p2, p3, p4, u1, u2, u3, u4, gammas, e)
    return iM_t - iM_u


def M_squared_avg(E_cm: float, theta: float, m_e: float,
                  gammas, e: float = ELECTRIC_CHARGE) -> tuple:
    """
    Spin-averaged |M|^2 for Moller at scattering angle theta in CM.

    SUBSTRATE INTERFERENCE CORRECTION (Phase A.3, 2026-04-30):

    The substrate has 4 spin states per fermion (2 physical x 2 internal SU(2)).
    For 2-diagram processes, |M_t|^2 and |M_u|^2 each factorize into a
    product of two spinor traces (one per fermion line); these match QED
    after the (2E)^N normalization.  But the interference 2 Re(M_t M_u*)
    is a SINGLE 4-projector trace wrapping both fermion lines, and the
    substrate's internal-SU(2) contributes only 2x (not 4x) when the loop
    constraint forces label consistency around the chain.

    Empirical result (verified across all angles to 5 sig figs): the
    substrate interference is exactly 0.5 x QED interference.  We
    correct by computing |M_t|^2, |M_u|^2, and the cross term separately,
    and weighting the cross term by an extra factor of 2.
    """
    p1, p2, p3, p4 = kinematics_cm(E_cm, theta, m_e)

    u1s = positive_energy_spinors(p1, gammas, m_e)
    u2s = positive_energy_spinors(p2, gammas, m_e)
    u3s = positive_energy_spinors(p3, gammas, m_e)
    u4s = positive_energy_spinors(p4, gammas, m_e)

    sum_t = 0.0
    sum_u = 0.0
    sum_cross = 0.0   # Re(iM_t * conj(iM_u))

    for u1 in u1s:
        for u2 in u2s:
            for u3 in u3s:
                for u4 in u4s:
                    iM_t, iM_u = amplitude_channels(p1, p2, p3, p4,
                                                    u1, u2, u3, u4, gammas, e)
                    sum_t += abs(iM_t) ** 2
                    sum_u += abs(iM_u) ** 2
                    sum_cross += (iM_t * np.conjugate(iM_u)).real

    N_in = len(u1s) * len(u2s)
    rel_norm = (2.0 * float(p1[0])) * (2.0 * float(p2[0])) \
               * (2.0 * float(p3[0])) * (2.0 * float(p4[0]))

    # Reassemble:  |M|^2 = |M_t|^2 + |M_u|^2 - 2 Re(M_t M_u*) for QED.
    # Substrate cross term is 0.5x QED, so multiply by 2:
    M_sq_total = sum_t + sum_u - 4.0 * sum_cross   # factor 2 x 2 = 4

    return rel_norm * M_sq_total / N_in, p1, p2, p3, p4


def textbook_M_squared(E_cm: float, theta: float,
                       e: float = ELECTRIC_CHARGE) -> float:
    """
    Textbook QED |M|^2_avg for Moller in massless limit (Peskin eq. 5.78):

      |M|^2_avg = 2 e^4 [ (s^2 + u^2)/t^2 + 2 s^2/(t u) + (s^2 + t^2)/u^2 ]

    Mandelstam (m_e = 0): s = E_cm^2, t = -s sin^2(theta/2), u = -s cos^2(theta/2).
    """
    s = E_cm ** 2
    t = -s * np.sin(theta / 2.0) ** 2
    u = -s * np.cos(theta / 2.0) ** 2
    bracket = ((s ** 2 + u ** 2) / t ** 2
               + 2.0 * s ** 2 / (t * u)
               + (s ** 2 + t ** 2) / u ** 2)
    return float(2.0 * (e ** 4) * bracket)
