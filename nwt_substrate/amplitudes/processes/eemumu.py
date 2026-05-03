"""
e- + e+ -> mu- + mu-:  s-channel virtual photon.

Substrate-algebra expression:

    iM = (i e^2 / s)  [v_bar(p2) gamma^mu u(p1)]  [u_bar(p3) gamma_mu v(p4)]

where:
    p1 : incoming e-      (u spinor)
    p2 : incoming e+      (v spinor)
    p3 : outgoing mu-     (u spinor)
    p4 : outgoing mu+     (v spinor)
    s  = (p1 + p2)^2

Standard QED textbook result (spin-averaged):

    |M|^2_avg = (8 e^4 / s^2) * [ (p1.p3)(p2.p4) + (p1.p4)(p2.p3) + m_mu^2 (p1.p2) ]

In CM frame with m_e = 0 (high-energy limit):

    |M|^2_avg = e^4 * [(1 + 4 m_mu^2/s) + (1 - 4 m_mu^2/s) cos^2(theta)]

which reduces to e^4 * (1 + cos^2(theta)) when m_mu << sqrt(s)/2.
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
    r"$i\mathcal{M} = \frac{i\,e^2}{s}\,"
    r"[\bar{v}(p_2)\gamma^{\mu}u(p_1)]\,[\bar{u}(p_3)\gamma_{\mu}v(p_4)]$"
)


def kinematics_cm(E_cm: float, theta: float, m_e: float, m_mu: float):
    """
    Center-of-mass kinematics for e- e+ -> mu- mu+.

    Beam axis = z; outgoing mu- in xz plane at angle theta from z.

    Returns (p1, p2, p3, p4).
    """
    if E_cm < 2 * max(m_e, m_mu):
        raise ValueError(f"E_cm = {E_cm} below 2*max(m_e, m_mu)")
    E = E_cm / 2.0
    p_e = float(np.sqrt(E ** 2 - m_e ** 2))
    p_mu = float(np.sqrt(E ** 2 - m_mu ** 2))

    p1 = np.array([E, 0.0, 0.0, p_e])           # e-
    p2 = np.array([E, 0.0, 0.0, -p_e])          # e+
    p3 = np.array([E, p_mu * np.sin(theta), 0.0, p_mu * np.cos(theta)])    # mu-
    p4 = np.array([E, -p_mu * np.sin(theta), 0.0, -p_mu * np.cos(theta)])  # mu+
    return p1, p2, p3, p4


def amplitude(p1, p2, p3, p4,
              u_e: np.ndarray, v_eb: np.ndarray,
              u_mu: np.ndarray, v_mub: np.ndarray,
              gammas, e: float = ELECTRIC_CHARGE) -> complex:
    """
    Single-spin substrate-algebra amplitude  i M  for e- e+ -> mu- mu+.

    Implements
      i M = (i e^2 / s)  [v_bar(p2) gamma^mu u(p1)]  [u_bar(p3) gamma_mu v(p4)]
    with mostly-minus metric for the gamma_mu contraction.
    """
    p_total = p1 + p2
    s = float(p_total[0] ** 2 - p_total[1] ** 2 - p_total[2] ** 2 - p_total[3] ** 2)
    if abs(s) < 1e-15:
        raise ValueError("On-shell s-channel pole; choose nonzero s.")

    g0 = gammas[0]
    vbar_eb = adjoint_spinor(v_eb, g0)
    ubar_mu = adjoint_spinor(u_mu, g0)

    M = 0.0 + 0.0j
    for mu in range(4):
        sign = 1.0 if mu == 0 else -1.0   # mostly-minus
        J_e_mu = complex(vbar_eb @ gammas[mu] @ u_e)
        J_mu_mu = complex(ubar_mu @ gammas[mu] @ v_mub)
        M += sign * J_e_mu * J_mu_mu

    return 1j * (e ** 2 / s) * M


def M_squared_avg(E_cm: float, theta: float, m_e: float, m_mu: float,
                  gammas, e: float = ELECTRIC_CHARGE) -> tuple:
    """
    Spin-averaged |M|^2 for e- e+ -> mu- mu+ at scattering angle theta in CM.

    Convention matches Compton: sum over substrate-spin states (4 per
    fermion line), divide by 16 (initial states), multiply by the
    relativistic-normalization factor (2 E_e)(2 E_eb)(2 E_mu)(2 E_mub).

    Returns (|M|^2_avg, p1, p2, p3, p4).
    """
    p1, p2, p3, p4 = kinematics_cm(E_cm, theta, m_e, m_mu)

    u_es = positive_energy_spinors(p1, gammas, m_e)
    v_ebs = negative_energy_spinors(p2, gammas, m_e)
    u_mus = positive_energy_spinors(p3, gammas, m_mu)
    v_mubs = negative_energy_spinors(p4, gammas, m_mu)

    M_sq_total = 0.0
    for u_e in u_es:
        for v_eb in v_ebs:
            for u_mu in u_mus:
                for v_mub in v_mubs:
                    M = amplitude(p1, p2, p3, p4,
                                  u_e, v_eb, u_mu, v_mub, gammas, e)
                    M_sq_total += abs(M) ** 2

    # Initial-state average:  4 (u_e) * 4 (v_eb) = 16 states
    N_in = len(u_es) * len(v_ebs)
    avg = M_sq_total / N_in

    # Relativistic normalization:  product of (2 E_i) over external fermions
    rel_norm = (2.0 * float(p1[0])) * (2.0 * float(p2[0])) \
               * (2.0 * float(p3[0])) * (2.0 * float(p4[0]))

    return rel_norm * avg, p1, p2, p3, p4


def textbook_M_squared(E_cm: float, theta: float, m_mu: float,
                       e: float = ELECTRIC_CHARGE) -> float:
    """
    Standard QED textbook (Peskin eq. 5.13) spin-averaged |M|^2 with m_e = 0.

      |M|^2_avg = (e^4) * [(1 + 4 m_mu^2/s) + (1 - 4 m_mu^2/s) cos^2(theta)]
    """
    s = E_cm ** 2
    r = 4.0 * (m_mu ** 2) / s
    cos2 = float(np.cos(theta)) ** 2
    return (e ** 4) * ((1.0 + r) + (1.0 - r) * cos2)
