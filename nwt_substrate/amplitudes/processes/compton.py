"""
Compton scattering:  e- + gamma  ->  e- + gamma.

Substrate-algebra expression (s-channel + u-channel):

    iM = u_bar(p') [ (-ie eps'_{nu} gamma^{nu}) S_F(p+k) (-ie eps_{mu} gamma^{mu})
                   + (-ie eps_{mu} gamma^{mu}) S_F(p-k') (-ie eps'_{nu} gamma^{nu}) ] u(p)

       = -i e^2 u_bar(p') [ slash{eps'*} S_F(p+k) slash{eps}
                          + slash{eps}  S_F(p-k') slash{eps'*} ] u(p)

where:
  p, p'    : incoming, outgoing electron 4-momenta
  k, k'    : incoming, outgoing photon 4-momenta
  eps, eps': incoming, outgoing photon polarization 4-vectors
  S_F(P)   : fermion propagator (slash{P} + m)/(P^2 - m^2)

In the substrate algebra, the spinors are 8-dimensional and there are
4 positive-energy states per momentum (2 physical spin x 2 internal SU(2)).
The spin-averaged |M|^2 sum reproduces Klein-Nishina (verified to 4.8e-16
in walk-phase 3).
"""

from __future__ import annotations

import numpy as np

from ..spinors import positive_energy_spinors, adjoint_spinor
from ..polarizations import transverse_polarizations
from ..propagators import fermion_propagator
from ..vertices import ELECTRIC_CHARGE
from ...algebra.dirac import slash_p


# LaTeX expression for the substrate-algebra amplitude (for diagram captions)
LATEX_EXPRESSION = (
    r"$i\mathcal{M} = -ie^2\,\bar{u}(p')\left["
    r"\,\not{\!\epsilon}'^{*} S_F(p+k)\,\not{\!\epsilon}"
    r"\;+\;\not{\!\epsilon}\, S_F(p-k')\,\not{\!\epsilon}'^{*}"
    r"\right]u(p)$"
)


def kinematics(omega: float, theta: float, m: float):
    """
    Lab-frame kinematics for e- at rest with incoming photon along +z.

    Returns (p, k, p_prime, k_prime, omega_prime).
    """
    p = np.array([m, 0.0, 0.0, 0.0])
    k = np.array([omega, 0.0, 0.0, omega])
    omega_prime = omega / (1.0 + (omega / m) * (1.0 - np.cos(theta)))
    k_prime = np.array([omega_prime,
                        omega_prime * np.sin(theta), 0.0,
                        omega_prime * np.cos(theta)])
    p_prime = p + k - k_prime
    return p, k, p_prime, k_prime, omega_prime


def amplitude(p, k, p_prime, k_prime, eps_in, eps_out_conj,
              u_in, u_out, gammas, m, e=ELECTRIC_CHARGE):
    """
    Single-spin / single-polarization Compton amplitude (complex scalar).

    Implements the s-channel + u-channel formula above, returning iM
    (without the -ie^2 factor pulled into the function as 'e').
    """
    g0 = gammas[0]
    u_bar = adjoint_spinor(u_out, g0)

    # slash(eps) for input/output polarizations
    sl_eps = sum(gammas[i] * eps_in[i] * (1 if i == 0 else -1) for i in range(4))
    sl_eps_p = sum(gammas[i] * eps_out_conj[i] * (1 if i == 0 else -1)
                   for i in range(4))

    # Internal fermion propagators
    S_s = fermion_propagator(p + k, gammas, m)
    S_u = fermion_propagator(p - k_prime, gammas, m)

    Gamma_s = sl_eps_p @ S_s @ sl_eps
    Gamma_u = sl_eps @ S_u @ sl_eps_p

    M = u_bar @ (Gamma_s + Gamma_u) @ u_in
    # Multiply by the (-i e^2) overall factor:
    return -1j * (e ** 2) * complex(M)


def M_squared_avg(omega: float, theta: float, m: float, gammas,
                  e: float = ELECTRIC_CHARGE) -> tuple:
    """
    Spin/polarization-averaged |M|^2 for Compton at angle theta, in QED's
    relativistic normalization.

    Returns (|M|^2_avg, omega_prime).

    The substrate uses unit-norm Hilbert-space spinors (||u||^2 = 1), while
    QED uses relativistic normalization (u_bar u = 2m, equivalent to
    ||u||^2 = 2E).  We absorb the conversion factor (2E_in)(2E_out) here
    so the result is directly comparable to standard textbook formulas.

    Average over initial spins (4 substrate states) and polarizations (2),
    sum over final spins (4) and polarizations (2).
    """
    p, k, p_prime, k_prime, omega_prime = kinematics(omega, theta, m)

    u_ins = positive_energy_spinors(p, gammas, m)
    u_outs = positive_energy_spinors(p_prime, gammas, m)
    eps_ins = transverse_polarizations(k)
    eps_outs = transverse_polarizations(k_prime)

    M_sq_total = 0.0
    for u_in in u_ins:
        for u_out in u_outs:
            for eps_in in eps_ins:
                for eps_out in eps_outs:
                    M = amplitude(p, k, p_prime, k_prime,
                                  eps_in, np.conjugate(eps_out),
                                  u_in, u_out, gammas, m, e)
                    M_sq_total += abs(M) ** 2

    N_in = len(u_ins) * len(eps_ins)   # = 4 spin x 2 pol = 8

    # Convert unit-norm to QED relativistic normalization
    E_in = float(p[0])
    E_out = float(p_prime[0])
    rel_norm = (2.0 * E_in) * (2.0 * E_out)

    return rel_norm * M_sq_total / N_in, omega_prime


def klein_nishina_M_squared(omega: float, theta: float, m: float,
                            e: float = ELECTRIC_CHARGE) -> tuple:
    """
    Standard Klein-Nishina |M|^2 (averaged) for comparison.

    |M|^2_KN = 2 e^4 [omega'/omega + omega/omega' - sin^2(theta)]

    Returns (|M|^2_KN, omega_prime).
    """
    omega_prime = omega / (1.0 + (omega / m) * (1.0 - np.cos(theta)))
    sin2 = np.sin(theta) ** 2
    M_sq = 2.0 * (e ** 4) * (omega_prime / omega + omega / omega_prime - sin2)
    return M_sq, omega_prime
