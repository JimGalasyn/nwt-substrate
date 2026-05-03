"""
Cross-sections and decay rates from substrate-algebra |M|^2.

For 2->2 scattering in CM frame:

    dsigma/dOmega = (1/(64 pi^2 s)) * (|p_f|/|p_i|) * |M|^2_avg

    sigma_total = integral over solid angle

For 2->2 with one particle at rest (lab frame, like Compton):

    dsigma/dOmega_lab = (1/(64 pi^2 m_target^2)) * (omega_out/omega_in)^2 * |M|^2_avg

For 1->2 decay:

    Gamma = (1/(2 M)) * |M|^2_avg * |p_1| / (4 pi)        [2-body]

For 1->3 decay (e.g. muon decay):

    Gamma = (1/(256 pi^3 M^3)) * integral over Dalitz |M|^2_avg

Conventions: natural units (hbar = c = 1), GeV throughout.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.integrate import quad


def cm_momentum(s: float, m1: float, m2: float) -> float:
    """
    Magnitude of CM 3-momentum for two-body system with masses m1, m2 and
    invariant mass-squared s.

      |p|^2 = ( s - (m1+m2)^2 ) ( s - (m1-m2)^2 ) / (4 s)
    """
    if s < (m1 + m2) ** 2:
        return 0.0
    lam = (s - (m1 + m2) ** 2) * (s - (m1 - m2) ** 2)
    return float(np.sqrt(max(lam, 0.0))) / (2.0 * np.sqrt(s))


def sigma_2to2_cm(M_sq_func: Callable[[float], float],
                  s: float,
                  m_in: tuple[float, float],
                  m_out: tuple[float, float],
                  theta_min: float = 1e-4,
                  theta_max: float = np.pi - 1e-4) -> float:
    """
    Total 2->2 cross-section in the CM frame.

      sigma = (1/(32 pi s)) * (|p_f|/|p_i|) * integral_0^pi sin(theta) |M|^2(theta) d theta

    Args:
        M_sq_func   :  callable theta -> spin-averaged |M|^2 in QED rel-norm
        s           :  Mandelstam s = (p1 + p2)^2
        m_in        :  (m1, m2) incoming masses
        m_out       :  (m3, m4) outgoing masses
        theta_min, theta_max
                    :  cutoffs to avoid forward/backward singularities for
                       t-channel processes (e.g. Moller, Bhabha).

    Returns:
        cross-section in GeV^-2.  Multiply by (0.3894 mb / GeV^-2) to get mb.
    """
    p_i = cm_momentum(s, m_in[0], m_in[1])
    p_f = cm_momentum(s, m_out[0], m_out[1])
    if p_i == 0:
        return 0.0

    def integrand(theta):
        return M_sq_func(theta) * np.sin(theta)

    integral, _ = quad(integrand, theta_min, theta_max, limit=200)
    prefactor = (1.0 / (32.0 * np.pi * s)) * (p_f / p_i)
    return prefactor * integral


def dsigma_dOmega_cm(M_sq: float, s: float,
                     m_in: tuple[float, float],
                     m_out: tuple[float, float]) -> float:
    """Differential cross-section dsigma/dOmega in CM, GeV^-2 / sr."""
    p_i = cm_momentum(s, m_in[0], m_in[1])
    p_f = cm_momentum(s, m_out[0], m_out[1])
    if p_i == 0:
        return 0.0
    return M_sq * p_f / (p_i * 64.0 * np.pi ** 2 * s)


# ----------------------------------------------------------------------------
# Compton-style: lab frame, target at rest
# ----------------------------------------------------------------------------

def sigma_compton_lab(M_sq_at_theta: Callable[[float], tuple],
                      m_target: float,
                      theta_min: float = 1e-4,
                      theta_max: float = np.pi - 1e-4) -> float:
    """
    Total Compton cross-section in lab frame (target at rest).

      dsigma/dOmega_lab = (1/(64 pi^2 m^2)) * (omega'/omega)^2 * |M|^2_avg

    Args:
        M_sq_at_theta : callable theta -> (|M|^2_avg, omega_prime); takes
                        theta and returns the spin-averaged |M|^2 plus
                        the outgoing photon frequency.  Caller already
                        knows omega_in (lab-frame incoming photon energy).
        m_target      : target rest mass (e.g. m_e for Compton).

    Returns:
        sigma in GeV^-2.
    """
    # We need omega_in to form the (omega'/omega)^2 ratio.  Probe with theta=0.
    M_sq_0, omega_prime_0 = M_sq_at_theta(theta_min)
    # In the kinematics we use, omega_prime(theta=0) = omega_in.  So:
    omega_in = omega_prime_0

    def integrand(theta):
        M_sq, omega_prime = M_sq_at_theta(theta)
        return M_sq * (omega_prime / omega_in) ** 2 * np.sin(theta)

    integral, _ = quad(integrand, theta_min, theta_max, limit=200)
    return integral * 2.0 * np.pi / (64.0 * np.pi ** 2 * m_target ** 2)


def thomson_cross_section(m_target: float, alpha: float = 7.2973525693e-3) -> float:
    """
    Thomson cross-section (low-energy Compton limit):

        sigma_Thomson = (8 pi / 3) * r^2 = (8 pi / 3) * (alpha/m)^2

    in natural units (GeV^-2).
    """
    r = alpha / m_target
    return (8.0 * np.pi / 3.0) * r ** 2


# ----------------------------------------------------------------------------
# Standard reference cross-sections (textbook QED)
# ----------------------------------------------------------------------------

def sigma_eemumu_high_energy(s: float,
                             alpha: float = 7.2973525693e-3) -> float:
    """
    Textbook  sigma(e+e- -> mu+mu-)  in massless / high-energy limit:

        sigma = (4 pi alpha^2) / (3 s)

    in GeV^-2.
    """
    return 4.0 * np.pi * (alpha ** 2) / (3.0 * s)


def gev_inv_sq_to_mb(sigma_gev_inv_sq: float) -> float:
    """Convert GeV^-2 to millibarn (1 mb = 10^-27 cm^2)."""
    return sigma_gev_inv_sq * 0.3894


def gev_inv_sq_to_pb(sigma_gev_inv_sq: float) -> float:
    """Convert GeV^-2 to picobarn (1 pb = 10^-36 cm^2)."""
    return sigma_gev_inv_sq * 0.3894e9


# ----------------------------------------------------------------------------
# Decay rates (1 -> 2 and 1 -> 3 placeholders for Phase A.5)
# ----------------------------------------------------------------------------

def gamma_1to2(M_sq: float, M_parent: float,
               m1: float, m2: float) -> float:
    """
    Two-body decay rate:  Gamma = (1/(8 pi M)) * (|p|/M) * |M|^2_avg

    where |p| = cm_momentum(M^2, m1, m2).  Returns Gamma in GeV.
    """
    if M_parent < m1 + m2:
        return 0.0
    p = cm_momentum(M_parent ** 2, m1, m2)
    return M_sq * p / (8.0 * np.pi * M_parent ** 2)


# Standard reference:  muon decay rate
G_F = 1.1663787e-5   # Fermi constant, GeV^-2


def gamma_muon_decay_textbook(m_mu: float = 0.10566) -> float:
    """
    Standard textbook muon decay rate:

        Gamma_mu = G_F^2 * m_mu^5 / (192 pi^3)

    Lifetime: tau = 1/Gamma.  m_mu = 0.10566 GeV gives tau ~ 2.2 us.
    """
    return (G_F ** 2) * (m_mu ** 5) / (192.0 * np.pi ** 3)
