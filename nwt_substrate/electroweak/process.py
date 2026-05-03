"""
e⁺ e⁻ → f f̄ via γ + Z exchange (full neutral-current Drell-Yan).

The combined γ + Z amplitude squared, spin- and color-averaged in the
massless-fermion limit, is

    |M|²_avg(e⁺e⁻ → f f̄)  =  e⁴ × N_c^f × [ A_γ + A_int + A_Z ]

with

    A_γ   = Q_e² Q_f² × (1 + cos²θ)

    A_int = 2 Q_e Q_f × Re[χ(s)]
              × [ (g_V^e g_V^f + g_A^e g_A^f) (1 + cos²θ)
                 + 2 (g_V^e g_A^f + g_A^e g_V^f) cos θ ]

    A_Z   = |χ(s)|² × [ ((g_V^e)² + (g_A^e)²) ((g_V^f)² + (g_A^f)²) (1 + cos²θ)
                       + 8 g_V^e g_A^e g_V^f g_A^f cos θ ]

where  χ(s) = s / (4 sin²θ_W cos²θ_W) × 1 / (s - M_Z² + i M_Z Γ_Z)
            = s / sin²(2θ_W) × 1 / (s - M_Z² + i M_Z Γ_Z).

Reference: Peskin-Schroeder Chapter 20, Halzen-Martin Eq 13.55-13.59,
Schwartz Chapter 29.

The integrated cross-section is obtained by ∫dΩ; the cos θ terms in
A_int and A_Z give the famous forward-backward asymmetry (the Z's
parity-violating signature) but vanish in the total cross-section.
"""

from __future__ import annotations
import math

from .constants import (M_Z, GAMMA_Z, SIN2_THETA_W, COS2_THETA_W,
                          ALPHA_QED, E_CHARGE_Z, G_Z_SQ)
from .couplings import coupling


GEV2_TO_PB = 0.3894e9   # (GeV)⁻² → pb

# Effective α at the Z pole, derived from the G_F-based EW coupling.
# Using α(0)=1/137 here would underestimate by ~6% (the QED running from
# 0 to M_Z is absorbed into G_F).  This is the standard EW convention.
ALPHA_Z = E_CHARGE_Z ** 2 / (4.0 * math.pi)


def chi(s: float) -> complex:
    """
    Z-propagator dressing factor (Halzen-Martin Eq 13.55):

        χ(s) = s / (4 sin²θ_W cos²θ_W) × 1 / (s - M_Z² + i M_Z Γ_Z)
             = s / sin²(2θ_W)            × 1 / (s - M_Z² + i M_Z Γ_Z)

    The (4 sin²θ_W cos²θ_W) factor comes from g_Z² / 4 / e², where
    g_Z = e / (sin θ_W cos θ_W) is the Z-fermion coupling and the (1/4)
    is the (1/2)² from Z vertex's chiral projector.
    """
    return (s / (4.0 * SIN2_THETA_W * COS2_THETA_W)) \
        / complex(s - M_Z ** 2, M_Z * GAMMA_Z)


def M_squared_avg(s: float, cos_theta: float, fermion: str) -> float:
    """
    |M|²(e⁺e⁻ → f f̄) summed over polarisations, averaged over initial
    spin states, including color (N_c × 3 for quark final states).
    Massless-fermion limit; for m_f << √s = √(s) the only mass effect
    is the kinematic β_f and (1 + 2m²/s) factor applied separately in
    the cross-section.

    Parameters
    ----------
    s : float                  -- E_cm² in GeV²
    cos_theta : float          -- scattering cosine in CM frame
    fermion : str              -- final-state fermion species
    """
    e_coupling = coupling("e")
    f_coupling = coupling(fermion)

    Q_e, gV_e, gA_e = e_coupling.Q, e_coupling.g_V, e_coupling.g_A
    Q_f, gV_f, gA_f = f_coupling.Q, f_coupling.g_V, f_coupling.g_A
    N_c = f_coupling.n_color

    chi_s = chi(s)
    chi_re = chi_s.real
    chi_abs2 = abs(chi_s) ** 2

    one_plus_cos2 = 1.0 + cos_theta ** 2
    cos = cos_theta

    A_gamma = (Q_e ** 2) * (Q_f ** 2) * one_plus_cos2

    A_int = (2.0 * Q_e * Q_f * chi_re
             * ((gV_e * gV_f + gA_e * gA_f) * one_plus_cos2
                + 2.0 * (gV_e * gA_f + gA_e * gV_f) * cos))

    A_Z = (chi_abs2
           * ((gV_e ** 2 + gA_e ** 2) * (gV_f ** 2 + gA_f ** 2) * one_plus_cos2
              + 8.0 * gV_e * gA_e * gV_f * gA_f * cos))

    e4 = (4.0 * math.pi * ALPHA_Z) ** 2  # G_F-derived EW coupling at Z scale
    return e4 * N_c * (A_gamma + A_int + A_Z)


def dsigma_dcos(s: float, cos_theta: float, fermion: str,
                 m_f: float = 0.0) -> float:
    """
    dσ/d(cos θ)  [pb]  for e⁺e⁻ → f f̄ via γ + Z exchange.

    Phase-space:
        dσ/dΩ = |M|²_avg / (64 π² s) × β_f
        dσ/d(cos θ) = 2π × dσ/dΩ
    """
    if 4.0 * m_f ** 2 >= s:
        return 0.0
    beta_f = math.sqrt(1.0 - 4.0 * m_f ** 2 / s)
    M_sq = M_squared_avg(s, cos_theta, fermion)
    return (M_sq / (32.0 * math.pi * s)) * beta_f * GEV2_TO_PB


def sigma_total(sqrt_s: float, fermion: str, m_f: float = 0.0) -> float:
    """
    σ_total(e⁺e⁻ → f f̄, √s)  [pb]  via γ + Z exchange.

    Integrating over cos θ ∈ [-1, +1]:

        ∫ (1 + cos²θ) d(cos θ) = 8/3
        ∫ cos θ d(cos θ)        = 0   (forward-backward asymmetric → 0 in total)

    so σ_total picks up only the (1 + cos²θ) terms × (8/3).
    """
    s = sqrt_s ** 2
    if 4.0 * m_f ** 2 >= s:
        return 0.0
    beta_f = math.sqrt(1.0 - 4.0 * m_f ** 2 / s)

    e_coupling = coupling("e")
    f_coupling = coupling(fermion)

    Q_e, gV_e, gA_e = e_coupling.Q, e_coupling.g_V, e_coupling.g_A
    Q_f, gV_f, gA_f = f_coupling.Q, f_coupling.g_V, f_coupling.g_A
    N_c = f_coupling.n_color

    chi_s = chi(s)
    chi_re = chi_s.real
    chi_abs2 = abs(chi_s) ** 2

    # Coefficients of (1 + cos²θ); the cos-only terms vanish in the total
    A_gamma = (Q_e ** 2) * (Q_f ** 2)
    A_int = 2.0 * Q_e * Q_f * chi_re * (gV_e * gV_f + gA_e * gA_f)
    A_Z = chi_abs2 * (gV_e ** 2 + gA_e ** 2) * (gV_f ** 2 + gA_f ** 2)

    e4 = (4.0 * math.pi * ALPHA_Z) ** 2  # G_F-derived EW coupling at Z scale
    M_sq_integrated = e4 * N_c * (A_gamma + A_int + A_Z) * (8.0 / 3.0)

    # σ = ∫ dσ/d(cos θ) d(cos θ); d(cos θ) factor is already in M_sq_integrated
    # via the (8/3) Jacobian.
    return (M_sq_integrated / (32.0 * math.pi * s)) * beta_f * GEV2_TO_PB


def sigma_qed_only(sqrt_s: float, fermion: str, m_f: float = 0.0) -> float:
    """
    σ(e⁺e⁻ → f f̄) using only the photon channel (no Z) -- for comparison
    with the full γ + Z result.  Equivalent to the QED prediction used by
    R-ratio (Phase Q.8) and Drell-Yan partonic (Phase Q.9).
    """
    s = sqrt_s ** 2
    if 4.0 * m_f ** 2 >= s:
        return 0.0
    f_coupling = coupling(fermion)
    Q_f, N_c = f_coupling.Q, f_coupling.n_color
    beta_f = math.sqrt(1.0 - 4.0 * m_f ** 2 / s)
    sigma = (4.0 * math.pi * ALPHA_QED ** 2 / (3.0 * s)) * (Q_f ** 2) * N_c * beta_f
    return sigma * GEV2_TO_PB
