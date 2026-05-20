"""Abelian-Higgs condensate at the BPS critical point.

Paper 16 NWT Lagrangian (minimal):

    L_NWT = |D_μ ψ|²  -  (1/4) F_μν F^μν  -  (λ/4)(|ψ|² - v²)²

at the BPS critical coupling λ = e²/2 (Bogomolny 1976), with the
condensate identification v = m_e (Paper 16 §L_2).  Vortex cores are
zeros of ψ; far from cores, |ψ|² → v² = m_e².

Mass spectrum of small perturbations around ψ_0 = v
-----------------------------------------------------
Write ψ = (v + σ + iτ)/√2 with σ, τ real.  Substitute into L_NWT,
keep quadratic terms.  After the gauge field eats τ (Higgs mechanism)
the spectrum is:

    Higgs (radial) mode σ:    m_σ² = λ v² = (e²/2) m_e²  at BPS
    Gauge boson A_μ:           m_A² = e² v² = e² m_e²      (Higgs mech.)

Bogoliubov-like dispersion is the relativistic massive form:

    E_k² = ℏ² c² k²  +  m_σ² c⁴      (Higgs branch)
    E_k² = ℏ² c² k²  +  m_A² c⁴      (gauge branch)

NO phonon limit at low k (relativistic theory). The relevant "healing
length" is the condensate correlation length ξ = ℏ/(m_σ c) (Higgs
Compton wavelength) — sets the vortex-core radius.

Substrate-monism anchors (Paper 16 §L_2):
  v = m_e                                   condensate VEV identification
  e = √(4π α_NWT)                            U(1) gauge coupling, Paper 17 α
  λ = e²/2  =  2π α_NWT                      BPS critical coupling
  μ_BPS = 2π v² c³/ℏ = 2π m_e² c³/ℏ          line tension (Paper 6 mass scale)

Dimensional anchors derived:
  m_σ = √(λ) v = √(2π α) m_e ≈ 0.214 m_e     Higgs mode (substrate prediction!)
  m_A = e v   = √(4π α) m_e ≈ 0.303 m_e     gauge boson (substrate prediction!)
  ξ_σ = ℏ/(m_σ c) = λ̄_C / √(2π α)            Higgs correlation length
  ξ_A = ℏ/(m_A c) = λ̄_C / √(4π α)            gauge correlation length

These predict NEW substrate species at sub-electron mass scale.  Whether
m_σ and m_A correspond to known particles or to substrate-specific
modes is a Phase D question (σ-orbit projection).

References: Paper 16 §L_2 (BPS sector); Paper 6 (mass spectrum, line
tension); [[framework_bridge_density_and_condensate]];
[[framework_healing_length_principle]] (ξ_substrate = λ̄_C as foundational).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from nwt_substrate.gravity.constants import (
    HBAR_J_S as HBAR_SI,
    C_LIGHT_M_S as C_LIGHT_SI,
    M_ELECTRON_KG as M_E_SI,
)
from nwt_substrate.isa.constants import ALPHA_NWT

# Natural-units conversion factors
M_E_EV = 0.510998928e6                       # electron mass, eV
HBARC_EV_NM = 197.3269804                    # ℏc, eV·nm
HBARC_EV_M = HBARC_EV_NM * 1e-9              # ℏc, eV·m  ≈ 1.973e-7
HBARC_EV_FM = HBARC_EV_NM * 1e6              # ℏc, eV·fm ≈ 1.973e8
LAMBDA_C_M = HBARC_EV_M / M_E_EV             # λ̄_C = ℏ/(m_e c), m ≈ 3.86e-13
LAMBDA_C_FM = LAMBDA_C_M * 1e15              # in fm ≈ 386 fm


@dataclass(frozen=True)
class AbelianHiggsParams:
    """Parameters of the Paper 16 abelian-Higgs condensate at the BPS point.

    All quantities in NATURAL UNITS (eV, fm; ℏ = c = 1) unless flagged.

    Attributes
    ----------
    v_phi_eV : float
        Condensate VEV v, in eV (so ψ has mass dimension 1).
        Substrate identification: v = m_e = 5.11×10⁵ eV.
    e_gauge : float
        U(1) gauge coupling.  Substrate value: e = √(4πα_NWT).
    lambda_quartic : float
        λ in (λ/4)(|ψ|²-v²)² quartic.  At BPS, λ = e²/2.

    Derived properties (in natural units):
    --------------------------------------
    m_sigma_eV : Higgs mode mass = √λ · v
    m_A_eV     : gauge boson mass = e · v   (Higgs mechanism)
    xi_sigma_m : Higgs correlation length = ℏ/(m_σ c)
    xi_A_m     : gauge correlation length = ℏ/(m_A c)
    mu_BPS_eV_per_fm : BPS line tension μ = 2π v² (in units where ℏc = 1)
                       reproduces Paper 16 §L_2 result 8.31 keV/fm at v = m_e.
    """

    v_phi_eV: float
    e_gauge: float
    lambda_quartic: float | None = None  # default to BPS critical value

    def __post_init__(self) -> None:
        if self.lambda_quartic is None:
            object.__setattr__(self, "lambda_quartic",
                                 self.e_gauge ** 2 / 2.0)

    @property
    def at_bps(self) -> bool:
        return abs(self.lambda_quartic - self.e_gauge ** 2 / 2.0) < 1e-12

    @property
    def m_sigma_eV(self) -> float:
        """Higgs mode mass m_σ = √λ · v."""
        return math.sqrt(self.lambda_quartic) * self.v_phi_eV

    @property
    def m_A_eV(self) -> float:
        """Gauge boson mass m_A = e · v (Higgs mechanism)."""
        return self.e_gauge * self.v_phi_eV

    @property
    def xi_sigma_m(self) -> float:
        """Higgs correlation length ξ_σ = ℏ/(m_σ c)."""
        return HBARC_EV_M / self.m_sigma_eV

    @property
    def xi_A_m(self) -> float:
        """Gauge correlation length ξ_A = ℏ/(m_A c)."""
        return HBARC_EV_M / self.m_A_eV

    @property
    def mu_BPS_eV_per_fm(self) -> float:
        """BPS line tension μ_BPS = 2π v²/(ℏc) (in eV/fm).

        In natural units (ℏc = 1) μ = 2π v² has units of (energy)².
        Restoring units: μ [eV/fm] = 2π v² [eV²] / (ℏc) [eV·fm].
        """
        return 2.0 * math.pi * self.v_phi_eV ** 2 / HBARC_EV_FM

    @classmethod
    def substrate_natural(cls) -> "AbelianHiggsParams":
        """Paper 16 §L_2 substrate identification: v = m_e, e = √(4πα),
        λ = e²/2 (BPS).
        """
        e_natural = math.sqrt(4.0 * math.pi * ALPHA_NWT)
        return cls(v_phi_eV=M_E_EV, e_gauge=e_natural)


def bogoliubov_dispersion_higgs(k_per_m: np.ndarray | float,
                                 p: AbelianHiggsParams) -> np.ndarray | float:
    """Relativistic dispersion E_k for Higgs-branch perturbations.

    E_k² = (ℏ c k)² + m_σ² c⁴

    Returns E in eV.  No phonon limit (relativistic).  At k → 0:
    E_k → m_σ c² (gapped). At k → ∞: E_k → ℏ c k (photon-like).
    """
    k = np.asarray(k_per_m)
    kinetic = (HBARC_EV_M * k) ** 2
    return np.sqrt(kinetic + p.m_sigma_eV ** 2)


def bogoliubov_dispersion_gauge(k_per_m: np.ndarray | float,
                                 p: AbelianHiggsParams) -> np.ndarray | float:
    """Relativistic dispersion E_k for the gauge-boson branch.

    E_k² = (ℏ c k)² + m_A² c⁴

    At k → 0: E_k → m_A c². At k → ∞: E_k → ℏ c k.
    """
    k = np.asarray(k_per_m)
    kinetic = (HBARC_EV_M * k) ** 2
    return np.sqrt(kinetic + p.m_A_eV ** 2)


def sound_speed_relativistic() -> float:
    """For a relativistic abelian-Higgs theory, both branches asymptote
    to c at high k.  No distinct sound speed — perturbations propagate
    at c modulo gap.  Returns 1.0 (in c units)."""
    return 1.0


def healing_length_higgs(p: AbelianHiggsParams) -> float:
    """ξ_σ = ℏ/(m_σ c) — Higgs Compton wavelength (vortex core radius)."""
    return p.xi_sigma_m


def healing_length_gauge(p: AbelianHiggsParams) -> float:
    """ξ_A = ℏ/(m_A c) — gauge correlation length (magnetic penetration)."""
    return p.xi_A_m


def line_tension_BPS(p: AbelianHiggsParams) -> float:
    """μ_BPS = 2π v² in eV/fm.  Paper 16 §L_2 + Paper 6 line tension."""
    return p.mu_BPS_eV_per_fm


# Backwards-compatible aliases (Phase A names)
def sound_speed(p: AbelianHiggsParams) -> float:
    return sound_speed_relativistic()


def healing_length(p: AbelianHiggsParams) -> float:
    """Default healing length = Higgs correlation length ξ_σ."""
    return healing_length_higgs(p)


def bogoliubov_dispersion(k: np.ndarray | float, p: AbelianHiggsParams,
                           dimensionless: bool = False) -> np.ndarray | float:
    """Default dispersion = Higgs branch (relativistic massive)."""
    if dimensionless:
        # Express k in units of 1/ξ_σ, return E in units of m_σ c²
        k_arr = np.asarray(k)
        return np.sqrt(k_arr ** 2 + 1.0)
    return bogoliubov_dispersion_higgs(k, p)
