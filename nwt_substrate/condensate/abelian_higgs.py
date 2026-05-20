"""Abelian-Higgs condensate at the BPS critical point.

Paper 16 NWT Lagrangian (minimal):

  L_NWT = |D_μ ψ|²  -  (1/4) F_μν F^μν  -  (λ/4)(|ψ|² - v²)²

at the BPS critical coupling λ = e²/2 (Bogomolny 1976). The condensate
order parameter ψ has vacuum expectation value v.  Vortex cores are
zeros of ψ; far from cores, |ψ|² → v².

For the slow-cosmogenesis Bogoliubov derivation we need:
  - The healing length ξ that sets the vortex core size
  - The sound speed c_s of Bogoliubov modes
  - The Bogoliubov dispersion E_k for perturbations around ψ₀ = v
  - The line tension μ_BPS = 2π v² (Paper 6 electron mass via Derrick scaling)

In dimensionless GPE units (length = ξ, time = ξ/c_s, ψ = ψ/v), the
linearized perturbation equation has dispersion

  E_k²  =  ε_k (ε_k + 2)

with ε_k = k²/2.  This is the textbook scalar BEC Bogoliubov spectrum
recovered as a sanity check (Phase B target).

Substrate-monism anchors:
  - ξ_substrate = λ̄_C = ℏ/(m_e c)  →  Paper 17 trefoil α gives ψ-mass
  - μ_BPS = 2π v²  →  Paper 6 electron mass = 2π m_e²/(ℏ c) at v = m_e/√(2π)
  - λ = e²/2  →  BPS critical, fixed by Paper 17 α

[[framework_healing_length_principle]] holds ξ_substrate = λ̄_C as a
foundational substrate primitive; this module checks that the abelian
Higgs Lagrangian reproduces it.

References: Paper 16 §L_2 (BPS sector); Paper 6 (mass spectrum, line
tension); [[framework_bridge_density_and_condensate]].
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


@dataclass(frozen=True)
class AbelianHiggsParams:
    """Parameters of the Paper 16 abelian-Higgs condensate at the BPS point.

    All quantities in SI unless flagged.

    Attributes
    ----------
    v_phi : float
        Condensate VEV v (in field units chosen so |ψ_0| = v).
    e_gauge : float
        U(1) gauge coupling.
    m_psi : float
        Effective mass m_ψ of the condensate quanta, kg.  Sets ξ via
        ξ = ℏ / √(2 m_ψ λ v²).  In the substrate frame, m_ψ = m_e
        (electron mass) is the natural identification — vortex cores
        have core radius ξ = λ̄_C.
    lambda_quartic : float
        λ in (λ/4)(|ψ|²-v²)² quartic.  At BPS, λ = e²/2.
    """

    v_phi: float
    e_gauge: float
    m_psi: float = M_E_SI
    lambda_quartic: float | None = None  # default to BPS value e²/2

    def __post_init__(self) -> None:
        if self.lambda_quartic is None:
            object.__setattr__(self, "lambda_quartic",
                                 self.e_gauge ** 2 / 2.0)

    @property
    def at_bps(self) -> bool:
        return abs(self.lambda_quartic - self.e_gauge ** 2 / 2.0) < 1e-12

    @classmethod
    def substrate_natural(cls) -> "AbelianHiggsParams":
        """First-pass substrate-natural set: m_ψ = m_e, e = √(4πα).

        ★ Phase B TODO: the v_phi value here is a NAIVE SI placeholder
        and does NOT yet reproduce ξ = λ̄_C.  Dimensional anchoring
        requires extracting the field-units convention from Paper 16
        and matching the BPS line tension μ_BPS = 2π m_e²/(ℏc) =
        line-tension form (Eq. ~ Paper 16 §L_2 around line 333).

        The DIMENSIONLESS dispersion below is correct (phonon → particle
        crossover at k = 1/ξ); the dimensional ξ value is wrong by ~10¹¹
        until the Paper 16 field-units convention is locked in.
        """
        e_natural = math.sqrt(4 * math.pi * ALPHA_NWT)   # √(4πα) gauge coupling
        v_natural = M_E_SI / math.sqrt(HBAR_SI * C_LIGHT_SI)
        return cls(v_phi=v_natural, e_gauge=e_natural, m_psi=M_E_SI)


def sound_speed(p: AbelianHiggsParams) -> float:
    """Bogoliubov sound speed c_s = v · √(λ / m_ψ).

    In dimensionless GPE units c_s = 1; in physical units c_s sets the
    propagation speed of low-k perturbations.
    """
    return p.v_phi * math.sqrt(p.lambda_quartic / p.m_psi)


def healing_length(p: AbelianHiggsParams) -> float:
    """Healing length ξ = ℏ / √(2 m_ψ λ v²).

    Substrate-monism target: ξ = λ̄_C = ℏ / (m_ψ c).
    """
    return HBAR_SI / math.sqrt(2.0 * p.m_psi * p.lambda_quartic * p.v_phi ** 2)


def line_tension_BPS(p: AbelianHiggsParams) -> float:
    """Vortex line tension at BPS: μ_BPS = 2π v² (Paper 6 §III)."""
    return 2.0 * math.pi * p.v_phi ** 2


def bogoliubov_dispersion(k: np.ndarray | float,
                           p: AbelianHiggsParams,
                           dimensionless: bool = True) -> np.ndarray | float:
    """Bogoliubov dispersion E_k for perturbations δψ̂ around ψ_0 = v.

    Dimensionless form (k in 1/ξ, E in m_ψ c_s²):

        E_k²  =  (k² / 2) · (k² / 2  +  2)

    Physical form:

        E_k²  =  (ℏ²k² / 2m_ψ) · (ℏ²k² / 2m_ψ  +  2 λ v²)

    Low-k limit:   E_k  →  ℏ c_s k    (phonon)
    High-k limit:  E_k  →  ℏ²k²/2m_ψ + λv²    (free particle + mean-field shift)

    The crossover scale is k ~ 1/ξ (healing-length wavenumber).
    """
    k = np.asarray(k)
    if dimensionless:
        eps_k = k ** 2 / 2.0
        return np.sqrt(eps_k * (eps_k + 2.0))
    eps_k = (HBAR_SI ** 2 * k ** 2) / (2.0 * p.m_psi)
    mean_field = 2.0 * p.lambda_quartic * p.v_phi ** 2
    return np.sqrt(eps_k * (eps_k + mean_field))


def crossover_wavenumber(p: AbelianHiggsParams) -> float:
    """Wavenumber at which phonon→free-particle crossover happens, k_xover = 1/ξ."""
    return 1.0 / healing_length(p)
