"""Non-relativistic superfluid limit of the abelian-Higgs condensate.

Reconciles Paper 16 (relativistic abelian Higgs) with Paper 5 (Gross-
Pitaevskii / superfluid vortex picture).  Phase C established that
the framework's ξ_substrate = λ̄_C is NOT a linear projection invariant
of the relativistic Bogoliubov spectrum.  This module shows where
ξ = λ̄_C *does* come from: it is the **Lorentz-invariance constraint**
on the NR limit of the abelian-Higgs condensate.

Paper 5 §II derivation
----------------------
GP equation:    iℏ ∂_t ψ = -ℏ²/(2m*) ∇²ψ  +  g|ψ|² ψ
with effective boson mass m* = m_e/√2 and uniform-condensate density
n_0 = ρ_0/m*.  The healing length is:

    ξ = ℏ / √(2 m* g n_0)                          (1)

The substrate identification ξ = λ̄_C = ℏ/(m_e c) **requires**:

    2 m* g n_0 = m_e² c²                            (2)

With m* = m_e/√2, equation (2) gives:

    g n_0 = m_e c² / √2                             (3)

which is the chemical potential μ_chem of the condensate.  The Bogoliubov
sound speed is:

    c_s = √(g n_0 / m*) = √((m_e c²/√2) / (m_e/√2)) = c   (4)

So ξ = λ̄_C is enforced by c_s = c (Lorentz invariance at leading order
in the NR limit).  These are NOT independent — picking either fixes both.

Phase B's relativistic gap masses m_σ = √(2πα)·m_e and m_A = √(4πα)·m_e
sit ABOVE the NR limit's effective particle mass m* = m_e/√2 ≈ 0.707 m_e:
they are massive POLE excitations of the field theory, not the gapless
phonon branch of the superfluid.  The "correlation lengths" ξ_σ ≈ 4.67 λ̄_C
and ξ_A ≈ 3.30 λ̄_C are the Compton wavelengths of those gap modes, NOT
the superfluid healing length.

Two distinct length scales coexist:
  ξ_substrate = λ̄_C            superfluid healing (NR limit, c_s = c)
  ξ_σ ~ 1/m_σ ~ 4.67 λ̄_C       Higgs Compton (relativistic gap mode)
  ξ_A ~ 1/m_A ~ 3.30 λ̄_C       gauge Compton (Higgs-mechanism mass)

The framework's ξ_substrate is the FIRST.  Phase B/C computed the SECOND
+ THIRD.  No contradiction — just different observables.

For σ-orbit bound states (electron, proton, ...): Paper 5 uses vortex-
ring energy with core cutoff ξ.  The (p, q) torus-knot solitons of
Papers 11-13 are the bound-state mass eigenvalues.  This module supplies
the underlying GP framework; the σ-orbit → torus-knot map is in
[[framework_k7_heffter_embedding]] and Papers 11-13.

References: Paper 5 §II eqs.(1)-(3); Paper 16 §L_1 + §L_2; Paper 11
abelian-Higgs ↔ SM connection.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from nwt_substrate.condensate.abelian_higgs import (
    AbelianHiggsParams, HBARC_EV_M, HBARC_EV_FM, M_E_EV, LAMBDA_C_M,
    LAMBDA_C_FM,
)
from nwt_substrate.gravity.constants import (
    HBAR_J_S as HBAR_SI,
    C_LIGHT_M_S as C_LIGHT_SI,
    M_ELECTRON_KG as M_E_KG,
)
from nwt_substrate.isa.constants import ALPHA_NWT


@dataclass(frozen=True)
class NRSuperfluidParams:
    """NR superfluid parameters consistent with Paper 5 §II.

    Determined by THREE substrate identifications:
      (1) v = m_e (relativistic VEV, Paper 16)
      (2) ξ = λ̄_C  (healing length = electron Compton wavelength)
      (3) c_s = c  (sound speed = c, Lorentz invariance)

    Conditions (2) and (3) are equivalent given (1) and m* = m_e/√2.

    All quantities in SI.

    Attributes
    ----------
    m_star_kg : effective boson mass m* = m_e/√2
    g_SI : self-interaction coupling, J·m³ (depends on n_0 calibration)
    n_0_per_m3 : equilibrium condensate number density, 1/m³
    rho_0_kg_m3 : mass density ρ_0 = m* n_0, kg/m³
    """

    m_star_kg: float
    g_SI: float
    n_0_per_m3: float

    @property
    def rho_0_kg_m3(self) -> float:
        return self.m_star_kg * self.n_0_per_m3

    @property
    def mu_chem_J(self) -> float:
        """Chemical potential μ_chem = g n_0."""
        return self.g_SI * self.n_0_per_m3

    @property
    def healing_length_m(self) -> float:
        """ξ = ℏ / √(2 m* g n_0).  Should equal λ̄_C for substrate-natural."""
        return HBAR_SI / math.sqrt(2.0 * self.m_star_kg * self.mu_chem_J)

    @property
    def sound_speed_m_s(self) -> float:
        """c_s = √(g n_0 / m*). Should equal c for substrate-natural."""
        return math.sqrt(self.mu_chem_J / self.m_star_kg)

    @property
    def vortex_ring_energy_J(self) -> float:
        """Energy of a single straight quantum vortex per unit length,
        cut off at the healing length.

        E/L = (π ρ_0 / m*²) × (ℏ²) × ln(R/ξ)  for a vortex ring of size R.
        For R ~ ξ × π² (Derrick equilibrium, Paper 16): ln ~ ln(π²) ≈ 2.29.

        Returns energy per unit length (J/m).  Phase D refinement.
        """
        # Placeholder: full formula in Paper 5 §III
        return math.nan


def substrate_natural_NR() -> NRSuperfluidParams:
    """Build NR superfluid parameters from substrate identifications:

        m* = m_e/√2
        g n_0 = m_e c² / √2  (forces ξ = λ̄_C AND c_s = c simultaneously)

    The free parameter is the calibration ratio between g and n_0; Paper 5
    uses ρ_0 = 1.69×10⁵ kg/m³ (calibrated from m_e c² vortex-ring matching
    in Paper 4).  We adopt that here.
    """
    m_star = M_E_KG / math.sqrt(2.0)
    rho_0 = 1.69e5  # kg/m³, Paper 5 calibration
    n_0 = rho_0 / m_star
    mu_chem = M_E_KG * C_LIGHT_SI ** 2 / math.sqrt(2.0)  # = g n_0
    g_SI = mu_chem / n_0
    return NRSuperfluidParams(
        m_star_kg=m_star,
        g_SI=g_SI,
        n_0_per_m3=n_0,
    )


def verify_substrate_natural() -> dict:
    """Check Paper 5 §II self-consistency:
       ξ = λ̄_C  AND  c_s = c  must both hold to machine precision.
    """
    p = substrate_natural_NR()
    xi_ratio = p.healing_length_m / LAMBDA_C_M
    cs_ratio = p.sound_speed_m_s / C_LIGHT_SI
    return {
        "params": p,
        "xi_m": p.healing_length_m,
        "lambda_C_m": LAMBDA_C_M,
        "xi_over_lambda_C": xi_ratio,
        "xi_match": abs(xi_ratio - 1.0) < 1e-6,
        "c_s_m_s": p.sound_speed_m_s,
        "c_m_s": C_LIGHT_SI,
        "c_s_over_c": cs_ratio,
        "c_s_match": abs(cs_ratio - 1.0) < 1e-6,
        "mu_chem_eV": p.mu_chem_J / 1.602176634e-19,
        "rho_0_kg_m3": p.rho_0_kg_m3,
        "n_0_per_m3": p.n_0_per_m3,
        "g_SI": p.g_SI,
    }
