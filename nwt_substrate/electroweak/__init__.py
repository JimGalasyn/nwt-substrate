"""
nwt.electroweak
===============

The electroweak view of the substrate algebra: SU(2)_L × U(1)_Y broken
to U(1)_em via the Higgs mechanism, with the Z boson as the chiral
linear combination W₃ cos θ_W − B sin θ_W.

The substrate connection (and an open thread): Cl(0,7)'s internal SU(2)
bivector triplet is *vector* under Lorentz (not chiral), while the SM's
SU(2)_L is chiral.  The library handles the chiral structure at the
Lagrangian / vertex-factor level via V-A couplings (memory:
substrate-boson-count-so7, Walk-phase 4a).

This shim provides:
  - PDG constants (M_Z, Γ_Z, M_W, Γ_W, sin²θ_W, ...)
  - SM fermion couplings g_V, g_A and a coupling-table API
  - Z partial / total decay widths Γ(Z → f f̄)
  - σ(e⁺e⁻ → f f̄) via γ + Z + interference (closed-form, agrees with
    Phase Q.8 R-ratio / Q.9 DY in the photon-only limit)

Quick start::

    import nwt_substrate.electroweak as ew

    # Z couplings
    print(ew.coupling("u"))           # T_3=+0.5, Q=+0.667, g_V=+0.192, ...
    print(ew.coupling_summary())      # full SM coupling table

    # Z boson properties
    print(ew.M_Z, ew.GAMMA_Z, ew.SIN2_THETA_W)

    # Z partial widths (sum to total Γ_Z within ~1% of PDG)
    print(ew.width_summary())

    # Cross-sections
    sigma_mumu = ew.sigma_total(91.2, "mu")     # peak at Z pole
    sigma_qed_only = ew.sigma_qed_only(91.2, "mu")
    print(f"At Z pole: σ(γ+Z) = {sigma_mumu:.3e} pb,  σ(γ-only) = {sigma_qed_only:.3e}")
    # → Z resonance enhances σ by ~1000x
"""

from .constants import (
    M_Z, GAMMA_Z, M_W, GAMMA_W,
    SIN2_THETA_W, COS2_THETA_W,
    ALPHA_QED, E_CHARGE, E_CHARGE_Z, G_W, G_W_SQ, G_Z, G_Z_SQ, G_F_GEV,
    V_HIGGS_GEV, M_HIGGS,
)

from .couplings import (
    WeakCoupling,
    SM_COUPLINGS,
    coupling,
    coupling_summary,
)

from .decays import (
    FERMION_MASS_GEV,
    partial_width_Z,
    total_width_Z,
    branching_ratios_Z,
    width_summary,
)

from .process import (
    chi,
    M_squared_avg,
    dsigma_dcos,
    sigma_total,
    sigma_qed_only,
)


__all__ = [
    # constants
    "M_Z", "GAMMA_Z", "M_W", "GAMMA_W",
    "SIN2_THETA_W", "COS2_THETA_W",
    "ALPHA_QED", "E_CHARGE", "E_CHARGE_Z",
    "G_W", "G_W_SQ", "G_Z", "G_Z_SQ", "G_F_GEV",
    "V_HIGGS_GEV", "M_HIGGS",
    # couplings
    "WeakCoupling", "SM_COUPLINGS", "coupling", "coupling_summary",
    # decays
    "FERMION_MASS_GEV", "partial_width_Z", "total_width_Z",
    "branching_ratios_Z", "width_summary",
    # process
    "chi", "M_squared_avg", "dsigma_dcos", "sigma_total", "sigma_qed_only",
]
