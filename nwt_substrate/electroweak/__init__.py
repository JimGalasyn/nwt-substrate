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

from .substrate_gf import (
    ALPHA_SUBSTRATE,
    M_E_GEV,
    V_EW_PDG_GEV,
    v_ew_substrate,
    fermi_constant_substrate,
    fermi_constant_from_vev,
    precision_chain,
    verify_substrate_gf,
    precision_chain_summary,
)

from .form_factors import (
    FORM_FACTORS,
    cos_theta_C_substrate,
    sin_theta_C_substrate,
    f_plus_substrate,
    f_plus_for,
    f_plus_leading_order,
    precision_chain as form_factor_precision_chain,
    verify_form_factors,
    precision_chain_summary as form_factor_precision_chain_summary,
)

# Meson decay constants live in `nwt_substrate.particles.decay_constants` — the
# single canonical home for the unified light + heavy + vector + B_c sectors.
# Import directly from there (no re-export here, to avoid a circular import
# with `electroweak.substrate_gf`).

from .substrate_ckm import (
    N_VERTICES_K7 as CKM_N_VERTICES_K7,
    Q_B as CKM_Q_B,
    DIM_G2 as CKM_DIM_G2,
    DIM_SO7 as CKM_DIM_SO7,
    DELTA_CP_CKM,
    wolfenstein_lambda,
    wolfenstein_lambda_sq,
    wolfenstein_A,
    wolfenstein_A_sq,
    wolfenstein_apex_magnitude,
    wolfenstein_apex_magnitude_sq,
    wolfenstein_rho_bar,
    wolfenstein_eta_bar,
    V_us as ckm_V_us,
    V_cb as ckm_V_cb,
    V_ub as ckm_V_ub,
    V_td as ckm_V_td,
    jarlskog_ckm,
    jarlskog_coefficient,
    ckm_matrix,
    precision_chain as ckm_precision_chain,
    verify_substrate_ckm,
    precision_chain_summary as ckm_precision_chain_summary,
)

# ---- Substrate-ISA identities ----
from ..isa import (
    N_GENERATIONS,
    N_FERMION_TYPES_PER_GENERATION,
    N_SM_FERMION_FLAVORS,
    N_LORENTZ_GENERATORS,
    B_QED_SM,
    DIM_OCTONION,
    DIM_INTERNAL_SU2,
    N_EDGES_K7,
    RANK_SO7,
    N_C_SU3,
)


def verify_b_qed_sm(tol: float = 1e-10) -> float:
    """Compute b_QED^SM from the SM_COUPLINGS table and assert it equals
    isa.B_QED_SM = DIM_OCTONION = 8.

    b_QED^SM = Σ over SM fermions of N_c × Q²

    For an "all SM fermions active" run, this equals exactly 8 (modulo
    float rounding). Returns the computed value; raises AssertionError
    if the substrate identity is violated.
    """
    # Charged leptons (e, mu, tau) + up-type (u, c, t) + down-type (d, s, b)
    # Neutrinos contribute 0 since Q = 0.
    charged = ["e", "mu", "tau"]
    up_type = ["u", "c", "t"]
    down_type = ["d", "s", "b"]
    b = 0.0
    for name in charged + up_type + down_type:
        c = SM_COUPLINGS[name]
        b += c.n_color * c.Q ** 2
    if abs(b - B_QED_SM) > tol:
        raise AssertionError(
            f"Substrate identity violated: b_QED^SM = {b} != B_QED_SM = {B_QED_SM}"
        )
    return b


def substrate_breakdown() -> str:
    """Pretty-print the EW / SM-fermion-content substrate decomposition."""
    lines = ["Electroweak fermion content from K_7 substrate:"]
    lines.append("")
    lines.append(f"    so(7) decomposition (21 generators total):")
    lines.append(f"      6  Lorentz so(1,3)         (N_LORENTZ_GENERATORS)")
    lines.append(f"      3  internal SU(2)           (DIM_INTERNAL_SU2)")
    lines.append(f"     12  mixed                    (N_SM_FERMION_FLAVORS)")
    lines.append(f"     {N_LORENTZ_GENERATORS} + {DIM_INTERNAL_SU2} + "
                 f"{N_SM_FERMION_FLAVORS} = {N_EDGES_K7} ✓")
    lines.append("")
    lines.append(f"    SM fermion organization:")
    lines.append(f"      N_GENERATIONS = RANK_SO7 = {N_GENERATIONS}")
    lines.append(f"      N_FERMION_TYPES_PER_GENERATION = {N_FERMION_TYPES_PER_GENERATION}")
    lines.append(f"        (1 charged lepton + 1 neutrino + 1 up + 1 down)")
    lines.append(f"      N_SM_FERMION_FLAVORS = {N_GENERATIONS} × "
                 f"{N_FERMION_TYPES_PER_GENERATION} = {N_SM_FERMION_FLAVORS}")
    lines.append(f"        = the 12 'mixed' so(7) generators (hypothesis)")
    lines.append("")
    lines.append(f"    QED 1-loop β-function at full SM activation:")
    lines.append(f"      b_QED^SM = Σ N_c × Q²")
    lines.append(f"               = 3 (leptons × 1²)")
    lines.append(f"               + 4 (up × 3 × (2/3)²)")
    lines.append(f"               + 1 (down × 3 × (1/3)²)")
    lines.append(f"               = {B_QED_SM}")
    lines.append(f"      SUBSTRATE IDENTITY: B_QED_SM = DIM_OCTONION = "
                 f"{DIM_OCTONION}")
    lines.append("")
    lines.append(f"    Cross-shim ID:")
    lines.append(f"      The same 8 = DIM_OCTONION that gives QED its 8×8 γ^μ")
    lines.append(f"      and QCD its 8 gluons also equals the SM fermion")
    lines.append(f"      charge-squared sum that drives QED running.")
    return "\n".join(lines)


class _SubstrateNamespace:
    """Cross-shim K_7 substrate identities visible from the EW shim."""
    N_GENERATIONS = N_GENERATIONS
    N_FERMION_TYPES_PER_GENERATION = N_FERMION_TYPES_PER_GENERATION
    N_SM_FERMION_FLAVORS = N_SM_FERMION_FLAVORS
    N_LORENTZ_GENERATORS = N_LORENTZ_GENERATORS
    B_QED_SM = B_QED_SM
    DIM_OCTONION = DIM_OCTONION
    DIM_INTERNAL_SU2 = DIM_INTERNAL_SU2
    N_EDGES_K7 = N_EDGES_K7
    RANK_SO7 = RANK_SO7
    N_C_SU3 = N_C_SU3
    # Substrate G_F derivation chain (Galasyn 2026-05-23 P7b §7.4)
    ALPHA_SUBSTRATE = ALPHA_SUBSTRATE
    fermi_constant = staticmethod(fermi_constant_substrate)
    v_ew = staticmethod(v_ew_substrate)
    g_f_precision_chain = staticmethod(precision_chain)
    verify_g_f = staticmethod(verify_substrate_gf)

    # Substrate Wolfenstein CKM (Galasyn 2026-05-23 P6b)
    DELTA_CP_CKM = DELTA_CP_CKM
    ckm = staticmethod(ckm_matrix)
    wolfenstein_lambda = staticmethod(wolfenstein_lambda)
    wolfenstein_A = staticmethod(wolfenstein_A)
    wolfenstein_apex_magnitude = staticmethod(wolfenstein_apex_magnitude)
    jarlskog_ckm = staticmethod(jarlskog_ckm)
    ckm_precision_chain = staticmethod(ckm_precision_chain)
    verify_ckm = staticmethod(verify_substrate_ckm)

    # Meson decay constants now live in nwt_substrate.particles.decay_constants
    # (P7b §2-3, §7.5, §7.6 — light + heavy + vector + B_c sectors).
    # Removed from this namespace to avoid a circular import.

    # Substrate Dalitz form factors f_+(0) (Galasyn 2026-05-23 P7b §7.7)
    cos_theta_C = staticmethod(cos_theta_C_substrate)
    f_plus = staticmethod(f_plus_substrate)
    f_plus_for = staticmethod(f_plus_for)
    form_factor_chain = staticmethod(form_factor_precision_chain)
    verify_form_factors = staticmethod(verify_form_factors)

    breakdown = staticmethod(substrate_breakdown)
    verify_b_qed_sm = staticmethod(verify_b_qed_sm)


substrate = _SubstrateNamespace()


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
    # substrate G_F (P7b §7.4 closure)
    "ALPHA_SUBSTRATE", "M_E_GEV", "V_EW_PDG_GEV",
    "v_ew_substrate", "fermi_constant_substrate", "fermi_constant_from_vev",
    "precision_chain", "verify_substrate_gf", "precision_chain_summary",
    # substrate Wolfenstein CKM (P6b closure)
    "DELTA_CP_CKM",
    "wolfenstein_lambda", "wolfenstein_lambda_sq",
    "wolfenstein_A", "wolfenstein_A_sq",
    "wolfenstein_apex_magnitude", "wolfenstein_apex_magnitude_sq",
    "wolfenstein_rho_bar", "wolfenstein_eta_bar",
    "ckm_V_us", "ckm_V_cb", "ckm_V_ub", "ckm_V_td",
    "jarlskog_ckm", "jarlskog_coefficient", "ckm_matrix",
    "ckm_precision_chain", "verify_substrate_ckm",
    "ckm_precision_chain_summary",
    # NOTE: meson decay constants (heavy/light/vector + B_c, P7b §7.5/§2-3/§7.6)
    # were consolidated into nwt_substrate.particles.decay_constants in v0.3.1
    # and are NO LONGER re-exported here — import them from that canonical home.
    # substrate Dalitz form factors (P7b §7.7 closure)
    "FORM_FACTORS", "cos_theta_C_substrate", "sin_theta_C_substrate",
    "f_plus_substrate", "f_plus_for", "f_plus_leading_order",
    "form_factor_precision_chain", "verify_form_factors",
    "form_factor_precision_chain_summary",
    # substrate identities
    "substrate", "substrate_breakdown", "verify_b_qed_sm",
]
