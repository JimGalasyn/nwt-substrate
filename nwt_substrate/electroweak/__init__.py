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
    # substrate identities
    "substrate", "substrate_breakdown", "verify_b_qed_sm",
]
