"""
nwt.gravity
===========

The gravitation view of the substrate algebra.

The thesis: gravity in NWT is the macroscopic projection of the
substrate's K_7 Wilson amplitude.  Newton's `G` is structurally derived
as a function of α and m_e -- no free parameters.

The Einstein-Hilbert action arises via Sakharov-induced gravity from
substrate-condensate fluctuations (Paper 18 G1-G6: linearized AND full
nonlinear Einstein equations from the NWT condensate).

Key substrate-derived results, exposed by this shim:

  Newton's G:
    LO  (Paper 14):    G = (8/7)² α²¹ ℏc / m_e²       (-0.24% from CODATA)
    NNLO (Paper 17):   G = 1 / M_Pl² with
                       m_e/M_Pl = (8/7) α^(21/2) (1 + α/7 + 3α²)
                       (~29 ppm from CODATA)

  Each factor is a substrate primitive:
    8/7      = dim(Spin(7) spinor S) / dim(vector V) = 8/7
    α^(21/2) = K_7 Wilson amplitude on 21 edges (= K_7 edge count = dim Adj)
    α/7      = NLO per-K_7-vertex correction (7 = K_7 vertex count)
    3α²      = NNLO PSL(2,7) Fano correction (3 = rank so(7) = 21/7)

The same K_7 graph hosts the matter sector ((p,q)-windings → Paper 6
particle masses) and the gravity sector (21-edge Wilson tower → Newton's
G).  See nwt.string.K7_torus for the toroidal embedding view.

Quick start::

    import nwt_substrate.gravity as grav

    # Substrate-derived constants
    print(grav.G_substrate_SI())            # ≈ 6.6745 × 10⁻¹¹ m³/kg/s²
                                            # (29 ppm from CODATA 6.6743e-11)
    print(grav.M_PLANCK_GEV)                # 1.221e19 GeV
    print(grav.m_e_over_M_Pl_NNLO())        # ≈ 4.186 × 10⁻²³

    # Structural breakdown (8/7 × α^(21/2) × ...)
    print(grav.G_substrate_breakdown())

    # GR observables
    bh = grav.black_hole_summary("proton", use_substrate_G=True)
    print(f"r_s(proton) = {bh['r_s_m']:.3e} m")
    print(f"T_H(proton) = {bh['T_H_K']:.3e} K")

The substrate's gravity is *not* a separate force.  It is the macroscopic
expectation value of the K_7 Wilson amplitude on the substrate condensate,
which is also responsible (via different substrate readings) for the SM
particle masses and the Standard-Model gauge content.
"""

from .constants import (
    HBAR_J_S, C_LIGHT_M_S,
    G_NEWTON_SI,
    M_ELECTRON_KG, M_PROTON_KG, M_NEUTRON_KG,
    M_PLANCK_GEV, M_PLANCK_REDUCED_GEV,
    H_0_KM_S_MPC, H_0_GEV,
    LAMBDA_CC_OBSERVED_GEV4, RHO_CRIT_GEV4,
    ALPHA_QED, M_E_GEV, M_P_GEV,
)

from .coupling import (
    m_e_over_M_Pl_LO,
    m_e_over_M_Pl_NLO,
    m_e_over_M_Pl_NNLO,
    m_e_over_M_Pl_observed,
    G_substrate_LO_natural,
    G_substrate_NNLO_natural,
    G_substrate_LO_SI,
    G_substrate_SI,
    G_substrate_breakdown,
)

from .black_holes import (
    schwarzschild_radius_m,
    hawking_temperature_K,
    hawking_luminosity_W,
    evaporation_time_s,
    PARTICLE_MASSES_KG,
    black_hole_summary,
)

from .einstein import (
    M_Pl_over_m_e_squared_observed,
    M_Pl_over_m_e_squared_substrate,
    UV_IR_bridge_breakdown,
    einstein_hilbert_coefficient,
    sakharov_route_summary,
    verify_schwarzschild_vacuum_symbolic,
    linearized_graviton_propagator_text,
)

from .flrw_test import (
    Q0_from_density_variance,
    scaling_index_n_NWT,
    R_avg_scaling,
    K_BR_scaling,
    A_of_z_scaling,
    O_of_z_scaling,
    C_prime_scaling,
    FLRWTestPrediction,
    nwt_prediction,
    OBSERVED_O_BAND_KMSMPC2,
    in_observed_O_band,
    comparison_summary,
    amplitude_check,
    SIGMA_M_DEFAULT,
    SIGMA_8_PLANCK,
)


__all__ = [
    # constants
    "HBAR_J_S", "C_LIGHT_M_S", "G_NEWTON_SI",
    "M_ELECTRON_KG", "M_PROTON_KG", "M_NEUTRON_KG",
    "M_PLANCK_GEV", "M_PLANCK_REDUCED_GEV",
    "H_0_KM_S_MPC", "H_0_GEV",
    "LAMBDA_CC_OBSERVED_GEV4", "RHO_CRIT_GEV4",
    "ALPHA_QED", "M_E_GEV", "M_P_GEV",
    # substrate-derived G
    "m_e_over_M_Pl_LO", "m_e_over_M_Pl_NLO", "m_e_over_M_Pl_NNLO",
    "m_e_over_M_Pl_observed",
    "G_substrate_LO_natural", "G_substrate_NNLO_natural",
    "G_substrate_LO_SI", "G_substrate_SI",
    "G_substrate_breakdown",
    # black holes
    "schwarzschild_radius_m", "hawking_temperature_K",
    "hawking_luminosity_W", "evaporation_time_s",
    "PARTICLE_MASSES_KG", "black_hole_summary",
    # einstein equations
    "M_Pl_over_m_e_squared_observed", "M_Pl_over_m_e_squared_substrate",
    "UV_IR_bridge_breakdown", "einstein_hilbert_coefficient",
    "sakharov_route_summary", "verify_schwarzschild_vacuum_symbolic",
    "linearized_graviton_propagator_text",
    # FLRW null tests (Buchert backreaction predictions)
    "Q0_from_density_variance", "scaling_index_n_NWT",
    "R_avg_scaling", "K_BR_scaling", "A_of_z_scaling", "O_of_z_scaling",
    "C_prime_scaling",
    "FLRWTestPrediction", "nwt_prediction",
    "OBSERVED_O_BAND_KMSMPC2", "in_observed_O_band",
    "comparison_summary", "amplitude_check",
    "SIGMA_M_DEFAULT", "SIGMA_8_PLANCK",
]
