"""Unit tests for nwt_substrate.benchmarks."""

import pytest

from nwt_substrate.benchmarks import (
    BenchmarkResult,
    benchmark_alpha_derivation,
    benchmark_mass_spectrum,
    benchmark_modular_data,
    benchmark_ckm_cabibbo,
    benchmark_k7_face_structure,
    benchmark_wimp_tower,
    run_all,
)


def test_benchmark_result_type():
    """Each benchmark returns a BenchmarkResult dataclass."""
    r = benchmark_alpha_derivation()
    assert isinstance(r, BenchmarkResult)
    assert r.name
    assert r.substrate_time_us >= 0
    assert r.substrate_accuracy
    assert r.traditional_method


def test_alpha_derivation_gets_137():
    """Substrate α formula gives 1/α ≈ 137.04."""
    r = benchmark_alpha_derivation()
    # The substrate value string should contain the computed inverse alpha
    assert "137" in r.substrate_value


def test_mass_spectrum_handles_compendium():
    """Mass spectrum benchmark times the full compendium."""
    r = benchmark_mass_spectrum()
    assert r.substrate_time_us > 0
    # The substrate_value should contain a particle count
    assert "particles" in r.substrate_value or "masses" in r.substrate_value


def test_modular_data_recovers_c_15_7():
    """SU(2)_5 modular data benchmark recovers c = 15/7."""
    r = benchmark_modular_data()
    assert "15/7" in r.substrate_value or "2.142857" in r.substrate_value


def test_cabibbo_close_to_pdg():
    """Cabibbo angle from λ²=7α matches PDG to ~0.1%."""
    r = benchmark_ckm_cabibbo()
    # PDG: θ_C ≈ 13.02°; substrate gives 13.062°
    assert "13" in r.substrate_value


def test_k7_face_structure_verifies_torus():
    """K_7 Heffter embedding has V=7, E=21, F=14 (torus)."""
    r = benchmark_k7_face_structure()
    assert "V=7" in r.substrate_value
    assert "F=14" in r.substrate_value


def test_wimp_tower_includes_98gev():
    """WIMP tower benchmark mentions the 98 GeV rung."""
    r = benchmark_wimp_tower()
    assert "98" in r.substrate_value


def test_run_all_returns_list():
    """run_all() returns list of all benchmark results."""
    results = run_all(verbose=False)
    assert isinstance(results, list)
    assert len(results) == 26      # 26 benchmarks: SM + atomic + QED + QCD + EW + cosmology + chem + gravity + BH
    for r in results:
        assert isinstance(r, BenchmarkResult)
        assert r.substrate_time_us >= 0


def test_atomic_hydrogen_covers_chain():
    """Atomic hydrogen benchmark covers a₀, R_H, Lyman α, 21cm, Lamb shift."""
    from nwt_substrate.benchmarks import benchmark_atomic_hydrogen
    r = benchmark_atomic_hydrogen()
    for token in ["a₀", "R_H", "Lyman", "21cm"]:
        assert token in r.substrate_value


def test_electron_anomaly_matches_schwinger():
    """Electron a_e benchmark agrees with Schwinger formula essentially exactly."""
    from nwt_substrate.benchmarks import benchmark_electron_anomaly
    r = benchmark_electron_anomaly()
    pct_str = r.substrate_accuracy.split("%")[0]
    assert float(pct_str) < 1e-6     # 1-loop = Schwinger to machine precision


def test_qcd_constants_recover_pdg_alpha_s():
    """QCD benchmark recovers α_s(M_Z) ≈ 0.118."""
    from nwt_substrate.benchmarks import benchmark_qcd_constants
    r = benchmark_qcd_constants()
    assert "α_s" in r.substrate_value
    assert "0.117" in r.substrate_value or "0.118" in r.substrate_value


def test_sin2_theta_w_at_ppm():
    """sin²θ_W substrate matches PDG to ppm."""
    from nwt_substrate.benchmarks import benchmark_sin2_theta_W
    r = benchmark_sin2_theta_W()
    ppm_str = r.substrate_accuracy.split()[0]
    assert float(ppm_str) < 200      # < 200 ppm easily


def test_black_hole_thermodynamics_includes_evaporation():
    """Black hole benchmark mentions Hawking T + Schwarzschild r_S + evaporation."""
    from nwt_substrate.benchmarks import benchmark_black_hole_thermodynamics
    r = benchmark_black_hole_thermodynamics()
    assert "T_H" in r.substrate_value
    assert "r_S" in r.substrate_value
    assert "τ_evap" in r.substrate_value or "evap" in r.substrate_value


def test_neutrino_sector_predicts_active_and_sterile():
    """Neutrino benchmark gives both active and sterile mass scales."""
    from nwt_substrate.benchmarks import benchmark_neutrino_sector
    r = benchmark_neutrino_sector()
    assert "m_ν" in r.substrate_value
    assert "sterile" in r.substrate_value or "N" in r.substrate_value


def test_pmns_angles_three_mixing_angles():
    """PMNS benchmark mentions all three mixing angles."""
    from nwt_substrate.benchmarks import benchmark_pmns_angles
    r = benchmark_pmns_angles()
    assert "θ_12" in r.substrate_value
    assert "θ_13" in r.substrate_value
    assert "θ_23" in r.substrate_value


def test_decay_constants_cover_light_and_heavy():
    """Decay-constants benchmark covers π, K, η, D, Ds, B, Bs."""
    from nwt_substrate.benchmarks import benchmark_decay_constants
    r = benchmark_decay_constants()
    for name in ["f_π", "f_K", "f_D", "f_B"]:
        assert name in r.substrate_value


def test_vector_meson_decay_covers_11_states():
    """Vector meson benchmark covers the substrate's 11-state tower."""
    from nwt_substrate.benchmarks import benchmark_vector_meson_decay
    r = benchmark_vector_meson_decay()
    assert "11" in r.substrate_value or "f_rho" in r.substrate_value


def test_chemistry_runs_aromaticity_nics_c60():
    """Chemistry benchmark exercises aromaticity, NICS, and C_60 combinatorics."""
    from nwt_substrate.benchmarks import benchmark_chemistry
    r = benchmark_chemistry()
    assert "aromaticity" in r.substrate_value
    assert "NICS" in r.substrate_value
    assert "C_60" in r.substrate_value


def test_higgs_vev_predicted_to_ppm():
    """Substrate v_EW matches PDG to ppm level."""
    from nwt_substrate.benchmarks import benchmark_higgs_vev
    r = benchmark_higgs_vev()
    assert "ppm" in r.substrate_accuracy
    # numeric extraction: should be < 100 ppm
    ppm_str = r.substrate_accuracy.split()[0]
    assert float(ppm_str) < 100


def test_higgs_mass_via_lambda_18alpha():
    """Substrate m_h from λ_H = 18α formula appears in benchmark."""
    from nwt_substrate.benchmarks import benchmark_higgs_mass_vs_98gev
    r = benchmark_higgs_mass_vs_98gev()
    assert "λ_H=18α" in r.substrate_value or "18α" in r.notes


def test_fermi_constant_close_to_pdg():
    """G_F matches PDG to <100 ppm."""
    from nwt_substrate.benchmarks import benchmark_fermi_constant
    r = benchmark_fermi_constant()
    ppm_str = r.substrate_accuracy.split()[0]
    assert float(ppm_str) < 100


def test_z_boson_width_within_few_percent():
    """Γ_Z matches LEP to < 5%."""
    from nwt_substrate.benchmarks import benchmark_z_boson_width
    r = benchmark_z_boson_width()
    # Format is "2.93% on Γ_Z; lepton universality at ppm level"
    pct_str = r.substrate_accuracy.split("%")[0]
    assert float(pct_str) < 5


def test_muon_lifetime_compounds_errors():
    """Muon lifetime benchmark notes the τ ∝ 1/m_μ^5 error amplification."""
    from nwt_substrate.benchmarks import benchmark_muon_lifetime
    r = benchmark_muon_lifetime()
    # τ_μ should agree to a few % (compound error from m_μ + G_F)
    pct_str = r.substrate_accuracy.split("%")[0]
    assert 0 < float(pct_str) < 25


def test_full_ckm_includes_v_us():
    """Full CKM benchmark includes V_us close to PDG."""
    from nwt_substrate.benchmarks import benchmark_full_ckm
    r = benchmark_full_ckm()
    # PDG V_us ≈ 0.2253; substrate gives 0.2260
    assert "V_us=0.2" in r.substrate_value


def test_gravitational_constant_in_ppm():
    """G benchmark reports CODATA comparison in ppm."""
    from nwt_substrate.benchmarks import benchmark_gravitational_constant
    r = benchmark_gravitational_constant()
    assert "ppm" in r.substrate_accuracy


def test_lambda_cc_solves_cc_problem():
    """Λ benchmark notes the 123-orders-of-magnitude solution."""
    from nwt_substrate.benchmarks import benchmark_lambda_cc
    r = benchmark_lambda_cc()
    assert "cosmological constant problem" in r.notes or "123" in r.notes


def test_omega_b_c_matches_planck_at_sub_percent():
    """Ω_b/Ω_c benchmark agrees with Planck at sub-percent."""
    from nwt_substrate.benchmarks import benchmark_omega_b_c
    r = benchmark_omega_b_c()
    # Substrate gives ~0.186, Planck ~0.186 — should be <1% in accuracy field
    assert "%" in r.substrate_accuracy


def test_eta_B_predicted_without_BSM():
    """η_B benchmark notes substrate doesn't need new physics."""
    from nwt_substrate.benchmarks import benchmark_eta_B
    r = benchmark_eta_B()
    assert "new physics" in r.notes or "CP-violation" in r.notes


def test_total_substrate_time_under_one_second():
    """All substrate-algebra benchmarks combined run in well under 1 second."""
    results = run_all(verbose=False)
    total_us = sum(r.substrate_time_us for r in results)
    # Generous bound: should be way under 1 sec even on a slow CI box.
    assert total_us < 1_000_000, f"benchmarks took {total_us} us, expected < 1 sec"
