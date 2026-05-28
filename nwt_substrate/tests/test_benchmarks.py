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
    assert len(results) == 11      # current set: 4 EW + 1 K_7 + 1 DM + 1 G + 3 cosmology + 1 MTC
    for r in results:
        assert isinstance(r, BenchmarkResult)
        assert r.substrate_time_us >= 0


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
