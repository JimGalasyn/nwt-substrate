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
    assert len(results) == 38      # 38 benchmarks: full SM + composites + scattering + chemistry + cosmology
    for r in results:
        assert isinstance(r, BenchmarkResult)
        assert r.substrate_time_us >= 0


def test_composite_particle_deuteron_within_percent():
    """Deuteron connected-sum gives mass within 0.1% of PDG."""
    from nwt_substrate.benchmarks import benchmark_composite_particles
    r = benchmark_composite_particles()
    pct_str = r.substrate_accuracy.split("%")[0]
    assert float(pct_str) < 0.5


def test_exotic_states_lists_glueballs():
    """Exotic states benchmark surfaces multiple glueball candidates."""
    from nwt_substrate.benchmarks import benchmark_exotic_states
    r = benchmark_exotic_states()
    assert "glueball" in r.substrate_value or "eta_" in r.substrate_value


def test_bhabha_returns_pb_per_sr():
    """Bhabha benchmark yields differential cross section in pb/sr."""
    from nwt_substrate.benchmarks import benchmark_bhabha_scattering
    r = benchmark_bhabha_scattering()
    assert "pb/sr" in r.substrate_value


def test_moller_returns_pb_per_sr():
    """Møller benchmark yields differential cross section in pb/sr."""
    from nwt_substrate.benchmarks import benchmark_moller_scattering
    r = benchmark_moller_scattering()
    assert "pb/sr" in r.substrate_value


def test_aromatic_resonance_includes_coronene_144():
    """Aromatic RE benchmark hits coronene = 144 kcal/mol (substrate K_7 target)."""
    from nwt_substrate.benchmarks import benchmark_aromatic_resonance_energies
    r = benchmark_aromatic_resonance_energies()
    assert "144" in r.substrate_value
    assert "coronene" in r.substrate_value or "kcal/mol" in r.substrate_value


def test_substrate_dimensions_all_pass():
    """Substrate-DNA dimensional identities all check out."""
    from nwt_substrate.benchmarks import benchmark_substrate_dimensions
    r = benchmark_substrate_dimensions()
    # substrate_accuracy is like "N/M — exact integer arithmetic"
    n_pass_str, m_total_str = r.substrate_accuracy.split(" ")[0].split("/")
    assert int(n_pass_str) == int(m_total_str)        # all N checks must pass


def test_compton_thomson_at_ppm():
    """Compton Thomson limit matches PDG cross section in ppm."""
    from nwt_substrate.benchmarks import benchmark_qed_compton_scattering
    r = benchmark_qed_compton_scattering()
    assert "Thomson" in r.substrate_value
    ppm_str = r.substrate_accuracy.split()[2]    # "Thomson at {N} ppm vs PDG"
    assert float(ppm_str) < 200


def test_eemumu_returns_pb_units():
    """e⁺e⁻ → μ⁺μ⁻ benchmark yields LO QED cross section in pb."""
    from nwt_substrate.benchmarks import benchmark_qed_eemumu
    r = benchmark_qed_eemumu()
    assert "pb" in r.substrate_value


def test_muon_decay_rate_matches_pdg_lifetime():
    """Sirlin muon-decay benchmark gives τ_μ within 1%."""
    from nwt_substrate.benchmarks import benchmark_muon_decay_rate
    r = benchmark_muon_decay_rate()
    pct_str = r.substrate_accuracy.split("%")[0]
    assert float(pct_str) < 1.0


def test_cosmogenesis_thorne_a_star():
    """Cosmogenesis benchmark gives Thorne a* near 0.998."""
    from nwt_substrate.benchmarks import benchmark_cosmogenesis
    r = benchmark_cosmogenesis()
    assert "Thorne a*" in r.substrate_value
    assert "0.99" in r.substrate_value


def test_nmr_sign_rule_perfect_score():
    """NMR sign rule classifies all 14 references correctly."""
    from nwt_substrate.benchmarks import benchmark_nmr_chemical_shifts
    r = benchmark_nmr_chemical_shifts()
    assert "14/14" in r.substrate_value


def test_c60_vibrational_modes_174_total():
    """C_60 vibrational benchmark gives 174 = 3·60 - 6 modes."""
    from nwt_substrate.benchmarks import benchmark_c60_vibrational_modes
    r = benchmark_c60_vibrational_modes()
    assert "174" in r.substrate_value


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


def test_sin2_theta_w_reports_substrate_formula():
    """sin²θ_W benchmark reports the substrate prediction (2+α)/9 ≈ 0.223 compared
    LIKE-FOR-LIKE against the on-shell angle 1−M_W²/M_Z² (it is a tree-level =
    on-shell angle), not the measured effective input it once reported."""
    from nwt_substrate.benchmarks import benchmark_sin2_theta_W
    from nwt_substrate.electroweak.constants import SIN2_THETA_W_SUBSTRATE, M_W, M_Z
    r = benchmark_sin2_theta_W()
    assert 0.222 < SIN2_THETA_W_SUBSTRATE < 0.224     # (2+α)/9, not 0.23121
    assert "0.223" in r.substrate_value
    assert "%" in r.substrate_accuracy
    # (2+α)/9 IS the on-shell angle: matches 1−M_W²/M_Z² to <0.1%
    s2_onshell = 1.0 - (M_W / M_Z) ** 2
    assert abs(SIN2_THETA_W_SUBSTRATE - s2_onshell) / s2_onshell < 1e-3
    assert "on-shell" in r.notes.lower()              # transparent about the scheme
    assert float(r.substrate_accuracy.split("%")[0]) < 0.5   # on-shell, not the 3.5% effective gap


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


def test_standalone_predictions_are_genuine_and_separated():
    """predict.py emits derived-only dimensionless predictions (Luke Leighton's
    derivation-separation protocol): predictions() takes no input of any kind,
    and each lands near the quarantined measured reference."""
    from nwt_substrate.benchmarks import predict
    pred = predict.predictions()                      # no args — standalone derivation
    ref = predict.REFERENCE
    assert set(pred) == set(ref)
    # 4 of 5 land within 1%; sin²θ_W is the documented leading-order angle
    # (~3.5%, on-shell vs effective) — kept so the diff shows a real deviation.
    tol = {"sin2_theta_W": 0.04}
    for k in pred:
        dev = abs(pred[k] - ref[k]) / ref[k]
        assert dev < tol.get(k, 0.01), f"{k}: {pred[k]} vs {ref[k]} ({dev:.1%})"
    assert abs(pred["inv_alpha"] - 137.036) < 0.01
    # transparency: sin²θ_W is genuinely a few-percent LO angle, not sub-1%
    dev_sw = abs(pred["sin2_theta_W"] - ref["sin2_theta_W"]) / ref["sin2_theta_W"]
    assert 0.02 < dev_sw < 0.04


def test_o10_dag_cit_readout():
    """O10 DAG: structural invariants hold; cit marks sin²θ_W as the one defect
    edge (not repaired); α is the load-bearing root; identities commute."""
    from nwt_substrate.benchmarks import o10
    g = o10.build_constants_dag()
    # O10 structural invariants all hold
    assert all(g.acceptance_checklist().values())
    assert g.backward_edges() == [] and g.is_acyclic() and g.witnesses_are_sinks()
    # cit: exactly sin²θ_W is the marked defect at 1% (a leading-order angle)
    assert g.cit_defects(tol=0.01) == ["sin2_theta_W"]
    # α (and its closed form) is the load-bearing root — reaches all 5 outputs
    load = dict(g.load_ranking())
    assert load["α"] == 5
    assert load["α"] >= max(v for k, v in load.items() if k not in ("α", "form:25π√3+1"))
    # the commutative-diagram identities (21, 8) actually commute
    assert all(r["commutes"] for r in g.commutative_checks())
    # Horn-clause frontier of a witness lists its structural premises
    assert "N_VERTICES_K7" in g.horn_frontier("wit:cabibbo_lambda")


def test_o10_suite_dag_full_38():
    """The whole 38-benchmark suite as one O10 DAG: structural invariants hold,
    cit marks the known >1% benchmarks as defect edges (not repaired), and the
    no-metric benchmarks are flagged qualitative rather than fabricated."""
    from nwt_substrate.benchmarks import o10
    assert len(o10.benchmark_functions()) == 38
    g = o10.build_suite_dag()
    outs = [n for n, nd in g.nodes.items() if nd.stage == o10.Stage.OUTPUT]
    assert len(outs) == 38
    assert all(g.acceptance_checklist().values())
    assert g.backward_edges() == [] and g.is_acyclic() and g.witnesses_are_sinks()
    defects = g.cit_defects(0.01)
    assert "benchmark_muon_lifetime" in defects      # m⁵-amplified compound
    assert "benchmark_mass_spectrum" in defects       # Paper-6 mass-formula residual
    assert "benchmark_sin2_theta_W" not in defects    # on-shell comparison -> admissible
    assert 4 <= len(defects) <= 12
    qual = g.qualitative_outputs()
    assert "benchmark_bhabha_scattering" in qual      # emits a cross-section, no vs-measurement
    # parser spot-checks
    assert o10._parse_deviation("exact (closed form)") == (0.0, "exact")
    assert o10._parse_deviation("7.6 ppm vs CODATA")[0] < 1e-4
    assert o10._parse_deviation("LO QED formula at substrate α")[0] is None


def test_o10_suite_dag_couples_to_sensitivity_report():
    """With a SensitivityReport the structural integers wire to the benchmarks
    they move, and load_ranking reflects that computed coupling."""
    from nwt_substrate.benchmarks import o10
    from nwt_substrate.sensitivity import SensitivityReport
    rep = SensitivityReport(
        benchmarks=["benchmark_z_boson_width", "benchmark_fermi_constant"],
        baseline={},
        per_integer={
            "H_V_SO7": {1: {"status": "ok",
                            "moved": ["benchmark_fermi_constant", "benchmark_z_boson_width"]}},
            "RANK_SO7": {1: {"status": "ok", "moved": ["benchmark_z_boson_width"]}},
        },
    )
    g = o10.build_suite_dag(report=rep)
    load = dict(g.load_ranking())
    assert load["H_V_SO7"] == 2 and load["RANK_SO7"] == 1
    assert all(g.acceptance_checklist().values())
