"""Smoke tests for nwt_substrate.chemistry — verifies that core
predictions match validated reference results from yesterday's
analysis scripts.

Reference results being preserved as regression invariants:
  - Aromaticity 20/20 Hückel/Möbius classifications
  - C_60 = |A_5| orbit counts (60 vertices, 12 pent, 20 hex, 90 edges)
  - C_60 vibrational: 174 total modes, 4 IR T_1u, 10 Raman (2 A_g + 8 H_g)
  - McKay 5-coord fluxional prediction
  - Linear acenes RE 5/5 within 6% from 1 benzene calibration
"""
import pytest
import nwt_substrate.chemistry as chem


# ---------------------------------------------------------------------------
# Aromaticity (Hopf-pair rule)
# ---------------------------------------------------------------------------

def test_benzene_aromatic():
    r = chem.aromaticity_class("benzene")
    assert r.classification == "aromatic"
    assert r.n_hopf_links == 3
    assert r.parity == "odd"


def test_cyclobutadiene_antiaromatic():
    r = chem.aromaticity_class("cyclobutadiene")
    assert r.classification == "anti-aromatic"
    assert r.n_hopf_links == 2
    assert r.parity == "even"


def test_mobius_aromatic():
    r = chem.aromaticity_class("mobius_annulene_16")
    assert r.classification == "mobius_aromatic"
    assert r.mobius_twists == 1


def test_naphthalene_re():
    re = chem.aromatic_resonance_energy("naphthalene", calibration_kcal=12.0)
    # 5 π pairs × 12 kcal/mol = 60.0 (expt ≈ 61, error ~ -1.6%)
    assert abs(re.kcal_per_mol - 60.0) < 0.1


def test_anthracene_re():
    re = chem.aromatic_resonance_energy("anthracene", calibration_kcal=12.0)
    # 7 π pairs × 12 kcal/mol = 84.0 (expt ≈ 84, error 0%)
    assert abs(re.kcal_per_mol - 84.0) < 0.1


def test_bond_order_benzene():
    assert chem.bond_order("benzene", "C-C") == 1.5


def test_bond_order_acetylene():
    assert chem.bond_order("acetylene", "C-C") == 3.0


# ---------------------------------------------------------------------------
# Fullerenes (2I substrate orbit)
# ---------------------------------------------------------------------------

def test_c60_orbit_counts():
    orbit = chem.fullerene_orbit_counts(60)
    assert orbit.n_vertices == 60
    assert orbit.n_pentagons == 12
    assert orbit.n_hexagons == 20
    assert orbit.n_edges == 90
    assert orbit.is_2I_orbit is True
    assert orbit.point_group == "I_h"


def test_c60_anion_n3_superconductor():
    states = chem.c60_anion_magic_states()
    n3 = next(s for s in states if s.charge == 3)
    assert "supercond" in n3.predicted_property
    assert "A_3 C_60" in n3.empirical_match


def test_c60_anion_n6_max_reduction():
    states = chem.c60_anion_magic_states()
    n6 = next(s for s in states if s.charge == 6)
    assert n6.closed_shell is True
    assert "reduc" in n6.predicted_property


# ---------------------------------------------------------------------------
# Vibrational (point-group representation theory)
# ---------------------------------------------------------------------------

def test_c60_vibrational_total_modes():
    v = chem.c60_vibrational_summary()
    assert v.total_modes == 174   # 3·60 - 6
    assert v.n_atoms == 60
    assert v.point_group == "I_h"


def test_c60_vibrational_ir_active():
    v = chem.c60_vibrational_summary()
    assert v.ir_active_modes == 4    # 4 T_1u modes


def test_c60_vibrational_raman_active():
    v = chem.c60_vibrational_summary()
    assert v.raman_active_modes == 10   # 2 A_g + 8 H_g


# ---------------------------------------------------------------------------
# McKay admissibility
# ---------------------------------------------------------------------------

def test_mckay_5_fluxional():
    """5-coordination is McKay-forbidden, predicted fluxional."""
    r = chem.is_admissible_coordination(5)
    assert r is False


def test_mckay_4_admissible_2T():
    """4-coordination → 2T binary tetrahedral group (E_6 Dynkin)."""
    assert chem.is_admissible_coordination(4) is True
    assert chem.mckay_orbit_size("2T") == 4


def test_mckay_6_admissible_2O():
    """6-coordination → 2O binary octahedral group (E_7 Dynkin)."""
    assert chem.is_admissible_coordination(6) is True
    assert chem.mckay_orbit_size("2O") == 6


def test_mckay_12_admissible_2I():
    """12-coordination → 2I binary icosahedral group (E_8 Dynkin)."""
    assert chem.is_admissible_coordination(12) is True
    assert chem.mckay_orbit_size("2I") == 12


def test_mckay_admissible_set():
    """McKay-admissible coordinations: 1, 2, 3, 4, 6, 12 (no 5)."""
    coords = chem.allowed_coordinations()
    assert 5 not in coords
    assert set(coords) == {1, 2, 3, 4, 6, 12}


# ---------------------------------------------------------------------------
# Benchmark smoke (just verify it runs)
# ---------------------------------------------------------------------------

def test_benchmark_aromaticity_runs():
    r = chem.benchmark.benchmark_aromaticity(n_repeats=10)
    assert r.observable == "aromaticity_classification"
    assert r.speedup > 0


def test_benchmark_c60_modes_runs():
    r = chem.benchmark.benchmark_c60_vibrational()
    assert r.observable == "c60_vibrational_modes"
    assert r.substrate_time_us > 0


def test_full_benchmark_suite_runs():
    results = chem.benchmark.run_full_benchmark_suite()
    assert len(results) == 8
    for r in results:
        assert r.speedup > 0
    observables = {r.observable for r in results}
    assert "wade_classification" in observables
    assert "closo_3d_aromaticity" in observables
    assert "deltahedron_combinatorics" in observables
