"""
test_neutrino_shim.py
=====================

Tests for the nwt.neutrino shim — Paper 20 predictions from the
substrate algebra (K_8 extension + Spin(8) triality + π_1(PSU(3))
winding).
"""

from __future__ import annotations

import math
import pytest

import nwt_substrate.neutrino as nu
from nwt_substrate.isa.constants import (
    ALPHA_NWT,
    K8_ACTIVE_EDGE_COUNT,
    K8_PARTITION,
    K8_SEESAW_EDGE_DIFFERENCE,
    K8_STERILE_EDGE_COUNT,
    N_EDGES_K8,
    N_VERTICES_K8,
    RANK_SO7,
)


# ---------------------------------------------------------------------------
# K_8 structural identities (sanity — also enforced at import in constants.py)
# ---------------------------------------------------------------------------

def test_k8_partition_sums_to_edges():
    """K_8 partition (6, 3, 12, 1, 6) sums to N_EDGES_K8 = 28."""
    assert sum(K8_PARTITION) == N_EDGES_K8 == 28


def test_k8_active_minus_sterile_is_nine():
    """Seesaw edge difference is 9 = 12 (SM-flavor) − 3 (internal SU(2))."""
    assert K8_ACTIVE_EDGE_COUNT - K8_STERILE_EDGE_COUNT == 9
    assert K8_SEESAW_EDGE_DIFFERENCE == 9


def test_k8_vertices_eight():
    """K_8 has 8 vertices."""
    assert N_VERTICES_K8 == 8


# ---------------------------------------------------------------------------
# Active neutrino masses (Paper 19 W3.3-J II / Paper 20 §3)
# ---------------------------------------------------------------------------

def test_m1_matches_paper20():
    """m_1 from K_8 Wilson amplitude ≈ 14.84 meV (Paper 20 abstract)."""
    m1_eV = nu.active_mass_lightest()
    m1_meV = m1_eV * 1.0e3
    # Paper 20 states 14.84 meV; allow 0.1 meV tolerance (substrate
    # constants vs NLO/NNLO numerical precision)
    assert m1_meV == pytest.approx(14.84, abs=0.1)


def test_active_masses_normal_hierarchy_ordering():
    """m_1 < m_2 < m_3 (normal hierarchy)."""
    m1, m2, m3 = nu.active_masses()
    assert m1 < m2 < m3


def test_active_masses_squared_differences_match_pdg():
    """m_2² − m_1² = Δm²_21 and m_3² − m_1² = Δm²_31 by construction."""
    m1, m2, m3 = nu.active_masses()
    assert (m2 ** 2 - m1 ** 2) == pytest.approx(nu.DELTA_M_SQ_21_eV2, rel=1e-9)
    assert (m3 ** 2 - m1 ** 2) == pytest.approx(nu.DELTA_M_SQ_31_eV2, rel=1e-9)


def test_sum_active_masses_below_cosmological_bound():
    """Σ m_ν ≈ 85 meV, comfortably below the 120 meV Planck CMB+BAO bound."""
    sigma = nu.sum_active_masses()
    sigma_meV = sigma * 1.0e3
    assert 80.0 < sigma_meV < 90.0


def test_active_mass_order_consistency():
    """LO < NLO < NNLO (NLO and NNLO corrections are both positive)."""
    m_LO = nu.active_mass_lightest(order="LO")
    m_NLO = nu.active_mass_lightest(order="NLO")
    m_NNLO = nu.active_mass_lightest(order="NNLO")
    assert m_LO < m_NLO < m_NNLO


def test_active_mass_invalid_order_raises():
    """Unknown order raises ValueError."""
    with pytest.raises(ValueError):
        nu.active_mass_lightest(order="NNNLO")


# ---------------------------------------------------------------------------
# Sterile neutrino masses (Paper 20 §6)
# ---------------------------------------------------------------------------

def test_sterile_masses_paper20_window():
    """Sterile masses fall in the νMSM 50-250 MeV window of Paper 20."""
    m_N1, m_N2, m_N3 = nu.sterile_masses()
    for m in (m_N1, m_N2, m_N3):
        assert 5.0e7 < m < 2.5e8, (
            f"Sterile mass {m:.3e} eV outside νMSM window [50, 250] MeV"
        )


def test_sterile_mass_n1_matches_paper20():
    """m_N1 ≈ 61.3 MeV (Paper 20 abstract)."""
    m_N1, _, _ = nu.sterile_masses()
    m_N1_MeV = m_N1 / 1.0e6
    assert m_N1_MeV == pytest.approx(61.3, abs=0.5)


def test_sterile_mass_n2_matches_paper20():
    """m_N2 ≈ 70.8 MeV (Paper 20 abstract)."""
    _, m_N2, _ = nu.sterile_masses()
    m_N2_MeV = m_N2 / 1.0e6
    assert m_N2_MeV == pytest.approx(70.8, abs=0.5)


def test_sterile_mass_n3_in_paper20_range():
    """m_N3 ≈ 218.8 MeV (Paper 20 abstract); allow 5 MeV for Δm²_31 anchor choice."""
    _, _, m_N3 = nu.sterile_masses()
    m_N3_MeV = m_N3 / 1.0e6
    assert m_N3_MeV == pytest.approx(218.8, abs=5.0)


def test_sterile_active_mixing_sq_matches_paper20():
    """|U_α4|² ≈ 2.4 × 10⁻¹⁰ (Paper 20 abstract)."""
    u_sq = nu.sterile_active_mixing_sq()
    assert u_sq == pytest.approx(2.4e-10, rel=0.05)


def test_sterile_active_mixing_is_alpha_to_nine_halves():
    """|U_α4|² = α^(9/2) by construction."""
    expected = ALPHA_NWT ** 4.5
    assert nu.sterile_active_mixing_sq() == pytest.approx(expected, rel=1e-9)


def test_seesaw_ratio_substrate_identity():
    """seesaw_ratio = α^(K8_SEESAW_EDGE_DIFFERENCE / 2)."""
    expected = ALPHA_NWT ** (K8_SEESAW_EDGE_DIFFERENCE / 2.0)
    assert nu.seesaw_ratio() == pytest.approx(expected, rel=1e-12)


def test_seesaw_ratio_consistency():
    """m_sterile / m_active = 1 / seesaw_ratio = α^(-9/2) ≈ 4.13×10⁹."""
    active = nu.active_masses()
    sterile = nu.sterile_masses()
    ratios = [s / a for s, a in zip(sterile, active)]
    expected = 1.0 / nu.seesaw_ratio()
    for r in ratios:
        assert r == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# PMNS mixing angles
# ---------------------------------------------------------------------------

def test_theta_12_arctan_one_over_sqrt_two():
    """θ_12 = arctan(1/√2) ≈ 35.26°."""
    a = nu.pmns_angles_leading_order()
    assert math.degrees(a.theta_12) == pytest.approx(35.26, abs=0.01)


def test_theta_23_pi_over_four():
    """θ_23 = π/4 = 45°."""
    a = nu.pmns_angles_leading_order()
    assert math.degrees(a.theta_23) == pytest.approx(45.0, abs=1e-9)


def test_theta_13_arcsin_sqrt_three_alpha():
    """θ_13 = arcsin(√(rank · α)) with rank = RANK_SO7 = 3 (Paper 19 W3.3-K)."""
    a = nu.pmns_angles_leading_order()
    expected = math.asin(math.sqrt(RANK_SO7 * ALPHA_NWT))
    assert a.theta_13 == pytest.approx(expected, rel=1e-9)


def test_theta_13_near_paper20_target():
    """θ_13 ≈ 8.5° from the substrate (Paper 20 §7); experimental ≈ 8.6°."""
    a = nu.pmns_angles_leading_order()
    assert math.degrees(a.theta_13) == pytest.approx(8.5, abs=0.2)


# ---------------------------------------------------------------------------
# Dirac CP phase from π_1(PSU(3)) winding
# ---------------------------------------------------------------------------

def test_delta_cp_minus_two_pi_over_three():
    """δ_CP = -2π/3 (Paper 20 §7.5, Baez Fano convention)."""
    delta = nu.delta_cp_winding()
    assert delta == pytest.approx(-2.0 * math.pi / 3.0, rel=1e-9)
    assert math.degrees(delta) == pytest.approx(-120.0, abs=1e-9)


def test_delta_cp_magnitude_is_winding_invariant():
    """|δ_CP| = 2π/3 = 120° is the rephasing-invariant content (Paper 20 §7.5)."""
    assert math.degrees(abs(nu.delta_cp_winding())) == pytest.approx(120.0, abs=1e-9)


def test_jarlskog_matches_paper20():
    """|J_CP| ≈ 0.029 (Paper 20 abstract)."""
    j = nu.jarlskog_invariant()
    assert abs(j) == pytest.approx(0.029, abs=0.002)


def test_jarlskog_sign_tracks_sin_delta_cp():
    """sign(J_CP) tracks sign(sin(δ_CP))."""
    j = nu.jarlskog_invariant()
    delta = nu.delta_cp_winding()
    # both negative under Baez Fano convention
    assert j < 0
    assert math.sin(delta) < 0


# ---------------------------------------------------------------------------
# Unified Wilson recipe (Paper 17 + 19 + 20)
# ---------------------------------------------------------------------------

def test_wilson_recipe_recovers_electron_ratio():
    """Wilson recipe with (N_v, N_e) = (7, 21) reproduces Paper 17 m_e."""
    m_e_eV = nu.wilson_mass_eV(N_v=7, N_e=21)
    # Paper 17: m_e ≈ 0.511 MeV = 5.11×10⁵ eV
    assert m_e_eV == pytest.approx(5.11e5, rel=0.01)


def test_wilson_recipe_active_neutrino_matches_module():
    """Wilson recipe with (8, 28) matches active_mass_lightest()."""
    direct = nu.wilson_mass_eV(N_v=N_VERTICES_K8, N_e=K8_ACTIVE_EDGE_COUNT)
    via_helper = nu.active_mass_lightest()
    assert direct == pytest.approx(via_helper, rel=1e-12)


def test_wilson_recipe_sterile_neutrino_matches_via_seesaw():
    """Wilson recipe (8, 19) gives the same answer as active_mass / seesaw_ratio."""
    direct = nu.wilson_mass_eV(N_v=N_VERTICES_K8, N_e=K8_STERILE_EDGE_COUNT)
    m1_active = nu.active_mass_lightest()
    via_seesaw = m1_active / nu.seesaw_ratio()
    assert direct == pytest.approx(via_seesaw, rel=1e-9)


# ---------------------------------------------------------------------------
# Substrate breakdown text
# ---------------------------------------------------------------------------

def test_substrate_breakdown_runs():
    """substrate_breakdown() returns a non-empty string with Paper 20 markers."""
    s = nu.substrate_breakdown()
    assert isinstance(s, str)
    assert len(s) > 200
    # The breakdown should mention Paper 20's signature concepts
    for marker in ("K_8", "Wilson", "edge partition", "Spin(8)", "δ_CP", "Z_3"):
        assert marker in s, f"missing marker: {marker!r}"


def test_substrate_breakdown_lists_partition():
    """The substrate breakdown lists the (6, 3, 12, 1, 6) K_8 partition."""
    s = nu.substrate_breakdown()
    assert "(6, 3, 12, 1, 6)" in s


# ---------------------------------------------------------------------------
# Module API contract
# ---------------------------------------------------------------------------

def test_module_exports_paper20_predictions():
    """All Paper 20 predictions are exposed in the module __all__."""
    must_export = {
        "active_masses",
        "active_mass_lightest",
        "sum_active_masses",
        "sterile_masses",
        "sterile_active_mixing_sq",
        "seesaw_ratio",
        "pmns_angles_leading_order",
        "delta_cp_winding",
        "jarlskog_invariant",
        "wilson_mass_eV",
        "substrate_breakdown",
        "PMNSAngles",
        "DELTA_M_SQ_21_eV2",
        "DELTA_M_SQ_31_eV2",
    }
    for name in must_export:
        assert hasattr(nu, name), f"missing module export: {name}"
        assert name in nu.__all__, f"missing from __all__: {name}"
