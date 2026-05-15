"""Tests for the Kerr / cosmogenesis / disk-vortex extensions to
nwt.gravity (vortex-vision 2026-05-15 results)."""
from __future__ import annotations

import math

import pytest

import nwt_substrate.gravity as grav
from nwt_substrate.isa import (
    ALPHA_NWT,
    KAPPA_MACKEN,
    RANK_SO7,
    H_V_SO7,
    H_COXETER_SO7,
    DIM_OCTONION,
)


SQRT2 = math.sqrt(2.0)
SQRT3 = math.sqrt(3.0)


# ---------------------------------------------------------------------------
# Substrate Coxeter constants — ISA structural identities
# ---------------------------------------------------------------------------

def test_spin7_coxeter_decomposition_8_equals_5_plus_3():
    """8 = h_v + rank = 5 + 3 is the Spin(7) Coxeter decomposition."""
    assert H_V_SO7 + RANK_SO7 == DIM_OCTONION


def test_h_coxeter_is_twice_rank():
    """h(so(7)) = 2 × rank(so(7))."""
    assert H_COXETER_SO7 == 2 * RANK_SO7


def test_alpha_nwt_matches_codata_ppm():
    """α_NWT closed form (Paper 17 trefoil) within ~8 ppm of CODATA α_QED."""
    alpha_codata = 1.0 / 137.035999084
    assert abs(ALPHA_NWT - alpha_codata) / alpha_codata < 1e-5


def test_kappa_macken_closed_form():
    """κ_Macken² × √2 = 1/α (closed-form consistency)."""
    assert KAPPA_MACKEN ** 2 * SQRT2 == pytest.approx(1.0 / ALPHA_NWT, rel=1e-9)


# ---------------------------------------------------------------------------
# Kerr efficiency predictions
# ---------------------------------------------------------------------------

def test_m_irr_extremal_kerr_0p05_percent():
    """M_irr/M = 1 - 3√α - 5α matches GR 1/√2 at < 0.05 % (framework's
    tightest sub-1 % GR prediction)."""
    pred = grav.m_irreducible_over_m_extremal()
    gr = 1.0 / SQRT2
    assert abs(pred - gr) / gr < 5e-4


def test_penrose_extraction_0p1_percent():
    """Penrose 1 - 1/√2 = 3√α + 5α at < 0.1 %.
    Exact complement of M_irr — same substrate form, both sides match."""
    pred = grav.penrose_extraction_fraction()
    gr = 1.0 - 1.0 / SQRT2
    assert abs(pred - gr) / gr < 1e-3


def test_m_irr_and_penrose_sum_to_one():
    """The two halves of the mass partition must sum to 1 by construction."""
    m_irr = grav.m_irreducible_over_m_extremal()
    penrose = grav.penrose_extraction_fraction()
    assert m_irr + penrose == pytest.approx(1.0, rel=1e-12)


def test_bardeen_prograde_2_percent():
    """Bardeen extremal prograde η = 1 - 1/√3 ≈ 5√α at ~1 %."""
    pred = grav.bardeen_prograde_efficiency()
    gr = 1.0 - 1.0 / SQRT3
    assert abs(pred - gr) / gr < 0.02


def test_bardeen_retrograde_5_percent():
    """Retrograde extremal η = 1 - 5/(3√3) ≈ 5α at ~3-5 %."""
    pred = grav.bardeen_retrograde_efficiency()
    gr = 1.0 - 5.0 / (3.0 * SQRT3)
    assert abs(pred - gr) / gr < 0.05


def test_schwarzschild_isco_3_percent():
    """Schwarzschild ISCO η = 1 - √(8/9) ≈ 8α at ~2 %."""
    pred = grav.schwarzschild_isco_efficiency()
    gr = 1.0 - math.sqrt(8.0 / 9.0)
    assert abs(pred - gr) / gr < 0.03


def test_kerr_efficiency_table_keys():
    """Sanity check on the convenience table."""
    table = grav.kerr_efficiency_table()
    for key in (
        "M_irr_over_M_extremal",
        "penrose_extraction",
        "bardeen_prograde",
        "bardeen_retrograde",
        "schwarzschild_isco",
        "em_superradiance_max",
    ):
        assert key in table
        assert "substrate" in table[key]
        assert "gr" in table[key]


# ---------------------------------------------------------------------------
# Cosmogenesis: f_J, κ ratio, Thorne equilibrium
# ---------------------------------------------------------------------------

def test_f_J_cosmogenic_full_vs_leading():
    """(1-√α)³ expands to ≈ 1 - 3√α + 3α - α^(3/2). Full form > leading at
    O(α) and beyond."""
    full = grav.f_J_cosmogenic()
    leading = grav.f_J_cosmogenic_leading()
    assert full > leading
    assert full - leading == pytest.approx(
        3 * ALPHA_NWT - ALPHA_NWT ** 1.5, rel=1e-12
    )


def test_f_J_in_expected_range():
    """f_J = (1-√α)³ ≈ 0.765 (full); 1 - 3√α ≈ 0.744 (leading)."""
    assert grav.f_J_cosmogenic() == pytest.approx(0.765, abs=0.002)
    assert grav.f_J_cosmogenic_leading() == pytest.approx(0.744, abs=0.002)


def test_kappa_parent_over_daughter_is_3_times_one_plus_sqrt_alpha():
    """Substrate form: 3(1 + √α) ≈ 3.256."""
    pred = grav.kappa_parent_over_daughter()
    expected = RANK_SO7 * (1.0 + math.sqrt(ALPHA_NWT))
    assert pred == pytest.approx(expected, rel=1e-12)


def test_thorne_a_star_within_0p005_of_observed():
    """Cosmogenic-viability spin a*_eq matches Thorne 1974 thin-disk a* ≈ 0.998
    within 0.005."""
    a_star = grav.thorne_a_star_equilibrium()
    assert abs(a_star - 0.998) < 0.005


def test_mode_overlap_equal_kappa_is_unit_transmission():
    """When κ_parent == κ_daughter, every mode passes losslessly: T = 1."""
    for p, q in [(1, 0), (2, 1), (3, 2)]:
        T = grav.mode_overlap(p, q, KAPPA_MACKEN, KAPPA_MACKEN)
        assert T == pytest.approx(1.0, rel=1e-12)


def test_mode_overlap_unequal_kappa_below_unity():
    """When κ's differ and p > 0 (so the (p/κ)² term matters), T < 1."""
    T = grav.mode_overlap(1, 0, kappa_parent=20.0, kappa_daughter=KAPPA_MACKEN)
    assert 0.0 < T < 1.0


# ---------------------------------------------------------------------------
# Disk-ergosurface matching
# ---------------------------------------------------------------------------

def test_kappa_disk_over_parent_matches_qpo_median():
    """8(1 - √α) ≈ 7.32 matches QPO-implied κ_disk/κ_parent ≈ 7.28 at ~0.5 %."""
    pred = grav.kappa_disk_over_parent()
    qpo_median = 7.28
    assert abs(pred - qpo_median) / qpo_median < 0.01


def test_kappa_disk_over_daughter_telescopes_to_24_one_minus_alpha():
    """Compounded ratio 8(1-√α) × 3(1+√α) = 24(1 - α)."""
    pred = grav.kappa_disk_over_daughter()
    expected = 24.0 * (1.0 - ALPHA_NWT)
    assert pred == pytest.approx(expected, rel=1e-12)


def test_coxeter_decomposition_sums_to_eight():
    """Z_3 decomposition: 1 + 1 + 3 + 3 = 8 = DIM_OCTONION."""
    decomp = grav.coxeter_decomposition()
    parts = decomp["z3_decomp"]
    assert parts["identity"] + parts["z3_fixed_polar"] + \
           parts["z3_orbit_A"] + parts["z3_orbit_B"] == DIM_OCTONION


def test_substrate_energy_budget_5_plus_3_equals_8():
    """5√α (steady state) + 3√α (cosmogenic) = 8√α (total dissipation)."""
    budget = grav.substrate_energy_budget()
    assert budget["steady_state_loss"] + budget["cosmogenic_loss"] == \
           pytest.approx(budget["total_dissipation"], rel=1e-12)


# ---------------------------------------------------------------------------
# QPO torus eigenmodes
# ---------------------------------------------------------------------------

def test_twin_peak_ratio_exactly_3_over_2():
    """f(0,3) / f(0,2) = 3/2 exactly under the substrate-vortex picture."""
    assert grav.twin_peak_ratio() == pytest.approx(1.5, rel=1e-12)


def test_xte_j1550_564_92hz_qpo_at_zero_one_mode():
    """XTE J1550-564's mysterious 92 Hz HFQPO matches the predicted
    (p=0, q=1) substrate-vortex mode at 91.9 Hz (deviation < 1 %).

    The (0, 2) HFQPO at 184 Hz fixes r_disk; the (0, 1) prediction
    follows by f(0, 1) = (1/2) × f(0, 2) = 92 Hz."""
    r_disk = grav.predict_r_disk_from_qpo(184.0, 0, 2)
    f_01 = grav.torus_eigenmode_freq(0, 1, r_disk)
    assert abs(f_01 - 92.0) / 92.0 < 0.01


def test_qpo_mode_table_size():
    r_disk = grav.predict_r_disk_from_qpo(184.0, 0, 2)
    table = grav.qpo_mode_table(r_disk, p_max=2, q_max=3)
    # (p, q) ∈ [0..2] × [0..3] minus (0, 0) = 3 × 4 - 1 = 11 modes
    assert len(table) == 11


# ---------------------------------------------------------------------------
# Healing length
# ---------------------------------------------------------------------------

def test_xi_substrate_is_reduced_electron_compton():
    """ξ_substrate = ℏ / (m_e c) ≈ 3.86 × 10⁻¹³ m."""
    assert grav.XI_SUBSTRATE_M == pytest.approx(3.86e-13, rel=1e-2)


def test_k_chern_simons_paper_17_value():
    """Paper 17 Chern-Simons level k_CS ≈ 31.73 from κ_Macken."""
    k_cs = grav.k_chern_simons()
    assert k_cs == pytest.approx(31.73, abs=0.1)


def test_xi_cosmo_in_galactic_halo_band():
    """ξ_cosmo from k_CS-step amplification falls between dwarf-galaxy and
    cluster scales (the mesoscale band)."""
    xi_c = grav.xi_cosmo_m()
    parsec = 3.0857e16
    kpc = 1e3 * parsec
    Mpc = 1e6 * parsec
    # Allow a generous 1 kpc – 100 Mpc band, since the exponential
    # sensitivity to k_CS makes the precise value vary by orders of
    # magnitude. The key prediction is the qualitative regime.
    assert 1.0 * kpc < xi_c < 100.0 * Mpc
