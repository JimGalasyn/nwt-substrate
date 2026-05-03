"""Tests for the nwt.gravity shim."""

import pytest

import nwt_substrate.gravity as grav


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_G_NEWTON_codata():
    assert grav.G_NEWTON_SI == pytest.approx(6.6743e-11, rel=1e-4)


def test_M_PLANCK_GeV():
    assert grav.M_PLANCK_GEV == pytest.approx(1.221e19, rel=1e-3)


def test_M_PLANCK_REDUCED_GeV():
    """M̃_Pl = M_Pl / √(8π) ≈ 2.435 × 10¹⁸ GeV."""
    import math
    assert grav.M_PLANCK_REDUCED_GEV == pytest.approx(
        grav.M_PLANCK_GEV / math.sqrt(8 * math.pi), rel=1e-9
    )


# ---------------------------------------------------------------------------
# Substrate-derived m_e / M_Pl  (Paper 14 LO, Paper 17 NLO/NNLO)
# ---------------------------------------------------------------------------

def test_m_e_over_M_Pl_LO_within_0_3_percent():
    """Paper 14 LO: m_e/M_Pl = (8/7) α^(21/2)  →  -0.12% from CODATA."""
    pred = grav.m_e_over_M_Pl_LO()
    obs = grav.m_e_over_M_Pl_observed()
    assert abs(pred / obs - 1.0) < 0.003


def test_m_e_over_M_Pl_NLO_within_0_05_percent():
    """Paper 17 NLO: × (1 + α/7)  →  ~0.01% from CODATA."""
    pred = grav.m_e_over_M_Pl_NLO()
    obs = grav.m_e_over_M_Pl_observed()
    assert abs(pred / obs - 1.0) < 5e-4


def test_m_e_over_M_Pl_NNLO_within_50_ppm():
    """Paper 17 NNLO: × (1 + α/7 + 3α²)  →  ~14 ppm from CODATA."""
    pred = grav.m_e_over_M_Pl_NNLO()
    obs = grav.m_e_over_M_Pl_observed()
    assert abs(pred / obs - 1.0) < 5e-5


def test_m_e_over_M_Pl_progression():
    """Each correction order should bring us closer to CODATA."""
    obs = grav.m_e_over_M_Pl_observed()
    err_LO   = abs(grav.m_e_over_M_Pl_LO()   / obs - 1.0)
    err_NLO  = abs(grav.m_e_over_M_Pl_NLO()  / obs - 1.0)
    err_NNLO = abs(grav.m_e_over_M_Pl_NNLO() / obs - 1.0)
    assert err_LO > err_NLO > err_NNLO


# ---------------------------------------------------------------------------
# Substrate-derived G in SI
# ---------------------------------------------------------------------------

def test_G_substrate_LO_within_0_3_percent_of_codata():
    """Paper 14 G_LO matches CODATA to -0.24%."""
    G_pred = grav.G_substrate_LO_SI()
    assert abs(G_pred / grav.G_NEWTON_SI - 1.0) < 0.003


def test_G_substrate_NNLO_within_50_ppm_of_codata():
    """Paper 17 G_NNLO matches CODATA to ~30 ppm."""
    G_pred = grav.G_substrate_SI()
    assert abs(G_pred / grav.G_NEWTON_SI - 1.0) < 5e-5


def test_G_substrate_breakdown_mentions_substrate_factors():
    """The breakdown should explicitly identify each factor with its
    substrate origin."""
    breakdown = grav.G_substrate_breakdown()
    assert "8/7" in breakdown
    assert "α^(21/2)" in breakdown or "21/2" in breakdown
    assert "K_7" in breakdown
    assert "Spin(7)" in breakdown
    assert "rank(so(7))" in breakdown or "21/7" in breakdown


# ---------------------------------------------------------------------------
# Black-hole observables
# ---------------------------------------------------------------------------

def test_schwarzschild_radius_proton():
    """r_s = 2 G M_p / c² ≈ 2.5 × 10⁻⁵⁴ m for the proton."""
    M_p = grav.M_PROTON_KG
    r_s = grav.schwarzschild_radius_m(M_p)
    assert 2e-54 < r_s < 3e-54


def test_schwarzschild_radius_sun():
    """The Sun's Schwarzschild radius is ~3 km."""
    r_s = grav.schwarzschild_radius_m(grav.PARTICLE_MASSES_KG["Sun"])
    assert 2.0e3 < r_s < 4.0e3


def test_hawking_temperature_proton_is_extreme():
    """T_H(proton-mass BH) = ℏc³/(8π G M_p k_B) ≈ 7 × 10⁴⁹ K -- absurdly
    hot because the proton mass is so far below the Planck mass."""
    T = grav.hawking_temperature_K(grav.M_PROTON_KG)
    assert 1e49 < T < 1e50


def test_hawking_temperature_planck_mass_is_planck_temperature():
    """T_H(M_Pl) = ℏc³/(8π G M_Pl k_B) ≈ T_Planck/(8π) ≈ 5.6 × 10³⁰ K."""
    T = grav.hawking_temperature_K(grav.PARTICLE_MASSES_KG["M_Planck"])
    # Planck temperature T_Pl = √(ℏc⁵/Gk_B²) ≈ 1.4 × 10³² K, so T_H ≈ T_Pl/(8π)
    assert 1e30 < T < 1e32


def test_hawking_temperature_sun_is_micro_kelvin():
    """T_H(Sun) ≈ 6 × 10⁻⁸ K (extremely cold)."""
    T = grav.hawking_temperature_K(grav.PARTICLE_MASSES_KG["Sun"])
    assert 1e-8 < T < 1e-7


def test_evaporation_time_solar_mass_BH_exceeds_age_of_universe():
    """Solar-mass BH evaporates in ~10⁶⁷ years; age of universe ~10¹⁰ years."""
    t_evap = grav.evaporation_time_s(grav.PARTICLE_MASSES_KG["Sun"])
    age_universe_s = 4.4e17
    assert t_evap / age_universe_s > 1e40


def test_black_hole_summary_proton():
    bh = grav.black_hole_summary("proton")
    assert bh["particle"] == "proton"
    assert "r_s_m" in bh
    assert "T_H_K" in bh
    assert "L_H_W" in bh
    assert "t_evap_s" in bh


def test_black_hole_summary_with_substrate_G():
    """Using substrate G should give very nearly the same result as CODATA G
    (since substrate G matches CODATA to ~30 ppm)."""
    bh_codata = grav.black_hole_summary("proton", use_substrate_G=False)
    bh_substrate = grav.black_hole_summary("proton", use_substrate_G=True)
    rel_diff = abs(bh_substrate["r_s_m"] / bh_codata["r_s_m"] - 1.0)
    assert rel_diff < 1e-3


def test_black_hole_summary_unknown_particle_raises():
    with pytest.raises(KeyError):
        grav.black_hole_summary("xyzzy")


# ---------------------------------------------------------------------------
# Cross-shim consistency
# ---------------------------------------------------------------------------

def test_K7_cycles_match_string_K7_torus():
    """The K_7 graph used by gravity (21 edges = α^21 power, 7 vertices = NLO
    correction denominator) is the same K_7 in nwt.string."""
    import nwt_substrate.string as st
    # K_7 has 21 edges
    assert st.K7_torus.n_cycles == 21
    # And the (8/7) prefactor uses dim S = 8, dim V = 7 of Spin(7)
    # which lives in nwt.string.holonomy too
    h = st.holonomy_lookup("Spin(7)")
    assert h.target_dim == 8


# ---------------------------------------------------------------------------
# Einstein equations / UV-IR bridge / Sakharov / Schwarzschild vacuum
# ---------------------------------------------------------------------------

def test_M_Pl_over_m_e_squared_substrate_matches_codata():
    """The 10⁴⁵ Planck-electron hierarchy matches CODATA to ~30 ppm via
    α⁻²¹ × (8/7)⁻² × (1 + α/7 + 3α²)⁻²."""
    pred = grav.M_Pl_over_m_e_squared_substrate()
    obs = grav.M_Pl_over_m_e_squared_observed()
    assert abs(pred / obs - 1.0) < 5e-5


def test_UV_IR_bridge_breakdown_mentions_alpha_minus_21():
    """The breakdown should articulate 10⁴⁵ = α⁻²¹ × bracket⁻² explicitly."""
    text = grav.UV_IR_bridge_breakdown()
    assert "α⁻²¹" in text
    assert "K_7" in text
    assert "Spin(7)" in text


def test_einstein_hilbert_coefficient_is_1_over_16piG_natural():
    """The EH action prefactor 1/(16π G) in natural units."""
    import math
    coeff = grav.einstein_hilbert_coefficient()
    expected = grav.M_PLANCK_GEV ** 2 / (16.0 * math.pi)
    assert coeff == pytest.approx(expected, rel=1e-9)


def test_sakharov_route_summary_mentions_Paper_18():
    text = grav.sakharov_route_summary()
    assert "Sakharov" in text
    # Should reference the curvature expansion EH + R² + ... + Λ_cc
    assert "(1/16πG)" in text
    assert "Λ_cc" in text


def test_schwarzschild_is_vacuum_solution():
    """Symbolic verification that all 16 components of R_μν vanish for
    the Schwarzschild metric (Paper 18 G6 reproduced in-library via sympy)."""
    result = grav.verify_schwarzschild_vacuum_symbolic()
    assert result["vacuum"] is True
    R = result["R_munu"]
    assert R.shape == (4, 4)
    for mu in range(4):
        for nu in range(4):
            assert R[mu, nu] == 0


def test_linearized_graviton_propagator_text_mentions_spin2_projector():
    text = grav.linearized_graviton_propagator_text()
    assert "spin-2" in text or "P^{μν,αβ}" in text
    assert "8π G" in text or "8πG" in text or "8π G" in text
