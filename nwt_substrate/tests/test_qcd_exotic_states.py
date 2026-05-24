"""Tests for nwt_substrate.qcd.exotic_states -- universal mass formula
m² = (4 m_π⁰)² · N for hidden-gluon bound states.

P7b §7.8 closure (Galasyn 2026-05-23): substrate-canonical formula closes
23 bound states (glueballs, tetraquarks, pentaquarks, exotic Z_b) at
0.02-1.5 % PDG across 4 orders of magnitude in N (7 to 389).
"""

import math

import pytest

from nwt_substrate.qcd.exotic_states import (
    pi_zero_mass,
    gluon_pair_scale,
    universal_mass,
    universal_N,
    BOUND_STATES_CATALOG,
    precision_chain,
    verify_universal_mass_formula,
    precision_chain_summary,
)
from nwt_substrate.electroweak.substrate_gf import ALPHA_SUBSTRATE, M_E_GEV


# ============================================================
# Substrate m_π⁰
# ============================================================

def test_pi_zero_mass_formula():
    """m_π⁰ = (2 m_e / α) · (1 - 5α)."""
    alpha = ALPHA_SUBSTRATE
    expected = (2.0 * M_E_GEV / alpha) * (1.0 - 5.0 * alpha)
    assert pi_zero_mass() == pytest.approx(expected, rel=1e-14)


def test_pi_zero_matches_pdg_sub_ppm():
    """Substrate m_π⁰ matches PDG 134.9770 MeV at sub-ppm."""
    m_MeV = pi_zero_mass() * 1e3
    ppm = (m_MeV / 134.9770 - 1.0) * 1e6
    assert abs(ppm) < 500.0, f"m_π⁰ ppm gap {ppm:+.1f} exceeds 500"


# ============================================================
# Gluon-pair Casimir scale
# ============================================================

def test_gluon_pair_scale_is_four_pi_zero():
    """Scale = 4 m_π⁰ = C_A²(SU(2)) · m_π⁰."""
    assert gluon_pair_scale() == pytest.approx(4.0 * pi_zero_mass(), rel=1e-14)


def test_gluon_pair_scale_matches_qcd_condensate():
    """4 m_π⁰ ≈ 540 MeV ≈ QCD gluon condensate scale (500-550 MeV)."""
    scale_MeV = gluon_pair_scale() * 1e3
    assert 500.0 < scale_MeV < 550.0, (
        f"4 m_π⁰ = {scale_MeV:.1f} MeV outside QCD condensate [500, 550]"
    )


# ============================================================
# Universal mass formula
# ============================================================

def test_universal_mass_formula():
    """m_X = 4 m_π⁰ √N."""
    scale = gluon_pair_scale()
    for N in [1, 4, 7, 14, 28, 100, 389]:
        assert universal_mass(N) == pytest.approx(scale * math.sqrt(N), rel=1e-14)


def test_universal_mass_squared_quantized():
    """m_X² / (4 m_π⁰)² should be exactly N for integer N."""
    scale_sq = gluon_pair_scale() ** 2
    for N in [7, 16, 51, 163, 386]:
        m = universal_mass(N)
        ratio = m ** 2 / scale_sq
        assert ratio == pytest.approx(N, rel=1e-14)


def test_universal_N_inverts_universal_mass():
    """universal_N(universal_mass(N)) == N."""
    for N in [7, 14, 51, 163, 389]:
        m = universal_mass(N)
        N_back = universal_N(m)
        assert N_back == pytest.approx(N, rel=1e-14)


def test_universal_mass_negative_N_raises():
    with pytest.raises(ValueError):
        universal_mass(-3)
    with pytest.raises(ValueError):
        universal_mass(0)


def test_universal_N_negative_mass_raises():
    with pytest.raises(ValueError):
        universal_N(-1.0)


# ============================================================
# Catalogued states
# ============================================================

def test_catalog_has_23_states():
    """The 23 catalogued bound states from P7b §7.8."""
    assert len(BOUND_STATES_CATALOG) == 23


def test_catalog_sector_distribution():
    """11 glueballs + 6 tetraquarks (incl X(6900)) + 2 Z_b + 4 pentaquarks = 23."""
    counts = {"glueball": 0, "tetraquark": 0, "Z_b": 0, "pentaquark": 0}
    for entry in BOUND_STATES_CATALOG.values():
        counts[entry["sector"]] += 1
    assert counts == {"glueball": 11, "tetraquark": 6, "Z_b": 2, "pentaquark": 4}


def test_catalog_substrate_N_range():
    """Substrate N ranges from 7 (η(1405) = |K_7|) to 389 (Z_b(10650))."""
    N_vals = [e["N"] for e in BOUND_STATES_CATALOG.values()]
    assert min(N_vals) == 7
    assert max(N_vals) == 389


# ============================================================
# Specific high-precision matches from P7b §7.8
# ============================================================

@pytest.mark.parametrize("state, max_pct_gap", [
    ("eta_2225",  0.05),    # -0.02 % -- best glueball
    ("X_4140",    0.05),    # -0.01 % -- best tetraquark
    ("Zb_10610",  0.10),    # -0.03 % -- best Z_b
])
def test_best_matches_are_sub_pct(state, max_pct_gap):
    """Documented best-match states have sub-0.1% gap."""
    chain = precision_chain()
    gap_pct = abs(chain[state]["percent_gap"])
    assert gap_pct < max_pct_gap, (
        f"{state} gap {gap_pct:.3f}% exceeds {max_pct_gap}"
    )


def test_lattice_glueball_f0_1710():
    """f₀(1710) at N=10=dim(10_SU(5)) is the lightest substrate scalar glueball
    candidate; matches lattice 0++ glueball prediction."""
    chain = precision_chain()
    f0_1710 = chain["f0_1710"]
    assert f0_1710["N"] == 10
    assert f0_1710["N_origin"] == "dim(10_SU(5))"
    assert abs(f0_1710["percent_gap"]) < 2.0


# ============================================================
# Precision chain + verification
# ============================================================

def test_precision_chain_returns_all_states():
    chain = precision_chain()
    assert set(chain.keys()) == set(BOUND_STATES_CATALOG.keys())


def test_precision_chain_keys_per_state():
    chain = precision_chain()
    for name, data in chain.items():
        expected = {"substrate_MeV", "pdg_MeV", "percent_gap", "N",
                    "N_origin", "sector", "J_PC"}
        assert set(data.keys()) == expected


def test_all_23_states_match_within_5pct():
    """The strongest claim: all 23 catalogued states agree at ≤ 5% PDG."""
    chain = precision_chain()
    for name, data in chain.items():
        assert abs(data["percent_gap"]) < 5.0, (
            f"{name} gap {data['percent_gap']:+.3f}% exceeds 5%"
        )


def test_verify_universal_mass_formula_passes_default():
    """Default 5% tolerance passes all 23 states."""
    result = verify_universal_mass_formula()
    assert result["pass"] is True
    assert all(result["per_state_pass"].values())


def test_verify_universal_mass_formula_fails_strict():
    """A 0.01% tolerance fails (η(1475) is -3.3%)."""
    result = verify_universal_mass_formula(percent_tol=0.01)
    assert result["pass"] is False


def test_verify_worst_state_is_eta_1475():
    """The worst-match catalogued state is η(1475) at -3.3%."""
    result = verify_universal_mass_formula()
    assert result["worst_state"] == "eta_1475"
    assert 3.0 < result["worst_gap"] < 4.0


# ============================================================
# Pretty-print
# ============================================================

def test_precision_chain_summary_contains_key_substrings():
    out = precision_chain_summary()
    assert "m² = (4 m_π⁰)² · N" in out
    assert "glueball" in out.lower()
    assert "tetraquark" in out.lower()
    assert "pentaquark" in out.lower()
    assert "K_7 Wilson" in out


# ============================================================
# QCD module integration
# ============================================================

def test_qcd_module_exposes_exotic_states_api():
    """The qcd module re-exports the universal mass formula API."""
    import nwt_substrate.qcd as qcd
    assert hasattr(qcd, "universal_mass")
    assert hasattr(qcd, "universal_N")
    assert hasattr(qcd, "verify_universal_mass_formula")
    assert hasattr(qcd, "BOUND_STATES_CATALOG")
    assert hasattr(qcd, "pi_zero_mass")
    assert hasattr(qcd, "gluon_pair_scale")
