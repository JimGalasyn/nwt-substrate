"""Substrate decay constants from `particles.decay_constants`.

Single canonical home for the three sectors (P7b §2-3 light, §7.5 heavy,
§7.6 vector + B_c). Tests cover scales (f_π, m_τ, S_V), per-sector closed
forms, name-based lookups, ratio diagnostics, precision chains, the
unified verify, and the summary printer.
"""

from __future__ import annotations

import math
from fractions import Fraction

import pytest

from nwt_substrate.particles.decay_constants import (
    HEAVY_PSEUDOSCALARS,
    HeavyMesonSpec,
    LIGHT_PSEUDOSCALARS,
    LightPseudoscalarSpec,
    R_S_NONSTRANGE,
    R_S_STRANGE,
    SQRT_7_OVER_4,
    VECTOR_MESONS,
    VectorMesonSpec,
    c_ratio_precision_chain,
    cabibbo_scale_fX,
    f_pi_substrate,
    fX_ratio_strange_over_nonstrange,
    fibonacci_anomaly_fX,
    heavy_meson_fX,
    heavy_meson_fX_for,
    heavy_precision_chain,
    light_meson_fX_for,
    light_precision_chain,
    m_tau_substrate,
    precision_chain_summary,
    vector_meson_binding_scale,
    vector_meson_fX,
    vector_meson_fX_for,
    vector_precision_chain,
    verify_decay_constants,
    verify_heavy_meson_fX,
    verify_light_meson_fX,
    verify_vector_meson_fX,
)


# ---------------------------------------------------------------
# Substrate scales
# ---------------------------------------------------------------

def test_m_tau_substrate_matches_pdg_sub_percent():
    PDG_M_TAU = 1.77686  # GeV
    assert abs(m_tau_substrate() - PDG_M_TAU) / PDG_M_TAU < 0.01


def test_f_pi_substrate_within_3_pct_of_pdg():
    PDG_F_PI = 92.4  # MeV
    assert abs(f_pi_substrate() * 1e3 - PDG_F_PI) / PDG_F_PI < 0.03


def test_vector_binding_scale_matches_f_rho_m_rho():
    PDG = 0.215 * 0.77526  # f_ρ × m_ρ in GeV²
    assert abs(vector_meson_binding_scale() - PDG) / PDG < 0.04


def test_rs_strange_equals_sqrt_7_over_4():
    assert R_S_STRANGE == pytest.approx(math.sqrt(7.0 / 4.0))
    assert R_S_NONSTRANGE == 1.0
    assert SQRT_7_OVER_4 == R_S_STRANGE


# ---------------------------------------------------------------
# Light sector (K, η, π)
# ---------------------------------------------------------------

def test_cabibbo_scale_fX_formula():
    """f = m · √(7α). At m=1 GeV, substrate α: ~0.225 GeV."""
    f = cabibbo_scale_fX(1.0)
    assert f == pytest.approx(math.sqrt(7.0 * 0.00729927), rel=1e-3)


def test_fibonacci_anomaly_fX_formula():
    """f = m / walk_length^(1/4); default walk_length=5."""
    f5 = fibonacci_anomaly_fX(1.0)
    assert f5 == pytest.approx(1.0 / 5.0 ** 0.25, rel=1e-12)
    f8 = fibonacci_anomaly_fX(1.0, walk_length=8)
    assert f8 == pytest.approx(1.0 / 8.0 ** 0.25, rel=1e-12)


def test_light_meson_fX_for_K_within_2_pct():
    f_K_MeV = light_meson_fX_for("K") * 1e3
    PDG = 110.0
    assert abs(f_K_MeV - PDG) / PDG < 0.02


def test_light_meson_fX_for_pi_within_2_pct():
    f_pi_MeV = light_meson_fX_for("pi") * 1e3
    PDG = 92.4
    assert abs(f_pi_MeV - PDG) / PDG < 0.02


def test_light_meson_fX_for_unknown_raises():
    with pytest.raises(KeyError):
        light_meson_fX_for("nope")


# ---------------------------------------------------------------
# Heavy sector (D, D_s, B, B_s)
# ---------------------------------------------------------------

def test_heavy_meson_fX_D_within_3_pct():
    f_D = heavy_meson_fX_for("D") * 1e3
    PDG = 212.0
    assert abs(f_D - PDG) / PDG < 0.03


def test_heavy_meson_fX_Ds_within_3_pct():
    f_Ds = heavy_meson_fX_for("D_s") * 1e3
    PDG = 248.0
    assert abs(f_Ds - PDG) / PDG < 0.03


def test_heavy_meson_fX_strange_enhancement_is_7_over_4_quarter_power():
    """Strange enhances f_X by (7/4)^(1/4) at fixed N, m."""
    f_s = heavy_meson_fX(N=10, m_X_GeV=2.0, strange=True)
    f_n = heavy_meson_fX(N=10, m_X_GeV=2.0, strange=False)
    assert f_s / f_n == pytest.approx((7.0 / 4.0) ** 0.25, rel=1e-12)


def test_heavy_meson_fX_for_unknown_raises():
    with pytest.raises(KeyError):
        heavy_meson_fX_for("nope")


def test_fX_ratio_strange_over_nonstrange_Ds_over_D():
    """Substrate Ds/D ratio dominated by (7/4)^(1/4) ≈ 1.15."""
    ratio = fX_ratio_strange_over_nonstrange("D_s", "D")
    assert 1.05 < ratio < 1.30  # within a few % either side of 1.15


# ---------------------------------------------------------------
# Vector + B_c sector (11 states)
# ---------------------------------------------------------------

def test_vector_meson_fX_rho_reference_within_4_pct():
    f_rho = vector_meson_fX_for("rho") * 1e3
    PDG = 215.0
    assert abs(f_rho - PDG) / PDG < 0.04


def test_vector_meson_fX_accepts_Fraction_and_float():
    f_int = vector_meson_fX(C=3.0, m_X_GeV=2.0)
    f_frac = vector_meson_fX(C=Fraction(3, 1), m_X_GeV=2.0)
    assert f_int == pytest.approx(f_frac, rel=1e-15)


def test_vector_meson_fX_for_Bc_present():
    f_Bc = vector_meson_fX_for("B_c") * 1e3
    PDG = 427.0
    assert abs(f_Bc - PDG) / PDG < 0.08  # B_c has larger gap; 8 % envelope


def test_vector_meson_fX_for_unknown_raises():
    with pytest.raises(KeyError):
        vector_meson_fX_for("nope")


# ---------------------------------------------------------------
# Catalogs
# ---------------------------------------------------------------

def test_light_catalog_has_pi_K_eta():
    assert set(LIGHT_PSEUDOSCALARS) == {"pi", "K", "eta"}
    for name, spec in LIGHT_PSEUDOSCALARS.items():
        assert isinstance(spec, LightPseudoscalarSpec)
        assert spec.regime in {"cabibbo", "fibonacci"}


def test_heavy_catalog_has_D_Ds_B_Bs():
    assert set(HEAVY_PSEUDOSCALARS) == {"D", "D_s", "B", "B_s"}
    for name, spec in HEAVY_PSEUDOSCALARS.items():
        assert isinstance(spec, HeavyMesonSpec)


def test_vector_catalog_has_11_states_including_Bc():
    assert len(VECTOR_MESONS) == 11
    assert "B_c" in VECTOR_MESONS
    for name, spec in VECTOR_MESONS.items():
        assert isinstance(spec, VectorMesonSpec)
        assert isinstance(spec.C, Fraction)


# ---------------------------------------------------------------
# Per-sector precision chains + verify
# ---------------------------------------------------------------

def test_light_precision_chain_structure():
    chain = light_precision_chain()
    assert set(chain) == set(LIGHT_PSEUDOSCALARS)
    for name, data in chain.items():
        assert "substrate_MeV" in data
        assert "rel_err" in data
        assert "regime" in data


def test_heavy_precision_chain_carries_N_and_strange_flag():
    chain = heavy_precision_chain()
    for name, data in chain.items():
        assert "N" in data
        assert isinstance(data["N"], int)
        assert isinstance(data["strange"], bool)


def test_vector_precision_chain_carries_C_and_is_heavy_flag():
    chain = vector_precision_chain()
    for name, data in chain.items():
        assert isinstance(data["C"], float)
        assert isinstance(data["is_heavy"], bool)


def test_c_ratio_precision_chain_rho_is_unity():
    chain = c_ratio_precision_chain()
    assert chain["rho"]["C_substrate"] == 1.0
    assert chain["rho"]["C_empirical"] == pytest.approx(1.0, rel=1e-15)
    assert chain["rho"]["rel_err"] == pytest.approx(0.0, abs=1e-15)


def test_verify_light_meson_passes_5_pct():
    chain = verify_light_meson_fX()
    assert chain["pass"] is True


def test_verify_heavy_meson_passes_5_pct():
    chain = verify_heavy_meson_fX()
    assert chain["pass"] is True


def test_verify_vector_meson_passes_7_pct():
    chain = verify_vector_meson_fX()
    assert chain["pass"] is True


def test_verify_vector_meson_tight_tol_can_fail():
    """At 1 % tolerance some vector states should fail."""
    chain = verify_vector_meson_fX(tol=0.01)
    assert chain["pass"] is False


# ---------------------------------------------------------------
# Unified verify_decay_constants
# ---------------------------------------------------------------

def test_verify_decay_constants_all_sectors_within_default_tol():
    out = verify_decay_constants()
    assert out["pass"] is True
    assert out["worst_gap"] < 7.0
    assert out["worst_meson"]
    assert set(out) >= {"light", "heavy", "vector_Bc", "pass", "per_meson_pass",
                        "worst_gap", "worst_meson"}


def test_verify_decay_constants_tight_tol_fails_with_named_worst():
    out = verify_decay_constants(percent_tol=0.5)
    assert out["pass"] is False
    assert out["worst_gap"] > 0.5
    assert "/" in out["worst_meson"]  # sector/name slash


# ---------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------

def test_precision_chain_summary_contains_section_headers_and_keys():
    s = precision_chain_summary()
    assert "LIGHT PSEUDOSCALARS" in s
    assert "HEAVY PSEUDOSCALARS" in s
    assert "VECTORS + B_c" in s
    assert "Substrate f_π" in s
    assert "Substrate m_τ" in s
    assert "rho" in s
    assert "B_c" in s
    assert "pi" in s


