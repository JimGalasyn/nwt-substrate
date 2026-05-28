"""Unit tests for nwt_substrate.dark_sector.wimp_98gev."""

import math
import pytest

from nwt_substrate.dark_sector import (
    WIMP_98GeV,
    sigma_si_higgs_portal,
    lhc_production_cross_section,
    compare_with_lz_limit,
    predict_all,
)
from nwt_substrate.isa.constants import ALPHA_SUBSTRATE, M_PL_GEV


def test_wimp_mass_is_alpha8_mpl():
    """98 GeV WIMP mass = α^8 · M_Pl (K_8 tower N_e=16)."""
    w = WIMP_98GeV()
    assert w.mass_gev == pytest.approx(ALPHA_SUBSTRATE ** 8 * M_PL_GEV)
    assert 90 < w.mass_gev < 110             # near the 98 GeV target
    assert w.n_e == 16


def test_g_eff_is_g_times_v():
    """g_eff = g_Hχχ · v_Higgs."""
    g = 0.01
    w = WIMP_98GeV(g_hxx=g)
    assert w.g_eff == pytest.approx(g * 246.22)


def test_sigma_si_excluded_at_g_alpha2():
    """At g_Hχχ = α², σ_SI > LZ limit by O(10⁶)."""
    w = WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE ** 2)
    sigma = sigma_si_higgs_portal(w)
    cmp = compare_with_lz_limit(sigma)
    assert cmp["verdict"] == "EXCLUDED"
    assert cmp["ratio_pred_over_limit"] > 1e5


def test_sigma_si_allowed_at_g_alpha4():
    """At g_Hχχ = α^4, σ_SI well below LZ-2024 limit."""
    w = WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE ** 4)
    sigma = sigma_si_higgs_portal(w)
    cmp = compare_with_lz_limit(sigma)
    assert cmp["verdict"] == "ALLOWED"
    # σ_SI ~ 10⁻⁵⁰ cm², below LZ-G3 sensitivity ~10⁻⁴⁸ but possibly reachable
    assert 1e-52 < sigma < 1e-48


def test_sigma_si_scales_as_g_squared():
    """σ_SI ∝ g_eff² → σ_SI(α^4) / σ_SI(α^8) = (α^4 / α^8)² = 1/α^8."""
    w1 = WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE ** 4)
    w2 = WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE ** 8)
    ratio = sigma_si_higgs_portal(w1) / sigma_si_higgs_portal(w2)
    expected = (ALPHA_SUBSTRATE ** 4 / ALPHA_SUBSTRATE ** 8) ** 2
    assert ratio == pytest.approx(expected, rel=1e-9)


def test_lhc_cross_section_nonnegative():
    """LHC off-shell cross section should be nonneg + finite."""
    w = WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE)
    s = lhc_production_cross_section(w)
    assert s >= 0
    assert math.isfinite(s)


def test_predict_all_returns_dict():
    """The bundled predictor returns structured output."""
    result = predict_all(WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE ** 4), show_table=False)
    assert "wimp" in result
    assert "direct_detection" in result
    assert "lhc" in result
    assert result["wimp"]["mass_gev"] > 0
    assert result["direct_detection"]["sigma_si_cm2"] > 0
    assert result["direct_detection"]["lz_comparison"]["verdict"] in ("ALLOWED", "EXCLUDED")
