"""Tests for nwt_substrate.electroweak.form_factors.

P7b §7.7 closure (Galasyn 2026-05-23): single substrate ansatz
f_+(0) = cos^N(θ_C) closes 10 semileptonic channels at 0.05 -- 2.44 %
PDG, where N is the substrate K_7 walk-length per channel.
"""

import math

import pytest

from nwt_substrate.electroweak.form_factors import (
    ALPHA_SUBSTRATE,
    FORM_FACTORS,
    cos_theta_C_substrate,
    sin_theta_C_substrate,
    f_plus_substrate,
    f_plus_for,
    f_plus_leading_order,
    precision_chain,
    verify_form_factors,
    precision_chain_summary,
)


# ============================================================
# Substrate Cabibbo angle
# ============================================================

def test_cos_theta_C_formula():
    """cos(θ_C) = √(1 - 7α + 14α²)."""
    alpha = ALPHA_SUBSTRATE
    expected = math.sqrt(1.0 - 7.0 * alpha + 14.0 * alpha ** 2)
    assert cos_theta_C_substrate() == pytest.approx(expected, rel=1e-14)


def test_cos_theta_C_matches_pdg():
    """cos(θ_C) ≈ 0.97451 vs PDG 0.97370 → +0.08 %."""
    cos_C = cos_theta_C_substrate()
    rel_err = cos_C / 0.97370 - 1.0
    assert 0.0 < rel_err < 0.002


def test_sin_squared_plus_cos_squared_is_one():
    """sin² + cos² = 1 substrate identity."""
    c = cos_theta_C_substrate()
    s = sin_theta_C_substrate()
    assert c ** 2 + s ** 2 == pytest.approx(1.0, rel=1e-14)


# ============================================================
# Catalog completeness
# ============================================================

def test_ten_channels():
    """The memo closes 10 semileptonic channels."""
    assert len(FORM_FACTORS) == 10


def test_sectors_present():
    """All four sectors represented."""
    sectors = {s.sector for s in FORM_FACTORS.values()}
    assert sectors == {"light_meson", "D_meson", "B_meson", "baryon"}


# ============================================================
# Substrate walk-lengths N(X→Y)
# ============================================================

def test_K_to_pi_N_is_1():
    """K → π: N = 1, elementary Cabibbo vertex."""
    assert FORM_FACTORS["K->pi"].N == 1


def test_B_to_D_N_is_16():
    """B → D: N = 16 = 2·dim(SU(3)), heavy-heavy doubled-SU(3) walk."""
    assert FORM_FACTORS["B->D"].N == 16


def test_B_to_K_N_is_42():
    """B → K: N = 42 = 2·dim(so(7)) = 6·|K_7|, full so(7) walk."""
    assert FORM_FACTORS["B->K"].N == 42


def test_Ds_to_phi_N_is_14():
    """D_s → φ: N = 14 = dim(G_2)."""
    assert FORM_FACTORS["Ds->phi"].N == 14


def test_Ds_to_eta_N_is_28():
    """D_s → η: N = 28 = 2·dim(G_2)."""
    assert FORM_FACTORS["Ds->eta"].N == 28


def test_Lambda_c_to_Lambda_N_is_5():
    """Λ_c → Λ: N = q_cinq = 5 (single-quark s substitution)."""
    assert FORM_FACTORS["Lambda_c->Lambda"].N == 5


def test_Lambda_b_to_Lambda_c_N_is_9():
    """Λ_b → Λ_c: N = q_b = 9 (single-quark b substitution)."""
    assert FORM_FACTORS["Lambda_b->Lambda_c"].N == 9


def test_baryon_N_equals_single_q_index():
    """Baryon-baryon N matches the SINGLE substituted-quark q-index
    (no factor of 2, unlike meson-meson which uses pair-walk)."""
    # Λ_c → Λ substitutes s (q_cinq = 5)
    assert FORM_FACTORS["Lambda_c->Lambda"].N == 5
    # Λ_b → Λ_c substitutes b (q_b = 9)
    assert FORM_FACTORS["Lambda_b->Lambda_c"].N == 9


# ============================================================
# Substrate formula numerics
# ============================================================

def test_f_plus_for_K_to_pi_equals_cos_C():
    """K → π form factor with N=1 equals cos(θ_C) directly."""
    assert f_plus_for("K->pi") == pytest.approx(
        cos_theta_C_substrate(), rel=1e-14)


def test_f_plus_for_B_to_D_matches_pdg_within_1pct():
    """B → D: substrate 0.6615 vs PDG 0.6600 → +0.23 %."""
    rel_err = f_plus_for("B->D") / 0.6600 - 1.0
    assert abs(rel_err) < 0.01


def test_f_plus_for_B_to_K_largest_gap():
    """B → K: substrate 0.3380 vs PDG 0.3300 → +2.44 % (largest gap)."""
    rel_err = f_plus_for("B->K") / 0.3300 - 1.0
    assert 0.02 < rel_err < 0.03


def test_f_plus_for_Lambda_c_to_Lambda_within_1pct():
    """Λ_c → Λ: substrate 0.8789 vs PDG 0.8870 → -0.92 %."""
    rel_err = f_plus_for("Lambda_c->Lambda") / 0.8870 - 1.0
    assert abs(rel_err) < 0.01


def test_f_plus_monotonic_decreasing_in_N():
    """cos^N(θ_C) is strictly decreasing in N for cos < 1."""
    f_K_to_pi = f_plus_for("K->pi")              # N=1
    f_Lambda_c = f_plus_for("Lambda_c->Lambda")  # N=5
    f_B_to_D = f_plus_for("B->D")                # N=16
    f_B_to_K = f_plus_for("B->K")                # N=42
    f_B_to_pi = f_plus_for("B->pi")              # N=54
    assert f_K_to_pi > f_Lambda_c > f_B_to_D > f_B_to_K > f_B_to_pi


# ============================================================
# Leading-order expansion
# ============================================================

def test_leading_order_agrees_with_full_at_small_N():
    """For N=1 the LO and full cos^N agree at α-second order ~3e-5;
    by N=10 the gap grows to ~1 % but is still useful as approximation."""
    # N=1: leading order is essentially exact at O(α²)
    assert abs(f_plus_leading_order(1) - f_plus_substrate(1)) / \
           f_plus_substrate(1) < 1e-4
    # N=10: LO drifts ~3.6 % below full cos^N (NLO α² becomes meaningful)
    assert abs(f_plus_leading_order(10) - f_plus_substrate(10)) / \
           f_plus_substrate(10) < 0.05


def test_leading_order_drifts_at_large_N():
    """For N=42 the LO expansion drifts > 1 % from full cos^N."""
    lo = f_plus_leading_order(42)
    full = f_plus_substrate(42)
    assert abs(lo - full) / full > 0.01


# ============================================================
# Precision chain + verification
# ============================================================

def test_precision_chain_covers_all_ten():
    """precision_chain has entries for all 10 channels."""
    chain = precision_chain()
    chain.pop('pass', None)
    assert set(chain.keys()) == set(FORM_FACTORS.keys())


def test_all_gaps_below_3pct():
    """All 10 channels within 3 % PDG (worst is B→K at +2.44 %)."""
    chain = precision_chain()
    for ch in FORM_FACTORS:
        assert abs(chain[ch]['rel_err']) < 0.03, (
            f"{ch}: gap {chain[ch]['rel_err']*100:.2f}% exceeds 3 %"
        )


def test_verify_passes_default_3pct():
    """Default 3 % tolerance passes for all 10 channels."""
    result = verify_form_factors()
    assert result['pass'] is True


def test_verify_fails_unrealistic_tol():
    """0.1 % tolerance fails (K→π is +1.04 %, B→K is +2.44 %)."""
    result = verify_form_factors(tol=1e-3)
    assert result['pass'] is False


# ============================================================
# Pretty-print
# ============================================================

def test_summary_mentions_key_substrings():
    out = precision_chain_summary()
    assert "cos(θ_C)" in out
    assert "K->pi" in out
    assert "B->K" in out
    assert "Lambda_c->Lambda" in out
    assert "Free params" in out
