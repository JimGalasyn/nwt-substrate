"""Tests for nwt_substrate.electroweak.light_meson_decay_constants.

P7b §2-3 closure: light pseudoscalar f_X via two distinct substrate laws:
  - Cabibbo scale f_X = m_X · √(7α) for K, η (ordinary light pseudoscalars)
  - Fibonacci anomaly f_π = m_π / 5^(1/4) (Goldstone chiral special case)
"""

import math

import pytest

from nwt_substrate.electroweak.light_meson_decay_constants import (
    ALPHA_SUBSTRATE,
    LIGHT_PSEUDOSCALARS,
    cabibbo_scale_fX,
    fibonacci_anomaly_fX,
    light_meson_fX_for,
    precision_chain,
    verify_light_meson_fX,
    precision_chain_summary,
)


def test_cabibbo_scale_formula():
    """f = m · √(7α) directly."""
    m = 0.5
    f = cabibbo_scale_fX(m)
    assert f == pytest.approx(m * math.sqrt(7.0 * ALPHA_SUBSTRATE), rel=1e-14)


def test_fibonacci_anomaly_formula():
    """f = m / F^(1/4) with default F = 5."""
    m = 0.13957
    f = fibonacci_anomaly_fX(m)
    assert f == pytest.approx(m / math.pow(5.0, 0.25), rel=1e-14)


def test_f_K_matches_pdg_within_3pct():
    """f_K substrate ≈ 111.5 MeV vs PDG 110.0 → +1.4 %."""
    f_K_MeV = light_meson_fX_for("K") * 1e3
    rel_err = f_K_MeV / 110.0 - 1.0
    assert abs(rel_err) < 0.03


def test_f_eta_within_6pct():
    """f_η substrate ≈ 123.8 MeV vs PDG 130.0 → -4.7 %."""
    f_eta_MeV = light_meson_fX_for("eta") * 1e3
    rel_err = f_eta_MeV / 130.0 - 1.0
    assert abs(rel_err) < 0.06


def test_f_pi_within_2pct():
    """f_π substrate (Fibonacci) ≈ 93.4 MeV vs PDG 92.4 → +1.1 %."""
    f_pi_MeV = light_meson_fX_for("pi") * 1e3
    rel_err = f_pi_MeV / 92.4 - 1.0
    assert abs(rel_err) < 0.02


def test_pi_uses_fibonacci_regime():
    """π is in the Fibonacci regime, not Cabibbo."""
    assert LIGHT_PSEUDOSCALARS["pi"].regime == "fibonacci"


def test_K_eta_use_cabibbo_regime():
    assert LIGHT_PSEUDOSCALARS["K"].regime == "cabibbo"
    assert LIGHT_PSEUDOSCALARS["eta"].regime == "cabibbo"


def test_unknown_meson_raises():
    with pytest.raises(KeyError):
        light_meson_fX_for("nothing")


def test_precision_chain_covers_three_mesons():
    chain = precision_chain()
    chain.pop('pass', None)
    assert set(chain.keys()) == {"pi", "K", "eta"}


def test_verify_passes_default_5pct():
    """Default 5 % tolerance — η is the worst at -4.7 %, passes by ~0.3 %."""
    result = verify_light_meson_fX()
    assert result['pass'] is True


def test_verify_fails_unrealistic_tol():
    """1 % tolerance fails (K +1.4 %, η -4.7 %, π +1.1 %)."""
    result = verify_light_meson_fX(tol=1e-2)
    assert result['pass'] is False


def test_summary_mentions_key_substrings():
    out = precision_chain_summary()
    assert "Cabibbo" in out
    assert "Fibonacci" in out
    assert "f_K" in out
    assert "f_π" in out or "f_pi" in out
