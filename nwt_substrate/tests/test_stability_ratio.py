"""Tests for nwt_substrate.particles.stability_ratio.

P7b synthesis closure: substrate pattern-stability ratio ρ_X = m_X / Γ_X
for K_7 walks (Galasyn 2026-05-23).
"""

import math

import pytest

from nwt_substrate.particles.stability_ratio import (
    COMPENDIUM_STABILITY,
    stability_ratio,
    stability_ratio_for,
    log10_stability_ratio_for,
    classify_regime,
    all_k7_walks_are_passive_or_BPS,
    stability_summary,
)


# ============================================================
# Formula
# ============================================================

def test_stability_ratio_basic():
    """ρ = m / Γ directly."""
    assert stability_ratio(1.0, 0.1) == pytest.approx(10.0, rel=1e-14)


def test_stable_particle_has_infinite_rho():
    """Γ = 0 → ρ = ∞."""
    assert math.isinf(stability_ratio(1.0, 0.0))


def test_negative_mass_rejected():
    with pytest.raises(ValueError):
        stability_ratio(-1.0, 0.1)


def test_negative_width_rejected():
    with pytest.raises(ValueError):
        stability_ratio(1.0, -0.1)


# ============================================================
# Compendium values
# ============================================================

def test_proton_is_stable():
    assert math.isinf(stability_ratio_for("p"))


def test_electron_is_stable():
    assert math.isinf(stability_ratio_for("e"))


def test_muon_rho_within_order_of_magnitude():
    """μ⁻: ρ ≈ 3.5×10¹⁷ per memo §4.1."""
    rho = stability_ratio_for("mu")
    assert 2e17 < rho < 5e17


def test_pi0_rho_within_order_of_magnitude():
    """π⁰: ρ ≈ 1.7×10⁷ per memo §4.1."""
    rho = stability_ratio_for("pi0")
    assert 1e7 < rho < 3e7


def test_upsilon_rho_within_order_of_magnitude():
    """Υ: ρ ≈ 1.8×10⁵ per memo §4.1."""
    rho = stability_ratio_for("Upsilon")
    assert 1e5 < rho < 3e5


def test_log10_rho_finite_for_unstable():
    """log_10(ρ) is finite for any unstable particle."""
    assert math.isfinite(log10_stability_ratio_for("mu"))


def test_log10_rho_infinite_for_stable():
    assert math.isinf(log10_stability_ratio_for("p"))


# ============================================================
# Regime classification
# ============================================================

def test_regime_BPS():
    assert classify_regime(math.inf) == "BPS"


def test_regime_passive():
    assert classify_regime(1e10) == "passive"


def test_regime_critical():
    assert classify_regime(2.5) == "critical"


def test_regime_decaying():
    assert classify_regime(0.1) == "decaying"


def test_all_compendium_particles_passive_or_BPS():
    """Substrate-monism reading: K_7 walks live entirely in the
    re-inscription-dominated (or BPS) regime; no decaying-regime
    compendium particle exists."""
    assert all_k7_walks_are_passive_or_BPS()


def test_no_compendium_particle_is_critical_or_decaying():
    """Strengthen: explicitly verify each one is passive or BPS."""
    for name in COMPENDIUM_STABILITY:
        regime = classify_regime(stability_ratio_for(name))
        assert regime in {"passive", "BPS"}, (
            f"{name} regime is {regime} (ρ = {stability_ratio_for(name)})"
        )


# ============================================================
# Pretty-print
# ============================================================

def test_summary_mentions_proton_and_muon():
    out = stability_summary()
    assert "p " in out or "p\t" in out or "p  " in out
    assert "mu" in out
    assert "BPS" in out or "passive" in out
