"""Tests for nwt_substrate.gravity.flrw_test.

Validates the NWT prediction module for the Clifton-Heinesen /
Koksbang-Heinesen FLRW null-test signal (arXiv:2604.07244, 2604.05836,
2604.05822).
"""

import math
import pytest

import nwt_substrate.gravity as grav
from nwt_substrate.gravity.flrw_test import (
    Q0_from_density_variance,
    scaling_index_n_NWT,
    R_avg_scaling,
    K_BR_scaling,
    A_of_z_scaling,
    O_of_z_scaling,
    nwt_prediction,
    amplitude_check,
    SIGMA_M_DEFAULT,
    H_0_DEFAULT,
)


# ---------------------------------------------------------------------------
# FLRW limit: Q_0 = 0 → all backreaction observables collapse to FLRW
# ---------------------------------------------------------------------------

def test_FLRW_limit_A_equals_minusK():
    """A = -K identically in FLRW (any Q_0 → 0)."""
    for K in (-1.0, 0.0, +1.0):
        assert A_of_z_scaling(0.5, Q_0=0.0, n=-1.0, K=K) == pytest.approx(-K)
        assert A_of_z_scaling(2.0, Q_0=0.0, n=-1.0, K=K) == pytest.approx(-K)


def test_FLRW_limit_O_equals_minusK():
    """O = -K identically in FLRW."""
    for K in (-1.0, 0.0, +1.0):
        assert O_of_z_scaling(1.0, Q_0=0.0, n=-1.0, K=K) == pytest.approx(-K)


def test_FLRW_limit_K_BR_zero_for_K0():
    """K_BR vanishes in FLRW with K=0 and Q_0=0."""
    assert K_BR_scaling(0.5, Q_0=0.0, n=-1.0, K=0.0) == pytest.approx(0.0)
    assert K_BR_scaling(1.5, Q_0=0.0, n=-1.0, K=0.0) == pytest.approx(0.0)


def test_FLRW_limit_R_avg_scales_as_a_squared():
    """In FLRW, <R> = 6K(1+z)² for Q_0=0."""
    K = 1.0
    for z in (0.0, 0.5, 1.0, 2.0):
        assert R_avg_scaling(z, Q_0=0.0, n=-1.0, K=K) == pytest.approx(
            6.0 * K * (1.0 + z)**2
        )


# ---------------------------------------------------------------------------
# Q_0 amplitude formula
# ---------------------------------------------------------------------------

def test_Q0_proportional_to_sigma_M_squared():
    """Q_0 ∝ σ_M² as required by the variance formula."""
    Q1 = Q0_from_density_variance(sigma_M=0.5)
    Q2 = Q0_from_density_variance(sigma_M=1.0)
    # Q2 should be 4× Q1
    assert Q2 / Q1 == pytest.approx(4.0)


def test_Q0_six_H0_squared_at_sigma_unity():
    """Q_0 = 6 H_0² when σ_M = 1."""
    Q = Q0_from_density_variance(sigma_M=1.0, H_0_kms_Mpc=H_0_DEFAULT)
    assert Q == pytest.approx(6.0 * H_0_DEFAULT**2, rel=1e-12)


def test_Q0_zero_when_sigma_M_zero():
    """No matter clustering → no backreaction."""
    assert Q0_from_density_variance(sigma_M=0.0) == 0.0


# ---------------------------------------------------------------------------
# Scaling index n
# ---------------------------------------------------------------------------

def test_scaling_index_in_favored_window():
    """NWT matter-era n = -1 lies in Clifton-Heinesen window -2 < n < 0."""
    n = scaling_index_n_NWT(matter_era=True)
    assert -2.0 < n < 0.0


# ---------------------------------------------------------------------------
# Sign of A and O for the favored regime
# ---------------------------------------------------------------------------

def test_A_positive_for_isotropic_substrate_K0():
    """Isotropic substrate (Q_0 > 0), n = -1, K = 0 → A > 0 (Eq. 25)."""
    Q_0 = Q0_from_density_variance(sigma_M=SIGMA_M_DEFAULT)
    for z in (0.5, 1.0, 1.5, 2.0):
        assert A_of_z_scaling(z, Q_0, n=-1.0, K=0.0) > 0.0


def test_A_negative_for_shear_dominated_substrate_K0():
    """Shear-dominated substrate (Q_0 < 0) → A < 0, matches observed sign."""
    Q_0 = -Q0_from_density_variance(sigma_M=SIGMA_M_DEFAULT)  # negative
    for z in (0.5, 1.0, 1.5, 2.0):
        assert A_of_z_scaling(z, Q_0, n=-1.0, K=0.0) < 0.0


# ---------------------------------------------------------------------------
# Eq. 25 algebraic identity
# ---------------------------------------------------------------------------

def test_A_eq25_explicit_formula():
    """Hand-check: A = -K - n/(n+2) Q_0 (1+z)^(-n-2) at z=0 with K=0."""
    Q_0 = 1000.0
    n = -1.0
    A = A_of_z_scaling(z=0.0, Q_0=Q_0, n=n, K=0.0)
    # -(-1/1) * 1000 / 1 = 1000
    assert A == pytest.approx(1000.0, rel=1e-12)


def test_K_BR_equals_negative_A_for_K0():
    """A = -K_BR (Clifton-Heinesen Eq. 21) for K = 0."""
    Q_0 = 5000.0
    for z in (0.0, 0.5, 1.0, 2.0):
        A = A_of_z_scaling(z, Q_0, n=-1.0, K=0.0)
        K_BR = K_BR_scaling(z, Q_0, n=-1.0, K=0.0)
        assert A == pytest.approx(-K_BR, rel=1e-10)


# ---------------------------------------------------------------------------
# n = -2 singularity guard
# ---------------------------------------------------------------------------

def test_n_minus_2_raises():
    """n = -2 is the FLRW-curvature special case (Eq. 24 footnote)."""
    with pytest.raises(ValueError):
        A_of_z_scaling(1.0, Q_0=1.0, n=-2.0, K=0.0)


# ---------------------------------------------------------------------------
# Amplitude lands in the observed range for plausible σ_M
# ---------------------------------------------------------------------------

def test_amplitude_check_consistency():
    """Amplitude check returns sensible numbers."""
    info = amplitude_check(sigma_M=SIGMA_M_DEFAULT)
    assert info["Q_0_kmsMpc2"] > 0
    assert info["Q_0_over_H0_squared"] == pytest.approx(6.0 * SIGMA_M_DEFAULT**2)
    # Required σ_M for observed band should be in (0.2, 0.6) -- plausible
    # for ~30-50 Mpc effective averaging scale.
    assert 0.2 < info["sigma_M_required_for_lower_obs"] < 0.6
    assert 0.2 < info["sigma_M_required_for_upper_obs"] < 0.6


def test_NWT_prediction_in_observed_band_for_intermediate_sigma():
    """For σ_M in [0.45, 0.70] (intermediate scale, ~30-50 Mpc effective
    averaging), |O| at z=1 lands in the observed [2500, 7500] band."""
    for sigma_M in (0.45, 0.5, 0.6, 0.7):
        Q_0 = Q0_from_density_variance(sigma_M=sigma_M)
        O_pred = abs(A_of_z_scaling(1.0, Q_0, n=-1.0, K=0.0))
        assert 2500.0 <= O_pred <= 7500.0, (
            f"sigma_M={sigma_M}: |O(z=1)|={O_pred:.0f} outside observed band"
        )


# ---------------------------------------------------------------------------
# Wrapper class
# ---------------------------------------------------------------------------

def test_nwt_prediction_isotropic_returns_positive_O():
    pred = nwt_prediction(z=1.0, sigma_M=SIGMA_M_DEFAULT, Q0_sign=+1.0)
    assert pred.O > 0
    assert pred.sign == "+"


def test_nwt_prediction_shear_returns_negative_O():
    pred = nwt_prediction(z=1.0, sigma_M=SIGMA_M_DEFAULT, Q0_sign=-1.0)
    assert pred.O < 0
    assert pred.sign == "-"


def test_nwt_prediction_R_avg_negative_for_isotropic_K0():
    """For -2 < n < 0 with positive Q, <R> is strongly negative
    (timescape-like: voids dominate)."""
    pred = nwt_prediction(z=1.0, sigma_M=SIGMA_M_DEFAULT, Q0_sign=+1.0, K=0.0)
    assert pred.R_avg < 0.0
