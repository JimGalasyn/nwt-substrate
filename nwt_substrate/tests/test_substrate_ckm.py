"""Tests for nwt_substrate.electroweak.substrate_ckm -- substrate Wolfenstein CKM.

P6b closure (Galasyn 2026-05-23): full CKM in 4 substrate integers + α + δ=π/3.

    λ² = 7α                       (Cabibbo from K_7)
    A² = 9/14                     (q_b / dim(G_2))
    |ρ̄ + iη̄|² = 21 α             (dim(so(7)) · α = |E(K_7)| · α)
    δ_CP^CKM = π/3                (Z_2 carrier conjecture)

Zero free Wolfenstein parameters.  Dominant CKM elements match PDG to
sub-percent precision; closed-form Jarlskog J = 875 α^(7/2).
"""

import math

import numpy as np
import pytest

from nwt_substrate.electroweak.substrate_ckm import (
    ALPHA_SUBSTRATE,
    N_VERTICES_K7,
    Q_B,
    DIM_G2,
    DIM_SO7,
    DELTA_CP_CKM,
    wolfenstein_lambda_sq,
    wolfenstein_lambda,
    wolfenstein_A_sq,
    wolfenstein_A,
    wolfenstein_apex_magnitude_sq,
    wolfenstein_apex_magnitude,
    wolfenstein_rho_bar,
    wolfenstein_eta_bar,
    V_us,
    V_cb,
    V_ub,
    V_td,
    jarlskog_ckm,
    jarlskog_coefficient,
    ckm_matrix,
    precision_chain,
    verify_substrate_ckm,
    precision_chain_summary,
)


# ============================================================
# Substrate integers
# ============================================================

def test_substrate_integers():
    """The four substrate integers (7, 9, 14, 21) are the canonical set."""
    assert N_VERTICES_K7 == 7
    assert Q_B == 9
    assert DIM_G2 == 14
    assert DIM_SO7 == 21


def test_so7_equals_k7_edge_count():
    """dim(so(7)) = 21 = |E(K_7)| = C(7,2) substrate identity."""
    assert DIM_SO7 == N_VERTICES_K7 * (N_VERTICES_K7 - 1) // 2


# ============================================================
# Wolfenstein λ (Cabibbo)
# ============================================================

def test_lambda_sq_formula():
    """λ² = 7α."""
    assert wolfenstein_lambda_sq() == pytest.approx(7.0 * ALPHA_SUBSTRATE,
                                                    rel=1e-14)


def test_lambda_matches_pdg_within_1pct():
    """λ = √(7α) ≈ 0.22601, PDG 0.2253 → 0.32%."""
    lam = wolfenstein_lambda()
    assert lam == pytest.approx(0.22601, abs=1e-4)
    rel_err = lam / 0.2253 - 1.0
    assert abs(rel_err) < 0.01


# ============================================================
# Wolfenstein A
# ============================================================

def test_A_sq_formula():
    """A² = q_b / dim(G_2) = 9/14."""
    assert wolfenstein_A_sq() == pytest.approx(9.0 / 14.0, rel=1e-14)


def test_A_matches_pdg_within_1pct():
    """A = 3/√14 ≈ 0.80178, PDG 0.802 → 0.027%."""
    A = wolfenstein_A()
    assert A == pytest.approx(3.0 / math.sqrt(14.0), rel=1e-14)
    rel_err = A / 0.802 - 1.0
    assert abs(rel_err) < 0.01


# ============================================================
# Wolfenstein |ρ̄+iη̄|
# ============================================================

def test_apex_magnitude_formula():
    """|ρ̄ + iη̄|² = 21α."""
    assert wolfenstein_apex_magnitude_sq() == pytest.approx(
        21.0 * ALPHA_SUBSTRATE, rel=1e-14)


def test_apex_magnitude_matches_pdg_within_1pct():
    """|ρ̄+iη̄| = √(21α) ≈ 0.391, PDG 0.391 → 0.12%."""
    apex = wolfenstein_apex_magnitude()
    rel_err = apex / 0.391 - 1.0
    assert abs(rel_err) < 0.01


def test_rho_eta_with_delta_pi_over_3():
    """δ=π/3 gives ρ̄≈0.196 (cos π/3 = 1/2), η̄≈0.339 (sin π/3 = √3/2)."""
    apex = wolfenstein_apex_magnitude()
    rho = wolfenstein_rho_bar()
    eta = wolfenstein_eta_bar()
    assert rho == pytest.approx(apex * 0.5, rel=1e-14)
    assert eta == pytest.approx(apex * math.sqrt(3.0) / 2.0, rel=1e-14)


def test_delta_cp_ckm_is_pi_over_3():
    """δ_CP^CKM = π/3 substrate conjecture."""
    assert DELTA_CP_CKM == pytest.approx(math.pi / 3.0, rel=1e-14)


# ============================================================
# Individual CKM elements
# ============================================================

def test_V_us_equals_lambda():
    """|V_us| = λ to leading Wolfenstein."""
    assert V_us() == pytest.approx(wolfenstein_lambda(), rel=1e-14)


def test_V_cb_equals_A_lambda_sq():
    """|V_cb| = A λ²."""
    assert V_cb() == pytest.approx(wolfenstein_A() * wolfenstein_lambda_sq(),
                                   rel=1e-14)


def test_V_cb_matches_pdg_within_1pct():
    """|V_cb| substrate 0.04096 vs PDG 0.0410 → 0.1%."""
    rel_err = V_cb() / 0.0410 - 1.0
    assert abs(rel_err) < 0.01


def test_V_ub_magnitude_substrate_formula():
    """|V_ub| = A λ³ |ρ̄+iη̄|."""
    expected_mag = (wolfenstein_A() * wolfenstein_lambda() ** 3 *
                    wolfenstein_apex_magnitude())
    assert abs(V_ub()) == pytest.approx(expected_mag, rel=1e-14)


def test_V_ub_carries_phase():
    """V_ub ∝ (ρ̄ - iη̄), i.e. negative imaginary part."""
    assert V_ub().imag < 0.0


def test_V_td_substrate_formula():
    """V_td = A λ³ (1 - ρ̄ - iη̄)."""
    expected = (wolfenstein_A() * wolfenstein_lambda() ** 3 *
                complex(1.0 - wolfenstein_rho_bar(), -wolfenstein_eta_bar()))
    assert V_td() == pytest.approx(expected, rel=1e-14)


# ============================================================
# Closed-form Jarlskog
# ============================================================

def test_jarlskog_coefficient_value():
    """Coefficient 875.0 = 9·343·3√7 / 28 in J = (coeff) · α^(7/2)."""
    expected = (9.0 * 343.0 * 3.0 * math.sqrt(7.0)) / 28.0
    assert jarlskog_coefficient() == pytest.approx(expected, rel=1e-14)
    # And specifically close to 875.0
    assert 870.0 < jarlskog_coefficient() < 880.0


def test_jarlskog_equals_coefficient_alpha_7over2():
    """J_CP^CKM = 875 · α^(7/2) substrate closed form."""
    alpha = ALPHA_SUBSTRATE
    closed = jarlskog_coefficient() * alpha ** 3.5
    direct = jarlskog_ckm(alpha)
    assert direct == pytest.approx(closed, rel=1e-12)


def test_jarlskog_matches_pdg_within_15pct():
    """J_CP^CKM substrate 2.91e-5 vs PDG 3.18e-5 → 8.6%."""
    J = jarlskog_ckm()
    rel_err = J / 3.18e-5 - 1.0
    assert abs(rel_err) < 0.15


# ============================================================
# Full CKM matrix
# ============================================================

def test_ckm_matrix_shape_and_dtype():
    V = ckm_matrix()
    assert V.shape == (3, 3)
    assert V.dtype == complex


def test_ckm_matrix_top_row_us_element():
    """V[0,1] = |V_us| = λ (real)."""
    V = ckm_matrix()
    assert V[0, 1].real == pytest.approx(wolfenstein_lambda(), rel=1e-14)
    assert V[0, 1].imag == pytest.approx(0.0, abs=1e-14)


def test_ckm_matrix_V_ub_phase():
    """V[0,2] = V_ub carries the substrate δ_CP^CKM phase."""
    V = ckm_matrix()
    expected = V_ub()
    assert V[0, 2] == pytest.approx(expected, rel=1e-12)


def test_ckm_matrix_approx_unitarity():
    """At Wolfenstein O(λ⁴), the matrix is unitary to O(λ⁶) ~ 5e-4."""
    V = ckm_matrix()
    VVdag = V @ V.conj().T
    # Diagonal entries close to 1
    for i in range(3):
        assert abs(VVdag[i, i] - 1.0) < 1e-3
    # Off-diagonal entries small
    for i in range(3):
        for j in range(3):
            if i != j:
                assert abs(VVdag[i, j]) < 1e-3


# ============================================================
# Precision chain
# ============================================================

def test_precision_chain_keys():
    """precision_chain returns all expected observables."""
    chain = precision_chain()
    expected = {'lambda', 'A', 'apex_magnitude', 'V_us', 'V_cb',
                'V_ub_mag', 'V_td_mag', 'jarlskog', 'delta_deg'}
    assert set(chain.keys()) == expected
    for k in chain:
        assert set(chain[k].keys()) == {'substrate', 'pdg', 'rel_err'}


def test_precision_chain_dominant_within_1pct():
    """λ, A, |ρ̄+iη̄|, V_us, V_cb all match PDG within 1%."""
    chain = precision_chain()
    for k in ['lambda', 'A', 'apex_magnitude', 'V_us', 'V_cb']:
        assert abs(chain[k]['rel_err']) < 0.01, (
            f"{k} rel err {chain[k]['rel_err']*100:.2f}% exceeds 1%"
        )


def test_precision_chain_delta_tension():
    """δ_CP^CKM substrate 60° vs PDG 66.4° → ~10% tension (known open piece)."""
    chain = precision_chain()
    # Substrate is below PDG by ~10%
    assert chain['delta_deg']['rel_err'] < 0.0
    assert abs(chain['delta_deg']['rel_err']) < 0.15


# ============================================================
# Verification harness
# ============================================================

def test_verify_substrate_ckm_passes_default_tols():
    """Default 1% dominant / 10% subdominant tolerances pass."""
    result = verify_substrate_ckm()
    assert result['pass'] is True


def test_verify_substrate_ckm_fails_unrealistic_dominant_tol():
    """1 ppm dominant tolerance fails (substrate λ is 0.32% off)."""
    result = verify_substrate_ckm(tol_dominant=1e-6)
    assert result['pass'] is False


# ============================================================
# Pretty-print
# ============================================================

def test_precision_chain_summary_contains_key_substrings():
    """Summary mentions Wolfenstein params, CKM elements, Jarlskog."""
    out = precision_chain_summary()
    assert "λ" in out
    assert "Wolfenstein" in out
    assert "J_CP" in out
    assert "α^(7/2)" in out
    assert "K_7" in out


# ============================================================
# Integration via electroweak namespace
# ============================================================

def test_substrate_namespace_exposes_ckm():
    """ew.substrate.ckm() returns the substrate CKM matrix."""
    import nwt_substrate.electroweak as ew
    V = ew.substrate.ckm()
    assert V.shape == (3, 3)


def test_substrate_namespace_jarlskog():
    """ew.substrate.jarlskog_ckm() returns the substrate Jarlskog."""
    import nwt_substrate.electroweak as ew
    assert ew.substrate.jarlskog_ckm() == pytest.approx(jarlskog_ckm(),
                                                        rel=1e-14)
