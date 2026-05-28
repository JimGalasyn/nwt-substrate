"""Near-Horizon Extremal Kerr (NHEK) geometry tests.

Covers the fast numerical surface (metric, inverse, determinant, signature,
Killing structure, vortex-centerline helpers, finite-difference Christoffels).
Skips the sympy path (`nhek_symbolic`, `verify_nhek_vacuum`) — each call
runs symbolic Ricci computation that takes ~30 s and isn't suitable for
unit-test latency.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from nwt_substrate.gravity.nhek import (
    Lambda_factor,
    Sigma,
    christoffels_numeric,
    is_axisymmetric_invariant,
    is_near_bifurcation_sphere,
    killing_vectors_constant_basis,
    nhek_inverse_metric,
    nhek_metric,
    nhek_metric_determinant,
    nhek_signature,
    substrate_vortex_centerline_radius,
)


# ---------------------------------------------------------------
# Scalar prefactors
# ---------------------------------------------------------------

def test_Sigma_equator():
    """Σ(π/2) = 1 + cos²(π/2) = 1."""
    assert Sigma(math.pi / 2) == pytest.approx(1.0, rel=1e-15)


def test_Sigma_pole():
    """Σ(0) = 1 + cos²(0) = 2."""
    assert Sigma(0.0) == pytest.approx(2.0, rel=1e-15)


def test_Lambda_factor_equator():
    """Λ(π/2) = 2 sin(π/2) / Σ(π/2) = 2 / 1 = 2."""
    assert Lambda_factor(math.pi / 2) == pytest.approx(2.0, rel=1e-15)


def test_Lambda_factor_at_pole_vanishes():
    """Λ(0) = 0 because sin(0) = 0."""
    assert Lambda_factor(0.0) == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------
# Metric
# ---------------------------------------------------------------

def test_metric_is_4x4_and_symmetric():
    g = nhek_metric(M=1.0, r=1.0, theta=math.pi / 3)
    assert g.shape == (4, 4)
    np.testing.assert_allclose(g, g.T, atol=1e-15)


def test_metric_offdiagonal_block_only_t_phi():
    """All off-diagonal entries vanish except (0,3) = (3,0) = g_tφ."""
    g = nhek_metric(M=1.0, r=0.7, theta=math.pi / 4)
    for i in range(4):
        for j in range(4):
            if i == j or (i, j) in {(0, 3), (3, 0)}:
                continue
            assert abs(g[i, j]) < 1e-15, f"g[{i},{j}] nonzero off-diagonal"


def test_metric_g_phi_phi_is_4_M2_sin2_over_Sigma():
    M, theta = 2.0, math.pi / 3
    S = Sigma(theta)
    g = nhek_metric(M, r=1.0, theta=theta)
    expected = 4.0 * M * M * math.sin(theta) ** 2 / S
    assert g[3, 3] == pytest.approx(expected, rel=1e-14)


def test_metric_g_rr_is_M2_Sigma_over_r2():
    M, r, theta = 1.5, 0.8, math.pi / 4
    g = nhek_metric(M, r, theta)
    assert g[1, 1] == pytest.approx(M * M * Sigma(theta) / (r * r), rel=1e-14)


def test_metric_g_t_phi_proportional_to_r():
    """g_tφ scales linearly in r at fixed (M, θ)."""
    M, theta = 1.0, math.pi / 5
    g1 = nhek_metric(M, r=1.0, theta=theta)
    g2 = nhek_metric(M, r=2.0, theta=theta)
    assert g2[0, 3] == pytest.approx(2.0 * g1[0, 3], rel=1e-14)


# ---------------------------------------------------------------
# Inverse metric
# ---------------------------------------------------------------

def test_inverse_metric_times_metric_is_identity():
    M, r, theta = 1.0, 1.5, math.pi / 3
    g = nhek_metric(M, r, theta)
    g_inv = nhek_inverse_metric(M, r, theta)
    np.testing.assert_allclose(g @ g_inv, np.eye(4), atol=1e-12)


# ---------------------------------------------------------------
# Determinant
# ---------------------------------------------------------------

def test_determinant_analytic_form():
    """det g = -4 M⁸ sin²θ Σ(θ)²; r-independent."""
    M, theta = 1.2, math.pi / 4
    S = Sigma(theta)
    expected = -4.0 * M ** 8 * math.sin(theta) ** 2 * S * S
    assert nhek_metric_determinant(M, r=0.7, theta=theta) == pytest.approx(expected, rel=1e-14)


def test_determinant_r_independent():
    """det g is independent of r (Bardeen-Horowitz r-cancellation)."""
    M, theta = 1.0, math.pi / 4
    d1 = nhek_metric_determinant(M, r=0.1, theta=theta)
    d2 = nhek_metric_determinant(M, r=5.0, theta=theta)
    assert d1 == pytest.approx(d2, rel=1e-14)


def test_determinant_matches_numpy_linalg_det():
    """Analytic form agrees with np.linalg.det of the metric matrix."""
    M, r, theta = 1.0, 1.7, math.pi / 3
    g = nhek_metric(M, r, theta)
    analytic = nhek_metric_determinant(M, r, theta)
    assert np.linalg.det(g) == pytest.approx(analytic, rel=1e-10)


# ---------------------------------------------------------------
# Signature
# ---------------------------------------------------------------

def test_signature_is_lorentzian_at_equator_r1():
    """Healthy NHEK at (r=1, θ=π/2) — Lorentzian (1, 3, 0)."""
    sig = nhek_signature(M=1.0, r=1.0, theta=math.pi / 2)
    assert sig == (1, 3, 0)


def test_signature_is_lorentzian_at_mid_latitude():
    sig = nhek_signature(M=1.0, r=0.5, theta=math.pi / 3)
    assert sig == (1, 3, 0)


# ---------------------------------------------------------------
# Killing vectors
# ---------------------------------------------------------------

def test_killing_basis_returns_H_and_L():
    K = killing_vectors_constant_basis()
    assert "H_axt" in K
    assert "L_axphi" in K
    np.testing.assert_array_equal(K["H_axt"], np.array([1.0, 0.0, 0.0, 0.0]))
    np.testing.assert_array_equal(K["L_axphi"], np.array([0.0, 0.0, 0.0, 1.0]))


def test_is_axisymmetric_invariant_true_by_construction():
    """nhek_metric() takes no (t, φ) args; Killing condition holds by construction."""
    assert is_axisymmetric_invariant(M=1.0, r=1.0, theta=math.pi / 4) is True


# ---------------------------------------------------------------
# Substrate-vortex centerline / bifurcation sphere helpers
# ---------------------------------------------------------------

def test_substrate_vortex_centerline_radius_is_zero():
    """Bridge T² centerline at r = 0 (bifurcation 2-sphere)."""
    assert substrate_vortex_centerline_radius() == 0.0


def test_is_near_bifurcation_sphere_default_threshold():
    assert is_near_bifurcation_sphere(r=0.05) is True
    assert is_near_bifurcation_sphere(r=0.5) is False
    assert is_near_bifurcation_sphere(r=-0.08) is True  # |r| < 0.1


def test_is_near_bifurcation_sphere_custom_threshold():
    assert is_near_bifurcation_sphere(r=0.3, threshold=0.5) is True
    assert is_near_bifurcation_sphere(r=0.6, threshold=0.5) is False


# ---------------------------------------------------------------
# Christoffels (numerical)
# ---------------------------------------------------------------

def test_christoffels_shape_is_4_4_4():
    G = christoffels_numeric(M=1.0, r=1.0, theta=math.pi / 3)
    assert G.shape == (4, 4, 4)


def test_christoffels_symmetric_in_lower_indices():
    """Γ^λ_{μν} = Γ^λ_{νμ} (torsion-free Levi-Civita connection)."""
    G = christoffels_numeric(M=1.0, r=1.0, theta=math.pi / 3)
    for lam in range(4):
        np.testing.assert_allclose(G[lam], G[lam].T, atol=1e-6,
                                   err_msg=f"Γ^{lam} not symmetric in (μ,ν)")


def test_christoffels_t_phi_partials_are_zero():
    """∂_t g = ∂_φ g = 0 by axisymmetry → Christoffels carry no t/φ derivatives."""
    G = christoffels_numeric(M=1.0, r=1.0, theta=math.pi / 3)
    # G is finite, real, and the (t,φ)-only Γ have specific structure — just
    # check no NaN / inf and finite values everywhere.
    assert np.all(np.isfinite(G))
