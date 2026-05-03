"""Tests for the 3-gluon vertex Lorentz structure and Bose symmetry."""

import numpy as np
import pytest

from nwt_substrate.amplitudes.vertices import (
    three_gluon_vertex_lorentz,
    three_gluon_vertex_full,
    four_gluon_vertex_lorentz,
)
from nwt_substrate.algebra.su3 import structure_constants


def test_three_gluon_lorentz_shape():
    """V_Lorentz returns a (4, 4, 4) tensor."""
    k1 = np.array([1.0, 0.5, 0.0, 0.0])
    k2 = np.array([1.0, -0.3, 0.4, 0.0])
    k3 = -(k1 + k2)
    V = three_gluon_vertex_lorentz(k1, k2, k3)
    assert V.shape == (4, 4, 4)


def test_three_gluon_bose_symmetry_under_12_swap():
    """V(k2,k1,k3)[nu,mu,rho] = -V(k1,k2,k3)[mu,nu,rho]."""
    k1 = np.array([1.0, 0.5, 0.3, 0.1])
    k2 = np.array([2.0, -0.2, 0.5, 0.4])
    k3 = -(k1 + k2)

    V = three_gluon_vertex_lorentz(k1, k2, k3)
    V_swap = three_gluon_vertex_lorentz(k2, k1, k3)
    V_swap_relabeled = np.transpose(V_swap, (1, 0, 2))   # swap mu <-> nu

    np.testing.assert_array_almost_equal(
        V_swap_relabeled, -V, decimal=12,
    )


def test_three_gluon_bose_symmetry_under_23_swap():
    """V(k1,k3,k2)[mu,rho,nu] = -V(k1,k2,k3)[mu,nu,rho]."""
    k1 = np.array([1.0, 0.5, 0.3, 0.1])
    k2 = np.array([2.0, -0.2, 0.5, 0.4])
    k3 = -(k1 + k2)

    V = three_gluon_vertex_lorentz(k1, k2, k3)
    V_swap = three_gluon_vertex_lorentz(k1, k3, k2)
    V_swap_relabeled = np.transpose(V_swap, (0, 2, 1))   # swap nu <-> rho

    np.testing.assert_array_almost_equal(
        V_swap_relabeled, -V, decimal=12,
    )


def test_three_gluon_full_with_color():
    """V_full = g_s f^abc V_Lorentz; f^012 = 1 (cyclic 1,2,3)."""
    k1 = np.array([1.0, 0.0, 0.0, 1.0])
    k2 = np.array([1.0, 0.0, 0.0, -1.0])
    k3 = -(k1 + k2)

    V_lor = three_gluon_vertex_lorentz(k1, k2, k3)
    V_full = three_gluon_vertex_full(k1, k2, k3, a=0, b=1, c=2, gs=1.5)
    f = structure_constants()
    expected = 1.5 * f[0, 1, 2] * V_lor

    np.testing.assert_array_almost_equal(V_full, expected, decimal=12)


def test_three_gluon_color_factor_zero_for_diagonal():
    """f^aab = 0 (antisymmetric) -> V_full = 0 for a=b."""
    k1 = np.array([1.0, 0.0, 0.0, 1.0])
    k2 = np.array([1.0, 0.0, 0.0, -1.0])
    k3 = -(k1 + k2)

    V = three_gluon_vertex_full(k1, k2, k3, a=2, b=2, c=5, gs=1.0)
    assert np.all(np.abs(V) < 1e-15)


def test_four_gluon_lorentz_shape():
    """4-gluon returns 3 Lorentz tensors of shape (4, 4, 4, 4)."""
    L = four_gluon_vertex_lorentz()
    assert L.shape == (3, 4, 4, 4, 4)


def test_four_gluon_lorentz_antisymmetric_in_pair():
    """L[0] is antisymmetric under (mu,rho) <-> (nu,sigma) swap... wait, this
    is the mu<->nu (or rho<->sigma) antisymmetry. L[0] = eta^mu_rho eta^nu_sig
    - eta^mu_sig eta^nu_rho is antisymmetric under (rho <-> sigma)."""
    L = four_gluon_vertex_lorentz()
    # L[0]: antisymmetric under rho <-> sigma
    L_swap = np.transpose(L[0], (0, 1, 3, 2))   # rho <-> sigma
    np.testing.assert_array_almost_equal(L[0], -L_swap, decimal=12)
