"""Tests for nwt_substrate.algebra.octonions."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import (
    make_octonion_table,
    octo_mul,
    octo_assoc,
    octo_norm_sq,
    basis_vector,
    fano_partner,
    FANO_LINES,
)


@pytest.fixture(scope="module")
def T():
    return make_octonion_table()


def test_table_shape(T):
    assert T.shape == (8, 8, 8)


def test_identity_e0(T):
    """e_0 * e_i = e_i and e_i * e_0 = e_i for all i."""
    for i in range(8):
        e_i = basis_vector(i)
        e_0 = basis_vector(0)
        np.testing.assert_array_equal(octo_mul(e_0, e_i, T), e_i)
        np.testing.assert_array_equal(octo_mul(e_i, e_0, T), e_i)


def test_imaginary_squares(T):
    """e_i * e_i = -e_0 for i = 1..7."""
    e_0 = basis_vector(0)
    for i in range(1, 8):
        e_i = basis_vector(i)
        np.testing.assert_array_equal(octo_mul(e_i, e_i, T), -e_0)


def test_hurwitz_norm(T):
    """|a*b|^2 = |a|^2 * |b|^2."""
    rng = np.random.default_rng(42)
    for _ in range(10):
        a = rng.standard_normal(8)
        b = rng.standard_normal(8)
        ab = octo_mul(a, b, T)
        assert abs(octo_norm_sq(ab) - octo_norm_sq(a) * octo_norm_sq(b)) < 1e-10


def test_associator_antisymmetric(T):
    """[a,b,c] = -[b,a,c] = -[a,c,b]."""
    e1, e2, e4 = basis_vector(1), basis_vector(2), basis_vector(4)
    A_124 = octo_assoc(e1, e2, e4, T)
    A_214 = octo_assoc(e2, e1, e4, T)
    A_142 = octo_assoc(e1, e4, e2, T)
    np.testing.assert_array_almost_equal(A_124, -A_214, decimal=12)
    np.testing.assert_array_almost_equal(A_124, -A_142, decimal=12)


def test_fano_lines_associate(T):
    """Associator vanishes on Fano lines."""
    for (i, j, k) in FANO_LINES:
        e_i = basis_vector(i)
        e_j = basis_vector(j)
        e_k = basis_vector(k)
        A = octo_assoc(e_i, e_j, e_k, T)
        assert octo_norm_sq(A) < 1e-20


def test_associator_sum_112(T):
    """Sum |[e_i, e_j, e_k]|^2 over 28 non-Fano triples = 112 (Paper 17)."""
    fano_sets = [frozenset(line) for line in FANO_LINES]
    total = 0.0
    for i in range(1, 8):
        for j in range(i + 1, 8):
            for k in range(j + 1, 8):
                if frozenset((i, j, k)) in fano_sets:
                    continue
                A = octo_assoc(basis_vector(i), basis_vector(j),
                               basis_vector(k), T)
                total += octo_norm_sq(A)
    assert abs(total - 112.0) < 1e-10


def test_moufang_identity(T):
    """(zxz)y = z(x(zy))."""
    rng = np.random.default_rng(123)
    for _ in range(5):
        x = rng.standard_normal(8)
        y = rng.standard_normal(8)
        z = rng.standard_normal(8)
        lhs = octo_mul(octo_mul(octo_mul(z, x, T), z, T), y, T)
        rhs = octo_mul(z, octo_mul(x, octo_mul(z, y, T), T), T)
        np.testing.assert_array_almost_equal(lhs, rhs, decimal=10)


def test_fano_partner_signs():
    """e_1 * e_2 = +e_3 (line (1,2,3) cyclic)."""
    k, sign = fano_partner(1, 2)
    assert k == 3 and sign == +1
    k, sign = fano_partner(2, 1)
    assert k == 3 and sign == -1
