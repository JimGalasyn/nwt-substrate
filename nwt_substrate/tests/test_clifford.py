"""Tests for nwt_substrate.algebra.clifford."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table, basis_vector
from nwt_substrate.algebra.clifford import (
    left_mult_matrix,
    bivector_matrix,
    gamma_is_antisymmetric,
    kform_matrix,
    kform_symmetry,
    pauli_decomposition_dims,
    verify_substrate_pauli,
)


@pytest.fixture(scope="module")
def T():
    return make_octonion_table()


def test_L0_is_identity(T):
    """L_0 = I_8."""
    L0 = left_mult_matrix(0, T)
    np.testing.assert_array_equal(L0, np.eye(8))


def test_imaginary_squares_to_minus_I(T):
    """L_i^2 = -I_8 for i = 1..7."""
    for i in range(1, 8):
        L_i = left_mult_matrix(i, T)
        np.testing.assert_array_almost_equal(L_i @ L_i, -np.eye(8), decimal=12)


def test_clifford_anticommutation(T):
    """{L_i, L_j} = -2 delta_ij I_8 for i, j in 1..7."""
    Ls = [left_mult_matrix(i, T) for i in range(8)]
    I8 = np.eye(8)
    for i in range(1, 8):
        for j in range(1, 8):
            anti = Ls[i] @ Ls[j] + Ls[j] @ Ls[i]
            expected = -2.0 * (1 if i == j else 0) * I8
            np.testing.assert_array_almost_equal(anti, expected, decimal=12)


def test_bivector_is_single_L(T):
    """L_(e_a * e_b) = +/- L_c where c = Fano partner of (a, b)."""
    from nwt_substrate.algebra.octonions import fano_partner
    for a in range(1, 8):
        for b in range(1, 8):
            if a == b:
                continue
            biv = bivector_matrix(a, b, T)
            c, sign = fano_partner(a, b)
            L_c = left_mult_matrix(c, T)
            np.testing.assert_array_almost_equal(biv, sign * L_c, decimal=12)


# ============================================================
# Substrate Pauli theorem (Galasyn 2026-05-23 P9)
# ============================================================
#
# Tests that the antisymmetric two-fermion subspace L^2(8_spinor) is
# the substrate K_8 edge space (dim 28 = 7 + 21).


def test_imaginary_gammas_are_antisymmetric(T):
    """(L_i)^T = -L_i for all i in 1..7 (octonion construction property)."""
    for i in range(1, 8):
        assert gamma_is_antisymmetric(i, T), f"L_{i} should be antisymmetric"


def test_gamma_zero_is_not_antisymmetric(T):
    """L_0 = I_8 is symmetric, not antisymmetric."""
    assert not gamma_is_antisymmetric(0, T)


def test_kform_matrix_rank_0(T):
    """Empty product is the 8x8 identity."""
    M = kform_matrix((), T)
    np.testing.assert_array_equal(M, np.eye(8))


def test_kform_matrix_rank_1(T):
    """Rank-1 k-form is just L_i."""
    for i in range(1, 8):
        M = kform_matrix((i,), T)
        np.testing.assert_array_almost_equal(M, left_mult_matrix(i, T), decimal=12)


@pytest.mark.parametrize("k, expected", [
    (0, 'SYM'),       # scalar
    (1, 'ANTISYM'),   # vector (7)
    (2, 'ANTISYM'),   # 2-form / adjoint (21)
    (3, 'SYM'),       # 3-form (35)
    (4, 'SYM'),       # 4-form (35, Hodge dual to 3-form)
    (5, 'ANTISYM'),   # 5-form (21, Hodge dual to 2-form)
    (6, 'ANTISYM'),   # 6-form (7, Hodge dual to vector)
    (7, 'SYM'),       # volume (1, Hodge dual to scalar)
])
def test_bilinear_symmetry_table(T, k, expected):
    """Cl(0,7) bilinear-symmetry table -- verifies the substrate-Pauli table.

    For Spin(7) Majorana spinors with C = I_8, the bilinear
    bar(psi) gamma^{mu_1...mu_k} chi is exchange-SYM iff the matrix
    gamma^{mu_1...mu_k} is matrix-SYM.  The k-vs-symmetry pattern is
    determined by the octonion-Clifford structure.
    """
    sym_type = kform_symmetry(k, T)
    assert sym_type == expected, (
        f"Rank-{k} bilinear symmetry: got {sym_type}, expected {expected}"
    )


def test_substrate_pauli_lambda_squared_dim(T):
    """Substrate Pauli theorem: dim L^2(8_spinor) = 28 = |E(K_8)|.

    Antisymmetric two-fermion subspace = vector (7) + adjoint (21) of Spin(7)
    = K_8 edge count = 28.
    """
    result = pauli_decomposition_dims(T)
    assert result['antisym'] == 28, (
        f"L^2(8) dim = {result['antisym']}, expected 28 = |E(K_8)|"
    )
    assert result['k8_edge_match'] is True


def test_substrate_pauli_sym_squared_dim(T):
    """Sym^2(8_spinor) = 1 (scalar) + 35 (4-form) = 36."""
    result = pauli_decomposition_dims(T)
    assert result['sym'] == 36, (
        f"Sym^2(8) dim = {result['sym']}, expected 36"
    )


def test_substrate_pauli_decomposition(T):
    """L^2(8) decomposition is exactly {1: vector, 2: adjoint}; Sym^2 is {0: scalar, 3: 4-form}."""
    result = pauli_decomposition_dims(T)
    assert set(result['antisym_decomp'].keys()) == {1, 2}
    assert set(result['sym_decomp'].keys()) == {0, 3}


def test_total_decomposition_sums_to_64(T):
    """L^2 + Sym^2 = 28 + 36 = 64 = 8 * 8."""
    result = pauli_decomposition_dims(T)
    assert result['antisym'] + result['sym'] == 64


def test_verify_substrate_pauli_standalone():
    """Substrate-Pauli verification with no T provided (builds fresh)."""
    result = verify_substrate_pauli()
    assert result['k8_edge_match'] is True
    assert result['antisym'] == 28
    assert result['sym'] == 36
    assert result['all_gammas_antisymmetric'] is True
