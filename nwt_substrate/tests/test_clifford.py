"""Tests for nwt_substrate.algebra.clifford."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table, basis_vector
from nwt_substrate.algebra.clifford import (
    left_mult_matrix,
    bivector_matrix,
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
