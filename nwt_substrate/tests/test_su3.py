"""Tests for SU(3) color algebra."""

import numpy as np
import pytest

from nwt_substrate.algebra.su3 import (
    gell_mann_matrices, su3_generators,
    structure_constants, d_constants,
    fundamental_casimir, adjoint_casimir,
    C_F, C_A, N_C,
)


def test_gell_mann_count():
    """SU(3) has 8 generators."""
    assert len(gell_mann_matrices()) == 8


def test_gell_mann_hermitian():
    """All Gell-Mann matrices are Hermitian."""
    for L in gell_mann_matrices():
        np.testing.assert_array_almost_equal(L, L.conj().T, decimal=12)


def test_gell_mann_traceless():
    """Tr(lambda^a) = 0 for all a."""
    for L in gell_mann_matrices():
        assert abs(np.trace(L)) < 1e-12


def test_generators_orthonormal():
    """Tr(T^a T^b) = (1/2) delta^ab."""
    T = su3_generators()
    for a in range(8):
        for b in range(8):
            tr = np.trace(T[a] @ T[b]).real
            expected = 0.5 if a == b else 0.0
            assert abs(tr - expected) < 1e-12


def test_lie_algebra():
    """[T^a, T^b] = i f^abc T^c."""
    T = su3_generators()
    f = structure_constants()
    for a in range(8):
        for b in range(a + 1, 8):
            comm = T[a] @ T[b] - T[b] @ T[a]
            rhs = sum(1j * f[a, b, c] * T[c] for c in range(8))
            np.testing.assert_array_almost_equal(comm, rhs, decimal=12)


def test_fundamental_casimir():
    """T^a T^a = C_F * I = (4/3) * I_3."""
    casimir = fundamental_casimir()
    expected = C_F * np.eye(3, dtype=complex)
    np.testing.assert_array_almost_equal(casimir, expected, decimal=12)


def test_adjoint_casimir():
    """f^acd f^bcd = C_A * delta^ab = 3 * delta^ab."""
    adj = adjoint_casimir()
    expected = C_A * np.eye(8)
    np.testing.assert_array_almost_equal(adj, expected, decimal=12)


def test_structure_constants_antisymmetric():
    """f^abc is totally antisymmetric."""
    f = structure_constants()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                # f^abc = -f^bac = -f^acb
                assert abs(f[a, b, c] + f[b, a, c]) < 1e-12
                assert abs(f[a, b, c] + f[a, c, b]) < 1e-12


def test_d_constants_symmetric():
    """d^abc is totally symmetric."""
    d = d_constants()
    for a in range(8):
        for b in range(8):
            for c in range(8):
                assert abs(d[a, b, c] - d[b, a, c]) < 1e-12
                assert abs(d[a, b, c] - d[a, c, b]) < 1e-12


def test_known_structure_constants():
    """f^123 = 1, f^458 = f^678 = sqrt(3)/2."""
    f = structure_constants()
    assert abs(f[0, 1, 2] - 1.0) < 1e-12
    assert abs(f[3, 4, 7] - np.sqrt(3.0) / 2.0) < 1e-12
    assert abs(f[5, 6, 7] - np.sqrt(3.0) / 2.0) < 1e-12
    assert abs(f[0, 3, 6] - 0.5) < 1e-12   # f^147


def test_casimir_constants():
    """C_F = 4/3, C_A = N_c = 3."""
    assert C_F == 4.0 / 3.0
    assert C_A == 3.0
    assert N_C == 3
