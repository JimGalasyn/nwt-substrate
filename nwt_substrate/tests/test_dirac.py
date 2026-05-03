"""Tests for nwt_substrate.algebra.dirac."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import (
    make_lorentz_gammas,
    slash_p,
    make_gamma5,
    ETA_LORENTZ,
)


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


def test_lorentz_clifford(gammas):
    """{gamma^mu, gamma^nu} = 2 eta^{mu nu} I_8."""
    I8 = np.eye(8, dtype=complex)
    for mu in range(4):
        for nu in range(4):
            anti = gammas[mu] @ gammas[nu] + gammas[nu] @ gammas[mu]
            expected = 2.0 * ETA_LORENTZ[mu, nu] * I8
            np.testing.assert_array_almost_equal(anti, expected, decimal=12)


def test_slash_p_squared(gammas):
    """(slash{p})^2 = p^2 I_8 for arbitrary 4-momenta."""
    rng = np.random.default_rng(42)
    for _ in range(5):
        p = rng.standard_normal(4)
        sp = slash_p(p, gammas)
        p_sq = p[0] ** 2 - p[1] ** 2 - p[2] ** 2 - p[3] ** 2
        expected = p_sq * np.eye(8, dtype=complex)
        np.testing.assert_array_almost_equal(sp @ sp, expected, decimal=10)


def test_gamma5_anticommutes(gammas):
    """{gamma^5, gamma^mu} = 0 for mu = 0..3."""
    g5 = make_gamma5(gammas)
    for mu in range(4):
        anti = g5 @ gammas[mu] + gammas[mu] @ g5
        np.testing.assert_array_almost_equal(anti, 0, decimal=12)


def test_gamma5_squared_is_identity(gammas):
    """(gamma^5)^2 = I_8."""
    g5 = make_gamma5(gammas)
    np.testing.assert_array_almost_equal(g5 @ g5, np.eye(8, dtype=complex),
                                          decimal=12)


def test_gamma5_traceless(gammas):
    """tr(gamma^5) = 0."""
    g5 = make_gamma5(gammas)
    assert abs(np.trace(g5)) < 1e-12
