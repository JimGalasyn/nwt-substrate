"""Tests for nwt_substrate.amplitudes."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas, slash_p
from nwt_substrate.amplitudes.spinors import (
    positive_energy_spinors, negative_energy_spinors, adjoint_spinor,
)
from nwt_substrate.amplitudes.polarizations import (
    transverse_polarizations, massive_polarizations,
)
from nwt_substrate.amplitudes.propagators import (
    fermion_propagator, scalar_propagator, photon_propagator,
)
from nwt_substrate.amplitudes.processes import compton


M_E = 0.510998928e-3   # GeV


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


# ----- Spinors -----

def test_positive_energy_spinors_count(gammas):
    """slash p has 4 +m eigenstates."""
    p = np.array([1.0, 0.3, 0.4, 0.5])  # arbitrary
    m = np.sqrt(1.0**2 - 0.3**2 - 0.4**2 - 0.5**2)
    u = positive_energy_spinors(p, gammas, m)
    assert u.shape == (4, 8)


def test_positive_energy_eigenequation(gammas):
    """(slash p) u = m u for each row of u."""
    p = np.array([1.0, 0.0, 0.0, 0.5])
    m = np.sqrt(1.0**2 - 0.5**2)
    u = positive_energy_spinors(p, gammas, m)
    sp = slash_p(p, gammas)
    for u_s in u:
        np.testing.assert_array_almost_equal(sp @ u_s, m * u_s, decimal=8)


def test_negative_energy_spinors(gammas):
    """slash p has 4 -m eigenstates."""
    p = np.array([1.0, 0.0, 0.0, 0.5])
    m = np.sqrt(1.0**2 - 0.5**2)
    v = negative_energy_spinors(p, gammas, m)
    assert v.shape == (4, 8)
    sp = slash_p(p, gammas)
    for v_s in v:
        np.testing.assert_array_almost_equal(sp @ v_s, -m * v_s, decimal=8)


# ----- Polarizations -----

def test_transverse_polarizations_orthogonal_to_k():
    k = np.array([1.0, 0.0, 0.0, 1.0])  # massless along z
    eps1, eps2 = transverse_polarizations(k)
    # Lorenz: k.eps = k^0 eps^0 - k.spatial . eps.spatial
    for eps in (eps1, eps2):
        dot = k[0] * eps[0] - sum(k[i] * eps[i] for i in (1, 2, 3))
        assert abs(dot) < 1e-12


def test_massive_polarizations_three_modes():
    """Massive vector has 3 polarizations, all transverse to k."""
    k = np.array([2.0, 0.0, 0.0, np.sqrt(3.0)])  # k^2 = 4 - 3 = 1
    m = 1.0
    pols = massive_polarizations(k, m)
    assert len(pols) == 3
    for eps in pols:
        dot = k[0] * eps[0] - sum(k[i] * eps[i] for i in (1, 2, 3))
        assert abs(dot) < 1e-12


# ----- Propagators -----

def test_fermion_propagator_inverse(gammas):
    """S_F(p) (slash p - m) = 1 / (p^2 - m^2) * (slash p^2 - m^2) = I."""
    p = np.array([1.5, 0.3, 0.0, 0.0])
    m = 0.7
    S = fermion_propagator(p, gammas, m)
    sp = slash_p(p, gammas)
    I8 = np.eye(8, dtype=complex)
    # S_F = (slash p + m) / (p^2 - m^2)  ->  S_F (slash p - m) = (p^2 - m^2)/(p^2 - m^2) I = I
    product = S @ (sp - m * I8)
    np.testing.assert_array_almost_equal(product, I8, decimal=8)


def test_scalar_propagator():
    p = np.array([2.0, 0.5, 0.0, 0.0])
    m = 1.0
    val = scalar_propagator(p, m)
    p_sq = 2.0**2 - 0.5**2
    assert abs(val - 1.0 / (p_sq - m**2)) < 1e-12


def test_photon_propagator_shape():
    p = np.array([1.0, 0.5, 0.0, 0.0])
    D = photon_propagator(p)
    assert D.shape == (4, 4)


# ----- Compton (the key verification) -----

def test_compton_kinematics():
    """omega' = omega / (1 + (omega/m)(1 - cos theta))."""
    omega = 0.01
    m = M_E
    for theta in [0.5, 1.0, np.pi / 2, 2.5]:
        p, k, p_p, k_p, omega_p = compton.kinematics(omega, theta, m)
        expected = omega / (1.0 + (omega / m) * (1.0 - np.cos(theta)))
        assert abs(omega_p - expected) < 1e-15


def test_compton_matches_klein_nishina():
    """Substrate-algebra Compton |M|^2 reproduces Klein-Nishina to ~1e-12."""
    T = make_octonion_table()
    gammas = make_lorentz_gammas(T)
    omega = 0.01            # 10 MeV photon in lab
    for theta in [0.5, 1.0, 2.0]:
        M_sub, _ = compton.M_squared_avg(omega, theta, M_E, gammas)
        M_kn, _ = compton.klein_nishina_M_squared(omega, theta, M_E)
        # Substrate has 8-dim spinors -> averaging factor cancels;
        # the ratio of avg sums should equal exactly.
        # Memory states agreement to 4.8e-16 in walk-phase 3.
        ratio = M_sub / M_kn
        assert abs(ratio - 1.0) < 1e-9, (
            f"theta={theta}: ratio = {ratio}, |M|^2_sub = {M_sub}, "
            f"|M|^2_KN = {M_kn}"
        )
