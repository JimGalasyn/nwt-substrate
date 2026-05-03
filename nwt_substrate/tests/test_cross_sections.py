"""Tests for cross-section integration and standard observables."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas
from nwt_substrate.amplitudes.processes import compton, eemumu, moller, bhabha
from nwt_substrate.amplitudes.cross_sections import (
    cm_momentum,
    sigma_2to2_cm,
    sigma_compton_lab,
    thomson_cross_section,
    sigma_eemumu_high_energy,
    gamma_muon_decay_textbook,
)


M_E = 0.510998928e-3
M_MU = 0.10566


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


# ----- CM momentum helper -----

def test_cm_momentum_at_threshold():
    """At s = (m1 + m2)^2, p_CM = 0."""
    s = (M_E + M_MU) ** 2
    assert cm_momentum(s, M_E, M_MU) == pytest.approx(0.0, abs=1e-15)


def test_cm_momentum_high_energy():
    """At s >> max(m)^2, p_CM ~ sqrt(s)/2."""
    s = 1e6
    p = cm_momentum(s, M_E, M_E)
    assert p == pytest.approx(np.sqrt(s) / 2, rel=1e-12)


# ----- Compton -> Thomson -----

def test_compton_thomson_limit(gammas):
    """sigma(Compton) -> Thomson cross-section as omega -> 0."""
    sigma_thomson = thomson_cross_section(M_E)
    omega = 1e-7   # 0.1 keV — deep in Thomson regime

    def M_sq_at_theta(theta):
        return compton.M_squared_avg(omega, theta, M_E, gammas)

    sigma_sub = sigma_compton_lab(M_sq_at_theta, M_E)
    ratio = sigma_sub / sigma_thomson
    # Within 0.1% of Thomson value
    assert abs(ratio - 1.0) < 1e-3, f"ratio = {ratio}"


# ----- e+e- -> mu+mu- -----

def test_eemumu_total_cross_section(gammas):
    """sigma(e+e- -> mu+mu-) matches 4 pi alpha^2 / (3 s) exactly in massless limit."""
    for E_cm in [10.0, 30.0, 100.0]:
        s = E_cm ** 2

        def M_sq(theta):
            m, *_ = eemumu.M_squared_avg(E_cm, theta, M_E, M_MU, gammas)
            return m

        sigma_sub = sigma_2to2_cm(M_sq, s, (M_E, M_E), (M_MU, M_MU))
        sigma_text = sigma_eemumu_high_energy(s)
        assert abs(sigma_sub / sigma_text - 1.0) < 1e-7, (
            f"E_cm={E_cm}: ratio={sigma_sub/sigma_text}"
        )


# ----- Moller and Bhabha (with t-channel pole cutoff) -----

def test_moller_cross_section_matches_qed(gammas):
    """Integrated Moller |M|^2 (with cutoff) matches textbook integration."""
    E_cm = 10.0
    s = E_cm ** 2
    cutoff_deg = 30
    tmin = cutoff_deg * np.pi / 180.0
    tmax = np.pi - tmin

    def M_sq_sub(theta):
        m, *_ = moller.M_squared_avg(E_cm, theta, M_E, gammas)
        return m

    def M_sq_qed(theta):
        return moller.textbook_M_squared(E_cm, theta)

    sigma_sub = sigma_2to2_cm(M_sq_sub, s, (M_E, M_E), (M_E, M_E),
                              theta_min=tmin, theta_max=tmax)
    sigma_qed = sigma_2to2_cm(M_sq_qed, s, (M_E, M_E), (M_E, M_E),
                              theta_min=tmin, theta_max=tmax)
    assert abs(sigma_sub / sigma_qed - 1.0) < 1e-7


def test_bhabha_cross_section_matches_qed(gammas):
    """Integrated Bhabha |M|^2 (with cutoff) matches textbook integration."""
    E_cm = 10.0
    s = E_cm ** 2
    cutoff_deg = 30
    tmin = cutoff_deg * np.pi / 180.0
    tmax = np.pi - tmin

    def M_sq_sub(theta):
        m, *_ = bhabha.M_squared_avg(E_cm, theta, M_E, gammas)
        return m

    def M_sq_qed(theta):
        return bhabha.textbook_M_squared(E_cm, theta)

    sigma_sub = sigma_2to2_cm(M_sq_sub, s, (M_E, M_E), (M_E, M_E),
                              theta_min=tmin, theta_max=tmax)
    sigma_qed = sigma_2to2_cm(M_sq_qed, s, (M_E, M_E), (M_E, M_E),
                              theta_min=tmin, theta_max=tmax)
    assert abs(sigma_sub / sigma_qed - 1.0) < 1e-7


# ----- Decay rate reference -----

def test_muon_decay_textbook_lifetime():
    """Tree-level Gamma = G_F^2 m_mu^5 / (192 pi^3) gives tau_mu within 1% of 2.197 us."""
    Gamma_GeV = gamma_muon_decay_textbook(M_MU)
    tau_s = 1.0 / (Gamma_GeV * 1.519267e24)   # 1 GeV = 1.519e24 / s
    # Tree level predicts 2.187 us; observed 2.197 us -- 0.45% EW radiative correction
    assert abs(tau_s - 2.197e-6) / 2.197e-6 < 1e-2
