"""Tests for QCD q qbar -> q' qbar' substrate amplitude."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas
from nwt_substrate.amplitudes.processes import qcd_qqbar


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


def test_color_factor():
    """Color factor for q qbar -> q' qbar' s-channel is 2/9."""
    assert abs(qcd_qqbar.color_factor_qqbar() - 2.0 / 9.0) < 1e-12


def test_qqbar_matches_textbook(gammas):
    """Substrate q qbar -> q' qbar' matches QCD textbook to ~1e-9."""
    m_q = 2.16e-3
    m_qprime = 4.67e-3
    alpha_s = 0.1179
    E_cm = 100.0
    for theta in [0.3, 0.5, 1.0, 1.5, 2.0, 2.5]:
        M_sub, *_ = qcd_qqbar.M_squared_avg(E_cm, theta, m_q, m_qprime,
                                             gammas, alpha_s)
        M_qcd = qcd_qqbar.textbook_M_squared(E_cm, theta, m_qprime, alpha_s)
        assert abs(M_sub / M_qcd - 1.0) < 1e-7, (
            f"theta={theta}: ratio={M_sub/M_qcd}"
        )


def test_qcd_has_correct_color_relative_to_qed(gammas):
    """At equal couplings (alpha_s = alpha), |M|^2(qqbar->q'q'bar) = (2/9) * |M|^2(e+e-->mu+mu-)."""
    from nwt_substrate.amplitudes.processes import eemumu

    m_lep = 0.10566  # use muon mass for both substrate calls
    alpha_eq = 0.1
    e_eq = float(np.sqrt(4.0 * np.pi * alpha_eq))
    E_cm = 100.0
    theta = 1.0

    M_qed, *_ = eemumu.M_squared_avg(E_cm, theta, 0.510998928e-3, m_lep,
                                      gammas, e=e_eq)
    M_qcd, *_ = qcd_qqbar.M_squared_avg(E_cm, theta, 0.510998928e-3, m_lep,
                                         gammas, alpha_s=alpha_eq)
    ratio = M_qcd / M_qed
    expected = 2.0 / 9.0
    assert abs(ratio - expected) / expected < 1e-7
