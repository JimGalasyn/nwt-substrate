"""Tests for QCD same-flavor q q -> q q amplitude (with substrate doubling)."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas
from nwt_substrate.amplitudes.processes import qcd_qq


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


def test_qcd_qq_color_factors():
    """Diagonal color factors are 2/9, interference is -2/27."""
    assert abs(qcd_qq.COLOR_T - 2.0 / 9.0) < 1e-15
    assert abs(qcd_qq.COLOR_U - 2.0 / 9.0) < 1e-15
    assert abs(qcd_qq.COLOR_CROSS + 2.0 / 27.0) < 1e-15


def test_qcd_qq_matches_textbook(gammas):
    """qq->qq matches textbook QCD to ~1e-9 (using substrate interference doubling)."""
    m_q = 2.16e-3
    alpha_s = 0.1179
    E_cm = 100.0
    for theta in [0.3, 0.5, 1.0, 1.5, 2.0, 2.5]:
        M_sub, *_ = qcd_qq.M_squared_avg(E_cm, theta, m_q, gammas, alpha_s)
        M_qcd = qcd_qq.textbook_M_squared(E_cm, theta, alpha_s)
        assert abs(M_sub / M_qcd - 1.0) < 1e-7, (
            f"theta={theta}: ratio={M_sub/M_qcd}"
        )


def test_qcd_qq_interference_sign():
    """QCD interference is negative (anti-screening) at theta=pi/2."""
    # At theta=pi/2: t = u = -s/2, so the interference term -(2/3) s^2/(tu)
    # = -(2/3) * (s^2 / (s^2/4)) = -(2/3) * 4 = -8/3 (in units of g_s^4 * (4/9))
    # Total: |M|^2 = (4/9) * [(s^2 + s^2/4)/(s^2/4) + (s^2 + s^2/4)/(s^2/4) - 8/3]
    #             = (4/9) * [5 + 5 - 8/3] = (4/9) * 22/3 = 88/27
    alpha_s = 0.1
    E_cm = 100.0
    theta = np.pi / 2
    M_qcd = qcd_qq.textbook_M_squared(E_cm, theta, alpha_s)
    g_s_4 = (4.0 * np.pi * alpha_s) ** 2
    expected = (88.0 / 27.0) * g_s_4
    assert abs(M_qcd / expected - 1.0) < 1e-12
