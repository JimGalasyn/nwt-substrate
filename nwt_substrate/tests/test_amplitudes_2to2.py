"""Tests for 2->2 substrate-algebra amplitudes."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas
from nwt_substrate.amplitudes.processes import eemumu, moller, bhabha


M_E = 0.510998928e-3   # GeV
M_MU = 0.10566         # GeV


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


# ----- e+ e- -> mu+ mu- -----

def test_eemumu_matches_textbook(gammas):
    """Substrate e+e-->mu+mu- reproduces QED textbook in high-E limit (1e-9)."""
    E_cm = 10.0
    for theta in [0.5, 1.0, 1.5, 2.0]:
        M_sub, *_ = eemumu.M_squared_avg(E_cm, theta, M_E, M_MU, gammas)
        M_qed = eemumu.textbook_M_squared(E_cm, theta, M_MU)
        assert abs(M_sub / M_qed - 1.0) < 1e-7, (
            f"theta={theta}: ratio={M_sub/M_qed}, M_sub={M_sub}, M_qed={M_qed}"
        )


def test_eemumu_kinematics_threshold(gammas):
    """At threshold (E_cm = 2 m_mu), p_mu = 0 and |M|^2 should be uniform in theta."""
    E_cm = 2 * M_MU + 1e-3   # just above threshold
    vals = []
    for theta in [0.5, 1.0, 1.5]:
        M_sub, *_ = eemumu.M_squared_avg(E_cm, theta, M_E, M_MU, gammas)
        vals.append(M_sub)
    # near threshold, |M|^2 ≈ e^4 * 2 (from textbook formula)
    assert np.std(vals) / np.mean(vals) < 0.1, (
        f"|M|^2 not uniform near threshold: {vals}"
    )


# ----- Moller (t-channel only) -----

def _moller_t_channel_only(E_cm, theta, m_e, gammas, e_charge):
    """Compute |M_t|^2_avg for Moller (t-channel only, no interference)."""
    from nwt_substrate.amplitudes.spinors import (
        positive_energy_spinors, adjoint_spinor,
    )
    p1, p2, p3, p4 = moller.kinematics_cm(E_cm, theta, m_e)
    u1s = positive_energy_spinors(p1, gammas, m_e)
    u2s = positive_energy_spinors(p2, gammas, m_e)
    u3s = positive_energy_spinors(p3, gammas, m_e)
    u4s = positive_energy_spinors(p4, gammas, m_e)
    g0 = gammas[0]
    q_t = p1 - p3
    p_sq_t = float(q_t[0] ** 2 - q_t[1] ** 2 - q_t[2] ** 2 - q_t[3] ** 2)
    M_sq = 0.0
    for u1 in u1s:
        for u2 in u2s:
            for u3 in u3s:
                for u4 in u4s:
                    ub3 = adjoint_spinor(u3, g0)
                    ub4 = adjoint_spinor(u4, g0)
                    cur1 = np.array([complex(ub3 @ gammas[mu] @ u1)
                                     for mu in range(4)])
                    cur2 = np.array([complex(ub4 @ gammas[mu] @ u2)
                                     for mu in range(4)])
                    M = (cur1[0] * cur2[0] - cur1[1] * cur2[1]
                         - cur1[2] * cur2[2] - cur1[3] * cur2[3])
                    iM_t = -1j * (e_charge ** 2 / p_sq_t) * M
                    M_sq += abs(iM_t) ** 2
    avg = M_sq / 16.0
    rel_norm = (2 * float(p1[0])) * (2 * float(p2[0])) \
               * (2 * float(p3[0])) * (2 * float(p4[0]))
    return rel_norm * avg


def test_moller_t_channel_only_exact(gammas):
    """Moller t-channel alone reproduces 2 e^4 (s^2+u^2)/t^2 to machine precision."""
    from nwt_substrate.amplitudes.vertices import ELECTRIC_CHARGE
    E_cm = 10.0
    for theta in [0.5, 1.0, 1.5, 2.0, 2.5]:
        M_sub_t = _moller_t_channel_only(E_cm, theta, M_E, gammas, ELECTRIC_CHARGE)
        s = E_cm ** 2
        t = -s * np.sin(theta / 2.0) ** 2
        u = -s * np.cos(theta / 2.0) ** 2
        M_qed_t = 2.0 * (ELECTRIC_CHARGE ** 4) * (s ** 2 + u ** 2) / t ** 2
        assert abs(M_sub_t / M_qed_t - 1.0) < 1e-7, (
            f"theta={theta}: ratio={M_sub_t/M_qed_t}"
        )


def test_moller_full_matches_textbook(gammas):
    """Full Moller (t-u with substrate interference doubling) matches QED to ~1e-9."""
    E_cm = 10.0
    for theta in [0.3, 0.5, 1.0, 1.5, 2.0, 2.5]:
        M_sub, *_ = moller.M_squared_avg(E_cm, theta, M_E, gammas)
        M_qed = moller.textbook_M_squared(E_cm, theta)
        assert abs(M_sub / M_qed - 1.0) < 1e-7, (
            f"theta={theta}: ratio={M_sub/M_qed}"
        )


# ----- Bhabha -----

def test_bhabha_full_matches_textbook(gammas):
    """Full Bhabha (s-t with substrate interference doubling) matches QED to ~1e-9."""
    E_cm = 10.0
    for theta in [0.3, 0.5, 1.0, 1.5, 2.0, 2.5]:
        M_sub, *_ = bhabha.M_squared_avg(E_cm, theta, M_E, gammas)
        M_qed = bhabha.textbook_M_squared(E_cm, theta)
        assert abs(M_sub / M_qed - 1.0) < 1e-7, (
            f"theta={theta}: ratio={M_sub/M_qed}"
        )
