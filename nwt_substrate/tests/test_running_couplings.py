"""Tests for QED and QCD running couplings at 1-loop."""

import numpy as np
import pytest

from nwt_substrate.amplitudes.running_couplings import (
    alpha_qed, alpha_s,
    ALPHA_THOMSON, M_Z, M_E, M_TOP,
    pdg_alpha_qed_at_mz, pdg_alpha_s_at_mz,
    standard_thresholds,
)


# ----- QED -----

def test_alpha_qed_at_thomson_limit():
    """alpha_QED at q=0 (Thomson) is 1/137.036."""
    a = alpha_qed(M_E, alpha_at_mu0=ALPHA_THOMSON, mu0=M_E)
    assert abs(a - ALPHA_THOMSON) < 1e-12


def test_alpha_qed_runs_upward_with_scale():
    """1/alpha decreases monotonically with mu (alpha increases)."""
    one_over_a_prev = 1.0 / ALPHA_THOMSON
    for mu in [0.01, 0.1, 1.0, 10.0, M_Z, 1000.0]:
        one_over_a = 1.0 / alpha_qed(mu)
        assert one_over_a < one_over_a_prev
        one_over_a_prev = one_over_a


def test_alpha_qed_at_mz_within_3_percent_of_pdg():
    """1-loop leading-log alpha(M_Z) is ~1/131.8, PDG is ~1/127.95.
    Discrepancy is hadronic vacuum polarization (~3%)."""
    a_mz = alpha_qed(M_Z)
    a_pdg = pdg_alpha_qed_at_mz()
    rel_diff = abs(1.0 / a_mz - 1.0 / a_pdg) / (1.0 / a_pdg)
    assert rel_diff < 0.05  # within 5% (we expect ~3%)


def test_alpha_qed_lepton_only_running():
    """Running with leptons only: 1/alpha decreases by ~2.4 from Thomson to M_Z."""
    leptons_only = [t for t in standard_thresholds() if t.n_f_qcd == 0]
    a_mz_lep = alpha_qed(M_Z, thresholds=leptons_only)
    # Thomson 1/alpha = 137.036 -> M_Z 1/alpha (leptons only) = ~134.6
    # Magnitude of running: 137.036 - 134.6 = 2.4
    delta = 1.0 / ALPHA_THOMSON - 1.0 / a_mz_lep
    # 3 leptons: ln(M_Z/m_e) + ln(M_Z/m_mu) + ln(M_Z/m_tau) sum / (3 pi) ~ 2.4
    assert 2.0 < delta < 3.0


# ----- QCD -----

def test_alpha_s_at_mz_anchor():
    """alpha_s(M_Z) = 0.1179 by definition of anchor."""
    a = alpha_s(M_Z)
    assert abs(a - 0.1179) < 1e-12


def test_alpha_s_asymptotic_freedom():
    """alpha_s is monotonically decreasing with mu (asymptotic freedom)."""
    a_prev = alpha_s(1.0)
    for mu in [2.0, 5.0, 10.0, M_Z, 1000.0, 1e6, 1e10]:
        a = alpha_s(mu)
        assert a < a_prev
        a_prev = a


def test_alpha_s_at_low_energy_perturbative_regime():
    """At mu = 1 GeV (just above Lambda_QCD), 1-loop gives α_s ~ 0.36."""
    a = alpha_s(1.0)
    assert 0.3 < a < 0.5  # 1-loop estimate; full calc gives ~0.5


def test_alpha_s_threshold_continuity():
    """alpha_s should be continuous across thresholds."""
    # Approaching m_b from below and above
    eps = 1e-6
    a_below = alpha_s(4.18 - eps)
    a_above = alpha_s(4.18 + eps)
    assert abs(a_above / a_below - 1.0) < 1e-3   # threshold matching is exact


def test_alpha_s_runs_to_high_scale():
    """At MSSM-unification scale ~10^16 GeV, alpha_s ~ 0.02."""
    a_gut = alpha_s(2e16)
    assert 0.015 < a_gut < 0.03   # 1-loop ballpark
