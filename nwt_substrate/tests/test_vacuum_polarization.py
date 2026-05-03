"""Tests for substrate-derived 1-loop vacuum polarization and beta-functions."""

import numpy as np
import pytest

from nwt_substrate.algebra.octonions import make_octonion_table
from nwt_substrate.algebra.dirac import make_lorentz_gammas
from nwt_substrate.amplitudes import vacuum_polarization as vp


@pytest.fixture(scope="module")
def gammas():
    T = make_octonion_table()
    return make_lorentz_gammas(T)


# ---------------------------------------------------------------------------
# Structural identity
# ---------------------------------------------------------------------------

def test_substrate_trace_factor_is_8(gammas):
    """The Cl(1,3) loop trace gives factor 8 (= Tr_substrate[I_8]) times
    the textbook tensor (k^mu (k+q)^nu + k^nu (k+q)^mu - eta^{mu nu}(k.(k+q) - m^2))."""
    result = vp.verify_substrate_trace(gammas)
    assert result["trace_factor_matches_dim"]
    assert result["residual"] < 1e-9
    assert result["trace_factor"] == pytest.approx(8.0, abs=1e-9)


def test_structural_factor_equals_one():
    """Substrate (8/4) trace inflation x (1/2) U(1)_em charge projector = 1.
    This is THE clean substrate-derivation result: vacuum polarisation,
    despite the 8-dim spinor space, gives exactly the textbook 4-dim
    Dirac answer."""
    assert vp.SUBSTRATE_TRACE_DIM == 8
    assert vp.DIRAC_TRACE_DIM == 4
    assert vp.SUBSTRATE_PER_DIRAC == 2
    assert vp.U1_EM_CHARGE_PROJECTOR == 0.5
    assert vp.QED_LOOP_STRUCTURAL_FACTOR == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# QED 1-loop running
# ---------------------------------------------------------------------------

def test_qed_pi_R_vanishes_at_zero_q():
    """Pi_R is renormalised to zero at q^2 = 0."""
    val = vp.qed_pi_transverse(q_sq=0.0, m=0.511e-3, Q=1.0)
    assert abs(val) < 1e-12


def test_qed_pi_R_positive_for_spacelike_q():
    """For deep spacelike q^2 < 0: log argument 1 + x(1-x) |q_sq|/m^2 > 1,
    so Pi_R > 0."""
    pi_spacelike = vp.qed_pi_transverse(q_sq=-1.0, m=0.511e-3, Q=1.0)
    assert pi_spacelike > 0


def test_qed_substrate_matches_dirac_for_pi_R():
    """The substrate's (8/4) x (1/2) = 1 structural factor means
    use_substrate=True gives exactly the same answer as use_substrate=False
    (textbook Dirac)."""
    m = 0.511e-3
    for q_sq in [-1.0, -100.0, -1e6 * m ** 2, -1e10 * m ** 2]:
        pi_sub = vp.qed_pi_transverse(q_sq=q_sq, m=m, Q=1.0,
                                        use_substrate=True)
        pi_dir = vp.qed_pi_transverse(q_sq=q_sq, m=m, Q=1.0,
                                        use_substrate=False)
        assert pi_sub == pytest.approx(pi_dir, rel=1e-12)


def test_qed_high_q_asymptote_matches_full_integral():
    """At deep spacelike, the analytic asymptote
       Pi_R -> (alpha/pi) * 1 * [(1/3) log(|q^2|/m^2) - 5/9]
    matches the full Feynman-parameter integral."""
    m = 0.511e-3
    q_sq = -1e10 * m ** 2
    alpha = 1.0 / 137.0
    pi_full = vp.qed_pi_transverse(q_sq=q_sq, m=m, Q=1.0, alpha=alpha)
    pi_asymp = vp.qed_pi_high_q_asymptote(q_sq=q_sq, m=m, Q=1.0, alpha=alpha)
    assert abs(pi_full / pi_asymp - 1.0) < 1e-3


def test_qed_pi_log_slope_extracts_beta_0():
    """Numerically extract the log-derivative coefficient from Pi_R at two
    deep-spacelike scales and verify it equals (alpha Q^2 / pi) * (1/3),
    feeding the substrate-derived beta_0 = (2/3) Q^2 N_c per Dirac."""
    m = 0.511e-3
    alpha = 1.0 / 137.0
    Q = 1.0
    q1 = -1e6 * m ** 2
    q2 = -1e10 * m ** 2
    pi1 = vp.qed_pi_transverse(q_sq=q1, m=m, Q=Q, alpha=alpha)
    pi2 = vp.qed_pi_transverse(q_sq=q2, m=m, Q=Q, alpha=alpha)
    slope = (pi2 - pi1) / np.log(abs(q2) / abs(q1))
    expected = (alpha * Q ** 2 / np.pi) * (1.0 / 3.0)
    assert abs(slope / expected - 1.0) < 1e-4


# ---------------------------------------------------------------------------
# Beta_0^QED matches PDG
# ---------------------------------------------------------------------------

def test_qed_beta_0_per_lepton_is_two_thirds():
    """Substrate-derived beta_0 per charged Dirac lepton = (2/3) * 1 * 1 = 2/3."""
    val = vp.qed_beta_0_per_dirac(Q=1.0, n_color=1)
    assert val == pytest.approx(2.0 / 3.0, abs=1e-12)


def test_qed_beta_0_per_up_quark():
    """Per up-type quark (Q=2/3, N_c=3): beta_0 = (2/3) * 3 * (4/9) = 8/9."""
    val = vp.qed_beta_0_per_dirac(Q=2.0 / 3.0, n_color=3)
    assert val == pytest.approx(8.0 / 9.0, abs=1e-12)


def test_qed_beta_0_per_down_quark():
    """Per down-type quark (Q=1/3, N_c=3): beta_0 = (2/3) * 3 * (1/9) = 2/9."""
    val = vp.qed_beta_0_per_dirac(Q=1.0 / 3.0, n_color=3)
    assert val == pytest.approx(2.0 / 9.0, abs=1e-12)


def test_qed_beta_0_total_matches_pdg_running():
    """Sum over 9 SM Dirac fermions (3 charged leptons + 6 quarks) gives
    beta_0^QED = (2/3) * b_QED with b_QED = 8 (PDG running constant for
    fully-active SM matter at Lambda > m_t):

        3 leptons * 1 * 1   = 3
        3 up   * 4/9 * 3    = 4
        3 down * 1/9 * 3    = 1
        b_QED total          = 8
        beta_0^QED total     = 16/3
    """
    species = vp.standard_qed_species()
    total = vp.qed_beta_0_total(species)
    assert total == pytest.approx(16.0 / 3.0, abs=1e-12)


def test_qed_beta_0_consistent_with_running_couplings_module():
    """The vacuum_polarization module's beta_0 equals (2/3) * b_QED of the
    running_couplings.py thresholds.  Both modules see the same SM matter
    content; the difference is that running_couplings.py has b_QED inputted
    while vacuum_polarization.py derives it from substrate primitives."""
    from nwt_substrate.amplitudes import running_couplings as rc
    species = vp.standard_qed_species()
    beta_0_vp = vp.qed_beta_0_total(species)
    b_qed_rc = sum(t.qed_b for t in rc.standard_thresholds())
    beta_0_rc = (2.0 / 3.0) * b_qed_rc
    assert beta_0_vp == pytest.approx(beta_0_rc, abs=1e-12)


# ---------------------------------------------------------------------------
# Beta_0^QCD
# ---------------------------------------------------------------------------

def test_qcd_quark_loop_per_dirac_is_minus_2_3():
    """Substrate-derived QCD beta_0 per Dirac quark in the loop:
        -(8/4) * (1/2) * (4/3) * T_F  with T_F = 1/2
        = -1 * (4/3) * (1/2) = -2/3
    matching the textbook -(2 T_F) / 3 per quark flavour."""
    val = vp.qcd_beta_0_quark_per_dirac()
    assert val == pytest.approx(-2.0 / 3.0, abs=1e-12)


def test_qcd_beta_0_pure_glue():
    """Pure-glue SU(3) (no quarks): beta_0 = 11 C_A / 3 = 11."""
    val = vp.qcd_beta_0(n_f_dirac=0, N_c=3)
    assert val == pytest.approx(11.0, abs=1e-12)


def test_qcd_beta_0_5flavor_at_M_Z():
    """SU(3), n_f=5 (b, c, s, u, d active at M_Z): beta_0 = 11 - 10/3 = 23/3."""
    val = vp.qcd_beta_0(n_f_dirac=5, N_c=3)
    assert val == pytest.approx(23.0 / 3.0, abs=1e-12)


def test_qcd_beta_0_6flavor_above_top():
    """All 6 flavors active above m_t: beta_0 = 11 - 12/3 = 7."""
    val = vp.qcd_beta_0(n_f_dirac=6, N_c=3)
    assert val == pytest.approx(7.0, abs=1e-12)


def test_qcd_beta_0_consistent_with_running_couplings_module():
    """vacuum_polarization.qcd_beta_0(n_f) matches the (11 - 2 n_f / 3)
    formula hard-coded in running_couplings.py."""
    for n_f in range(7):
        beta_0_vp = vp.qcd_beta_0(n_f_dirac=n_f, N_c=3)
        beta_0_rc = 11.0 - 2.0 * n_f / 3.0
        assert beta_0_vp == pytest.approx(beta_0_rc, abs=1e-12), \
            f"n_f={n_f}: vp={beta_0_vp}, rc={beta_0_rc}"
