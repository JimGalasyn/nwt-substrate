"""Tests for the nwt.qcd shim."""

import numpy as np
import pytest

import nwt_substrate.qcd as qcd


# ---------- Constants ----------

def test_alpha_s_at_mz():
    assert abs(qcd.alpha_s - 0.1179) < 1e-6


def test_running_decreases_at_high_mu():
    """Asymptotic freedom: alpha_s falls with energy."""
    assert qcd.alpha_s_at(1000.0) < qcd.alpha_s_at(100.0)


def test_color_constants():
    assert qcd.C_F == 4.0 / 3.0
    assert qcd.C_A == 3.0
    assert qcd.N_c == 3


# ---------- Color algebra ----------

def test_gell_mann_count():
    assert len(qcd.gell_mann()) == 8


def test_T_orthonormal():
    Tlist = qcd.T()
    for a in range(8):
        for b in range(8):
            tr = np.trace(Tlist[a] @ Tlist[b]).real
            expected = 0.5 if a == b else 0.0
            assert abs(tr - expected) < 1e-12


def test_structure_constants_known_values():
    farr = qcd.f()
    # f^123 = 1
    assert abs(farr[0, 1, 2] - 1.0) < 1e-12
    # f^458 = sqrt(3)/2
    assert abs(farr[3, 4, 7] - np.sqrt(3) / 2) < 1e-12


# ---------- Processes ----------

def test_qqbar_event():
    e = qcd.qqbar.event(E_cm=10.0, theta=1.0)
    assert e.process == "qqbar"
    assert e.s == pytest.approx(100.0, rel=1e-6)
    assert e.dsigma_dOmega_pb > 0


def test_qqbar_sigma_at_10_GeV_finite_pb():
    sigma_pb = qcd.qqbar.sigma_total(E_cm=10.0)
    assert sigma_pb > 0
    assert 1e2 < sigma_pb < 1e6   # rough scale check


def test_qq_dsigma_finite():
    ds = qcd.qq.dsigma_dOmega(E_cm=10.0, theta=1.5)
    assert ds > 0


def test_gg_textbook_at_pi_over_2():
    """At theta = pi/2: |M|^2 = (243/8) g_s^4."""
    M = qcd.gg.textbook_M_squared(E_cm=100.0, theta=np.pi / 2, alpha_s=0.1)
    g_s_4 = (4 * np.pi * 0.1) ** 2
    expected = (243.0 / 8.0) * g_s_4
    assert abs(M / expected - 1.0) < 1e-12


def test_gg_dsigma_textbook():
    ds = qcd.gg.dsigma_dOmega_textbook(E_cm=100.0, theta=1.0)
    assert ds > 0


# ---------- Diagrams ----------

def test_qqbar_diagram_render():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = qcd.qqbar.diagrams.s_channel.render()
    plt.close(fig)


def test_qq_diagram_tikz():
    code = qcd.qq.diagrams.t_channel.to_tikz()
    assert "gluon" in code
    assert "feynman" in code


def test_gg_diagram_tikz_uses_gluon_lines():
    """gg diagrams should use 'gluon' line type, not 'photon'."""
    code = qcd.gg.diagrams.s_channel.to_tikz()
    assert "gluon" in code
    assert "photon" not in code


def test_all_qcd_processes_have_diagrams():
    for proc in [qcd.qqbar, qcd.qq, qcd.gg]:
        attrs = [a for a in dir(proc.diagrams) if not a.startswith("_")]
        assert len(attrs) >= 1


def test_gallery_all_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = qcd.gallery_all()
    assert fig is not None
    plt.close(fig)
