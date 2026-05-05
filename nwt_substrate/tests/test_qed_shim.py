"""Tests for the nwt.qed shim — community-style QED interface."""

import numpy as np
import pytest

import nwt_substrate.qed as qed


# ---------- Constants ----------

def test_alpha_thomson():
    assert abs(qed.alpha - 1.0 / 137.035999084) < 1e-15


def test_lepton_masses():
    assert abs(qed.m_e - 0.510998928e-3) < 1e-15
    assert abs(qed.m_mu - 0.10565837) < 1e-15


def test_alpha_running():
    # alpha grows with energy: alpha(M_Z) > alpha(0)
    assert qed.alpha_at(qed.m_Z) > qed.alpha


# ---------- Mandelstam helper ----------

def test_mandelstam_s_t_u_consistency():
    # Massless 2->2 in CM: s + t + u = 0
    E = 5.0
    theta = 1.0
    p1 = np.array([E, 0.0, 0.0, E])
    p2 = np.array([E, 0.0, 0.0, -E])
    p3 = np.array([E, E * np.sin(theta), 0.0, E * np.cos(theta)])
    p4 = p1 + p2 - p3
    s, t, u = qed.mandelstam(p1, p2, p3, p4)
    assert abs(s + t + u) < 1e-9


# ---------- Compton ----------

def test_compton_cross_section_in_pb():
    """sigma in pb at low omega ~ Thomson cross-section ~ 6.65e11 pb (= 0.665 barn)."""
    sigma_pb = qed.compton.sigma_total(omega=1e-7)
    assert sigma_pb > 0
    assert 1e11 < sigma_pb < 1e12   # Thomson is 6.65e11 pb = 0.665 b


def test_compton_thomson_limit_property():
    """Thomson sigma = (8 pi/3) r_e^2 ≈ 6.65×10^11 pb."""
    sigma = qed.compton.thomson_limit
    assert 1e11 < sigma < 1e12


def test_compton_event():
    event = qed.compton.event(omega=0.01, theta=1.0)
    assert event.process == "compton"
    assert event.M_sq_avg > 0
    assert event.dsigma_dOmega_pb > 0


# ---------- e+e- -> mu+mu- ----------

def test_eemumu_sigma_at_10_GeV():
    """At E_cm = 10 GeV, sigma ~ 870 pb (textbook 4 pi alpha^2 / 3s value)."""
    sigma_pb = qed.eemumu.sigma_total(E_cm=10.0)
    # Textbook: 4 pi (1/137)^2 / (3 * 100 GeV^2) * 0.3894e9 pb/GeV^-2
    expected_pb = 4 * np.pi * qed.alpha**2 / (3 * 100.0) * 3.894e8
    assert abs(sigma_pb / expected_pb - 1.0) < 1e-6


def test_eemumu_event_returns_proper_object():
    event = qed.eemumu.event(E_cm=10.0, theta=1.0)
    assert event.process == "eemumu"
    assert event.s == pytest.approx(100.0, rel=1e-6)
    assert event.sigma_total_pb is not None
    assert event.sigma_total_pb > 0


# ---------- Møller and Bhabha ----------

def test_moller_dsigma_pb_positive():
    ds = qed.moller.dsigma_dOmega(E_cm=10.0, theta=np.pi / 2)
    assert ds > 0


def test_bhabha_sigma_finite():
    sigma_pb = qed.bhabha.sigma_total(E_cm=10.0, cutoff_deg=30.0)
    assert sigma_pb > 0
    assert sigma_pb < 1e10   # sanity bound


# ---------- Muon decay ----------

def test_muon_lifetime():
    tau_us = qed.muon_decay.lifetime()
    assert abs(tau_us - 2.187) < 0.01  # ~2.187 us tree-level


def test_muon_michel_spectrum_normalized():
    """Integral of dGamma/dx_e = Gamma."""
    from scipy.integrate import quad
    Gamma = qed.muon_decay.Gamma()
    integral, _ = quad(lambda x: qed.muon_decay.michel_spectrum(x), 0, 1)
    assert abs(integral / Gamma - 1.0) < 1e-9


# ---------- Diagrams ----------

def test_diagram_attributes():
    d = qed.compton.diagrams.s_channel
    assert d.process_name == "compton"
    assert d.channel == "s_channel"
    assert "mathcal{M}" in d.expression


def test_diagram_renders():
    """diagram.render() returns a matplotlib Figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = qed.compton.diagrams.s_channel.render()
    assert fig is not None
    plt.close(fig)


def test_diagram_to_tikz():
    """diagram.to_tikz() returns valid TikZ-Feynman code."""
    code = qed.compton.diagrams.s_channel.to_tikz()
    assert "feynman" in code
    assert "fermion" in code
    assert "photon" in code
    assert "tikz-feynman" in code   # in the comment header


def test_compton_diagrams_have_color_mapping():
    """Compton s/u channels expose Term lists + element_colors that
    map each algebraic chunk of iM to a diagram element."""
    for diag in (qed.compton.diagrams.s_channel,
                 qed.compton.diagrams.u_channel):
        assert diag.expression_terms, (
            f"{diag.process_name}/{diag.channel} has no expression_terms"
        )
        # Each element_id appearing in a Term must have a color
        # registered in element_colors (excluding the empty-id
        # prefactor/bracket terms).
        for term in diag.expression_terms:
            if term.element_id:
                assert term.element_id in diag.element_colors, (
                    f"term {term.element_id!r} has no color in"
                    f" {diag.process_name}/{diag.channel}"
                )
        # The five Compton element_ids must all be present and
        # external-electron / external-photon colors must agree across
        # the in/out pair (same particle line, same color).
        ec = diag.element_colors
        for key in ("in_e", "in_g", "internal_e", "out_e", "out_g"):
            assert key in ec, f"missing color for {key!r}"
        assert ec["in_e"] == ec["out_e"], (
            "incoming and outgoing electron should share the same color"
        )
        assert ec["in_g"] == ec["out_g"], (
            "incoming and outgoing photon should share the same color"
        )


def test_compton_render_color_mapped_returns_figure():
    """Diagram.render_color_mapped returns a Figure with two axes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = qed.compton.diagrams.s_channel.render_color_mapped()
    assert fig is not None
    # Two stacked axes: diagram on top, expression on bottom
    assert len(fig.axes) == 2
    plt.close(fig)


def test_render_color_mapped_raises_when_no_terms():
    """A Diagram without expression_terms must raise on color-mapped
    render rather than silently producing an unlabeled figure."""
    import pytest
    from nwt_substrate.qed.diagram import Diagram
    bare = Diagram(
        process_name="bare", channel="test",
        expression="$x$",
        _render_fn=lambda ax, color_map=None: None,
    )
    with pytest.raises(ValueError, match="expression_terms"):
        bare.render_color_mapped()


def test_diagram_to_tikz_with_file(tmp_path):
    out = tmp_path / "compton_s.tex"
    qed.compton.diagrams.s_channel.to_tikz(file=out)
    assert out.exists()
    content = out.read_text()
    assert "feynman" in content


def test_all_processes_have_diagrams():
    """Every canonical process exposes at least one diagram."""
    for proc in [qed.compton, qed.eemumu, qed.moller, qed.bhabha,
                 qed.muon_decay]:
        attrs = [a for a in dir(proc.diagrams) if not a.startswith("_")]
        assert len(attrs) >= 1, f"{proc.name} has no diagrams"


def test_gallery_all_produces_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = qed.gallery_all()
    assert fig is not None
    plt.close(fig)
