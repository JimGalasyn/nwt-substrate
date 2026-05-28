"""Tests for the Feynman-diagram helpers in nwt_substrate.qed.diagram and
nwt_substrate.qcd.diagram.

`qed.diagram` exposes the `Diagram` dataclass + TikZ templates; `qcd.diagram`
exposes per-channel matplotlib renderers + LaTeX expression dict + TikZ
templates. Tests use the matplotlib `Agg` backend so they run headless and
the figures land in memory only.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless — must come before any pyplot import

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from nwt_substrate.qcd.diagram import (
    EXPRESSIONS,
    TIKZ_GG_4VERTEX,
    TIKZ_GG_S,
    TIKZ_QQ_T,
    TIKZ_QQ_U,
    TIKZ_QQBAR_S,
    render_gg_4vertex,
    render_gg_s,
    render_qq_t,
    render_qq_u,
    render_qqbar_s,
)
from nwt_substrate.qed.diagram import (
    TIKZ_BHABHA_S,
    TIKZ_BHABHA_T,
    TIKZ_COMPTON_S,
    TIKZ_COMPTON_U,
    TIKZ_EEMUMU_S,
    TIKZ_MOLLER_T,
    TIKZ_MOLLER_U,
    TIKZ_MUON_DECAY,
    TIKZ_PREAMBLE_NOTE,
    Diagram,
    Term,
)


# ---------------------------------------------------------------
# qed.diagram.Term + Diagram dataclass
# ---------------------------------------------------------------

def test_term_defaults():
    t = Term(latex=r"\bar u")
    assert t.latex == r"\bar u"
    assert t.color == "#1a1a1a"
    assert t.element_id == ""


def test_diagram_repr_contains_process_and_channel():
    d = Diagram(process_name="compton", channel="s_channel", expression="")
    assert "compton" in repr(d)
    assert "s_channel" in repr(d)


def test_diagram_render_without_render_fn_raises():
    d = Diagram(process_name="x", channel="y", expression="")
    with pytest.raises(NotImplementedError, match="No matplotlib renderer"):
        d.render()


def test_diagram_render_with_render_fn_returns_figure():
    """Custom _render_fn is invoked; returns a Figure."""
    calls = []

    def _fn(ax):
        calls.append(ax)
        ax.plot([0, 1], [0, 1])

    d = Diagram(process_name="t", channel="t", expression="",
                _render_fn=_fn)
    fig = d.render()
    assert isinstance(fig, Figure)
    assert len(calls) == 1
    plt.close(fig)


def test_diagram_render_uses_supplied_axes():
    """When ax is provided, the parent Figure is returned (no new fig)."""
    calls = []

    def _fn(ax):
        calls.append(ax)

    d = Diagram(process_name="t", channel="t", expression="", _render_fn=_fn)
    fig0, ax0 = plt.subplots()
    fig = d.render(ax=ax0)
    assert fig is fig0
    plt.close(fig)


def test_diagram_save_writes_png(tmp_path):
    """save() emits a PNG with expression caption."""
    def _fn(ax):
        ax.plot([0, 1], [0, 1])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    d = Diagram(process_name="t", channel="t",
                expression=r"$i\mathcal{M} = 1$",
                _render_fn=_fn)
    out = d.save(tmp_path / "x.png", dpi=72)
    assert out.exists()
    assert out.stat().st_size > 100  # non-trivial PNG


def test_diagram_save_without_expression(tmp_path):
    """include_expression=False skips the caption branch."""
    def _fn(ax):
        ax.plot([0, 1], [0, 1])

    d = Diagram(process_name="t", channel="t", expression="some-expr",
                _render_fn=_fn)
    out = d.save(tmp_path / "y.png", include_expression=False, dpi=72)
    assert out.exists()


def test_diagram_to_tikz_without_template_raises():
    d = Diagram(process_name="t", channel="t", expression="")
    with pytest.raises(NotImplementedError, match="No TikZ template"):
        d.to_tikz()


def test_diagram_to_tikz_returns_template():
    d = Diagram(process_name="t", channel="t", expression="",
                _tikz_template="hello tikz")
    assert d.to_tikz() == "hello tikz"


def test_diagram_to_tikz_writes_to_file(tmp_path):
    d = Diagram(process_name="t", channel="t", expression="",
                _tikz_template="hello-on-disk")
    out_path = tmp_path / "diag.tikz"
    d.to_tikz(file=out_path)
    assert out_path.read_text() == "hello-on-disk"


def test_diagram_render_color_mapped_no_terms_raises():
    """render_color_mapped requires expression_terms populated."""
    d = Diagram(process_name="t", channel="t", expression="")
    with pytest.raises(ValueError, match="no expression_terms"):
        d.render_color_mapped()


def test_diagram_render_color_mapped_partial_axes_raises():
    d = Diagram(process_name="t", channel="t", expression="",
                expression_terms=[Term(latex="x")])
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="both ax_diagram and ax_expression"):
        d.render_color_mapped(ax_diagram=ax, ax_expression=None)
    plt.close(fig)


def test_diagram_render_color_mapped_with_legacy_renderer(tmp_path):
    """Legacy _render_fn that doesn't accept color_map kwarg → fallback path."""
    def _legacy_fn(ax):  # no color_map kwarg
        ax.plot([0, 1], [0, 1])

    d = Diagram(
        process_name="t", channel="t", expression="",
        expression_terms=[Term(latex="x"), Term(latex=r"\gamma_\mu")],
        element_colors={},
        _render_fn=_legacy_fn,
    )
    fig = d.render_color_mapped()
    assert isinstance(fig, Figure)
    plt.close(fig)


# ---------------------------------------------------------------
# qed.diagram TikZ templates carry the preamble note
# ---------------------------------------------------------------

@pytest.mark.parametrize("template", [
    TIKZ_COMPTON_S, TIKZ_COMPTON_U, TIKZ_EEMUMU_S,
    TIKZ_MOLLER_T, TIKZ_MOLLER_U,
    TIKZ_BHABHA_S, TIKZ_BHABHA_T, TIKZ_MUON_DECAY,
])
def test_qed_tikz_templates_have_preamble_and_feynman(template):
    assert TIKZ_PREAMBLE_NOTE in template
    assert r"\begin{feynman}" in template
    assert r"\end{feynman}" in template


# ---------------------------------------------------------------
# qcd.diagram render_* functions (matplotlib smoke)
# ---------------------------------------------------------------

@pytest.mark.parametrize("renderer", [
    render_qqbar_s, render_qq_t, render_qq_u,
    render_gg_s, render_gg_4vertex,
])
def test_qcd_renderer_standalone_returns_figure(renderer):
    """Each renderer called without an axes creates and returns a Figure."""
    fig = renderer()
    assert isinstance(fig, Figure)
    plt.close(fig)


@pytest.mark.parametrize("renderer", [
    render_qqbar_s, render_qq_t, render_qq_u,
    render_gg_s, render_gg_4vertex,
])
def test_qcd_renderer_with_axes_returns_axes(renderer):
    """Each renderer called with axes plots into them and returns the axes."""
    fig, ax = plt.subplots()
    out = renderer(ax)
    assert out is ax
    plt.close(fig)


# ---------------------------------------------------------------
# qcd.diagram constants
# ---------------------------------------------------------------

def test_qcd_expressions_has_three_channels():
    assert set(EXPRESSIONS) == {"qqbar", "qq", "gg"}
    for key, expr in EXPRESSIONS.items():
        assert expr.startswith("$") and expr.endswith("$")
        assert r"\mathcal{M}" in expr


@pytest.mark.parametrize("template", [
    TIKZ_QQBAR_S, TIKZ_QQ_T, TIKZ_QQ_U, TIKZ_GG_S, TIKZ_GG_4VERTEX,
])
def test_qcd_tikz_templates_have_preamble_and_feynman(template):
    assert TIKZ_PREAMBLE_NOTE in template
    assert r"\begin{feynman}" in template
    assert r"\end{feynman}" in template
