"""
QCD diagram rendering — matplotlib (with gluon as helical coil) +
TikZ-Feynman output (with `gluon` line type).
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from ..amplitudes.diagrams import (
    draw_fermion, draw_gluon, draw_photon, draw_vertex,
)
from ..qed.diagram import Diagram, TIKZ_PREAMBLE_NOTE


# ===========================================================================
# Matplotlib renderers for QCD processes (gluons as helical coils)
# ===========================================================================

def render_qqbar_s(ax=None):
    """q qbar -> q' qbar' s-channel via gluon."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))
    v1, v2 = (0.5, 0.65), (0.5, 0.35)
    in_q     = (0.05, 0.85)
    in_qbar  = (0.05, 0.15)
    out_qp   = (0.95, 0.85)
    out_qpbar = (0.95, 0.15)
    draw_fermion(ax, in_q, v1, label=r"$q$")
    draw_fermion(ax, v1, in_qbar, label=r"$\bar{q}$")
    draw_fermion(ax, v2, out_qp, label=r"$q'$")
    draw_fermion(ax, out_qpbar, v2, label=r"$\bar{q}'$")
    draw_gluon(ax, v1, v2, label=r"$g^{*}$", n_loops=4)
    draw_vertex(ax, v1)
    draw_vertex(ax, v2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal"); ax.axis("off")
    if standalone:
        return fig
    return ax


def render_qq_t(ax=None):
    """q q -> q q t-channel."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))
    v1, v2 = (0.4, 0.7), (0.4, 0.3)
    in_v1, in_v2 = (0.05, 0.7), (0.05, 0.3)
    out_v1, out_v2 = (0.95, 0.7), (0.95, 0.3)
    draw_fermion(ax, in_v1, v1, label=r"$q$")
    draw_fermion(ax, v1, out_v1, label=r"$q$")
    draw_fermion(ax, in_v2, v2, label=r"$q$")
    draw_fermion(ax, v2, out_v2, label=r"$q$")
    draw_gluon(ax, v1, v2, label=r"$g^{*}$", n_loops=4)
    draw_vertex(ax, v1); draw_vertex(ax, v2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal"); ax.axis("off")
    if standalone:
        return fig
    return ax


def render_qq_u(ax=None):
    """q q -> q q u-channel (crossed final state)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))
    v1, v2 = (0.4, 0.7), (0.4, 0.3)
    in_v1, in_v2 = (0.05, 0.7), (0.05, 0.3)
    out_v1 = (0.95, 0.3)
    out_v2 = (0.95, 0.7)
    draw_fermion(ax, in_v1, v1, label=r"$q$")
    draw_fermion(ax, v1, out_v1, label=r"$q$")
    draw_fermion(ax, in_v2, v2, label=r"$q$")
    draw_fermion(ax, v2, out_v2, label=r"$q$")
    draw_gluon(ax, v1, v2, label=r"$g^{*}$", n_loops=4)
    draw_vertex(ax, v1); draw_vertex(ax, v2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal"); ax.axis("off")
    if standalone:
        return fig
    return ax


def render_gg_s(ax=None):
    """gg -> gg s-channel (3-gluon × 3-gluon)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))
    v1, v2 = (0.5, 0.65), (0.5, 0.35)
    in1, in2 = (0.05, 0.85), (0.05, 0.15)
    out1, out2 = (0.95, 0.85), (0.95, 0.15)
    draw_gluon(ax, in1, v1, label=r"$g$", n_loops=4)
    draw_gluon(ax, in2, v1, label=r"$g$", n_loops=4)
    draw_gluon(ax, v1, v2, label=r"$g^{*}$", n_loops=3)
    draw_gluon(ax, v2, out1, label=r"$g$", n_loops=4)
    draw_gluon(ax, v2, out2, label=r"$g$", n_loops=4)
    draw_vertex(ax, v1); draw_vertex(ax, v2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal"); ax.axis("off")
    if standalone:
        return fig
    return ax


def render_gg_4vertex(ax=None):
    """gg -> gg 4-gluon contact term (single vertex)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))
    v = (0.5, 0.5)
    in1, in2 = (0.1, 0.85), (0.1, 0.15)
    out1, out2 = (0.9, 0.85), (0.9, 0.15)
    draw_gluon(ax, in1, v, label=r"$g$", n_loops=4)
    draw_gluon(ax, in2, v, label=r"$g$", n_loops=4)
    draw_gluon(ax, v, out1, label=r"$g$", n_loops=4)
    draw_gluon(ax, v, out2, label=r"$g$", n_loops=4)
    draw_vertex(ax, v, size=18)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal"); ax.axis("off")
    if standalone:
        return fig
    return ax


# ===========================================================================
# LaTeX expressions for the substrate-algebra QCD amplitudes
# ===========================================================================

EXPRESSIONS = {
    "qqbar": (
        r"$i\mathcal{M} = \frac{-i\,g_s^2}{s}\,(T^a)_{i'i}(T^a)_{jj'}\,"
        r"[\bar{v}_2\gamma^{\mu}u_1][\bar{u}_3\gamma_{\mu}v_4]$"
    ),
    "qq": (
        r"$i\mathcal{M} = \frac{-i\,g_s^2}{t}\,(T^a)_{i'i}(T^a)_{j'j}\,"
        r"[\bar{u}_3\gamma^{\mu}u_1][\bar{u}_4\gamma_{\mu}u_2]"
        r"\;-\;(t \leftrightarrow u, \; 3 \leftrightarrow 4)$"
    ),
    "gg": (
        r"$i\mathcal{M} = i\mathcal{M}_s + i\mathcal{M}_t + i\mathcal{M}_u + i\mathcal{M}_4$"
    ),
}


# ===========================================================================
# TikZ-Feynman templates (gluon = helical coil)
# ===========================================================================

TIKZ_QQBAR_S = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [above left=of a]  (i1) {\(q\)};
    \vertex [below left=of a]  (i2) {\(\bar{q}\)};
    \vertex [above right=of b] (f1) {\(q'\)};
    \vertex [below right=of b] (f2) {\(\bar{q}'\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (i2),
      (a)  -- [gluon, edge label=\(g^{*}\)] (b),
      (f2) -- [fermion] (b) -- [fermion] (f1),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_QQ_T = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [left=of a]  (i1) {\(q\)};
    \vertex [right=of a] (f1) {\(q\)};
    \vertex [left=of b]  (i2) {\(q\)};
    \vertex [right=of b] (f2) {\(q\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (f1),
      (a)  -- [gluon, edge label=\(g^{*}\)] (b),
      (i2) -- [fermion] (b) -- [fermion] (f2),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_QQ_U = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [left=of a]  (i1) {\(q\)};
    \vertex [right=of a] (f1) {\(q\)};
    \vertex [left=of b]  (i2) {\(q\)};
    \vertex [right=of b] (f2) {\(q\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (f2),
      (a)  -- [gluon, edge label=\(g^{*}\)] (b),
      (i2) -- [fermion] (b) -- [fermion] (f1),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_GG_S = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [right=2cm of a] (b);
    \vertex [above left=of a]  (i1) {\(g\)};
    \vertex [below left=of a]  (i2) {\(g\)};
    \vertex [above right=of b] (f1) {\(g\)};
    \vertex [below right=of b] (f2) {\(g\)};
    \diagram* {
      (i1) -- [gluon] (a),
      (i2) -- [gluon] (a),
      (a)  -- [gluon, edge label=\(g^{*}\)] (b),
      (b)  -- [gluon] (f1),
      (b)  -- [gluon] (f2),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_GG_4VERTEX = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (v);
    \vertex [above left=of v]  (i1) {\(g\)};
    \vertex [below left=of v]  (i2) {\(g\)};
    \vertex [above right=of v] (f1) {\(g\)};
    \vertex [below right=of v] (f2) {\(g\)};
    \diagram* {
      (i1) -- [gluon] (v),
      (i2) -- [gluon] (v),
      (v)  -- [gluon] (f1),
      (v)  -- [gluon] (f2),
    };
  \end{feynman}
\end{tikzpicture}
"""
