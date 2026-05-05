"""
Diagram class — a uniform interface over the matplotlib renderers in
`nwt_substrate.amplitudes.diagrams`, plus TikZ-Feynman output for papers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import matplotlib.pyplot as plt


@dataclass
class Term:
    """
    A single chunk of a Feynman amplitude expression, paired with the
    diagram element it represents.

    `latex` is the LaTeX fragment of the chunk (no surrounding $...$);
    `color` is a matplotlib color string (hex or named); `element_id`
    is a stable identifier matching one of the keys understood by the
    matplotlib renderer's `color_map` argument
    (e.g., `"in_e"`, `"internal_e"`, `"out_g"`).  Terms whose
    `element_id` is empty (the default, e.g. for prefactors and
    brackets) carry color but are not tied to any diagram part; they
    are used to give the expression a consistent visual rhythm.
    """

    latex: str
    color: str = "#1a1a1a"
    element_id: str = ""


@dataclass
class Diagram:
    """
    A Feynman diagram in the QED shim.

    Wraps a matplotlib renderer (from nwt_substrate.amplitudes.diagrams) with
    a uniform interface and TikZ-Feynman output for paper inclusion.

    Attributes
    ----------
    process_name : str
        Canonical process name ("compton", "eemumu", etc.).  Useful for
        looking up the corresponding amplitude calculator.
    channel : str
        Channel label ("s_channel", "t_channel", "tree", etc.).
    expression : str
        LaTeX string of the substrate-algebra amplitude i M (single
        unbroken chunk; used for plain `save()` rendering).
    feynman_rules : dict
        Metadata: which propagators / vertices appear in this diagram.
    expression_terms : list[Term]
        Optional term-by-term decomposition of the amplitude, with
        each term tagged by its diagram element (for the
        `render_color_mapped()` view that maps algebra to picture).
    element_colors : dict[str, str]
        Optional mapping from element_id to matplotlib color, used by
        the renderer to color each diagram part to match its
        corresponding term in `expression_terms`.
    _render_fn : callable
        Internal: matplotlib renderer that draws the diagram on an Axes.
        Renderers that support color mapping accept a `color_map=...`
        keyword argument; legacy renderers ignore it.
    _tikz_template : str
        Internal: TikZ-Feynman LaTeX template (idiomatic, lets tikz-feynman
        do auto-layout).
    """

    process_name: str
    channel: str
    expression: str
    feynman_rules: dict = field(default_factory=dict)
    expression_terms: list = field(default_factory=list)
    element_colors: dict = field(default_factory=dict)
    _render_fn: Callable = field(default=None, repr=False)
    _tikz_template: str = field(default="", repr=False)

    # ---------- Matplotlib rendering ----------

    def render(self, ax=None, figsize=(5, 4)) -> "plt.Figure":
        """
        Render the diagram into a matplotlib Axes.  If `ax` is None, creates
        a new Figure and returns it; otherwise plots into the given axes
        and returns the parent Figure.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        if self._render_fn is None:
            raise NotImplementedError(
                f"No matplotlib renderer attached to {self.process_name}/{self.channel}"
            )
        self._render_fn(ax)
        return fig

    def save(self, path: str | Path, dpi: int = 150,
             include_expression: bool = True, **kwargs) -> Path:
        """
        Render and save to PNG/PDF/SVG.  If `include_expression`, prints the
        substrate-algebra LaTeX underneath the diagram as a caption.
        """
        fig = self.render()
        if include_expression and self.expression:
            fig.text(0.5, 0.02, self.expression, ha="center", fontsize=10)
            fig.subplots_adjust(bottom=0.18)
        path = Path(path)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
        plt.close(fig)
        return path

    # ---------- Color-mapped rendering: algebra <-> picture ----------

    def render_color_mapped(
        self,
        ax_diagram=None,
        ax_expression=None,
        figsize: tuple[float, float] = (6.0, 5.6),
        height_ratios: tuple[float, float] = (8.0, 1.0),
        expression_fontsize: int = 12,
    ) -> "plt.Figure":
        """
        Render the diagram with each part colored by its `element_id`,
        and below it render the amplitude as a horizontal chain of
        colored chunks (one per `Term` in `expression_terms`), so the
        algebra-to-picture correspondence is visually explicit.

        If `ax_diagram` and `ax_expression` are both `None`, creates a
        new Figure with two stacked axes and returns it.  Otherwise
        the caller supplies axes (e.g., from a custom GridSpec) and
        the existing Figure is returned.

        Requires that `expression_terms` and `element_colors` are
        populated (raises `ValueError` if not), and that the
        `_render_fn` accepts a `color_map` keyword argument.
        """
        from matplotlib.offsetbox import (
            TextArea, HPacker, AnchoredOffsetbox,
        )

        if not self.expression_terms:
            raise ValueError(
                f"Diagram {self.process_name}/{self.channel} has no "
                f"expression_terms; populate them in the process "
                f"definition before calling render_color_mapped()."
            )

        if ax_diagram is None and ax_expression is None:
            fig = plt.figure(figsize=figsize)
            gs = fig.add_gridspec(
                2, 1, height_ratios=list(height_ratios),
                hspace=0.02, left=0.02, right=0.98,
                top=0.99, bottom=0.02,
            )
            ax_diagram = fig.add_subplot(gs[0])
            ax_expression = fig.add_subplot(gs[1])
        elif ax_diagram is None or ax_expression is None:
            raise ValueError(
                "Provide both ax_diagram and ax_expression, or neither."
            )
        else:
            fig = ax_diagram.get_figure()

        # Render the diagram with the per-element color map.  The
        # renderer must accept a `color_map=...` keyword.
        try:
            self._render_fn(ax_diagram, color_map=self.element_colors)
        except TypeError:
            # Legacy renderer: fall back to default colors.
            self._render_fn(ax_diagram)
        ax_diagram.set_aspect("equal")
        ax_diagram.axis("off")

        # Render the expression as a chain of colored chunks.  We use
        # an AnchoredOffsetbox with an HPacker so each Term can have
        # its own color (matplotlib MathText does not support inline
        # \color, but separate Text artists can be packed horizontally).
        children = [
            TextArea(
                t.latex if t.latex.startswith("$") else f"${t.latex}$",
                textprops=dict(color=t.color, fontsize=expression_fontsize),
            )
            for t in self.expression_terms
        ]
        hpacker = HPacker(children=children, align="baseline",
                          pad=0, sep=2)
        anchored = AnchoredOffsetbox(
            loc="center", child=hpacker, pad=0, frameon=False,
            bbox_to_anchor=(0.5, 0.5),
            bbox_transform=ax_expression.transAxes,
        )
        ax_expression.add_artist(anchored)
        ax_expression.axis("off")

        return fig

    def save_color_mapped(
        self, path: str | Path, dpi: int = 150, **kwargs,
    ) -> Path:
        """Render via `render_color_mapped()` and save to disk."""
        fig = self.render_color_mapped(**{
            k: v for k, v in kwargs.items()
            if k in ("figsize", "height_ratios", "expression_fontsize")
        })
        savefig_kwargs = {
            k: v for k, v in kwargs.items()
            if k not in ("figsize", "height_ratios",
                         "expression_fontsize")
        }
        path = Path(path)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", **savefig_kwargs)
        plt.close(fig)
        return path

    # ---------- TikZ-Feynman output ----------

    def to_tikz(self, file: str | Path | None = None) -> str:
        """
        Emit idiomatic TikZ-Feynman code for this diagram.

        Output drops directly into a paper preamble that has
            \\usepackage{tikz}
            \\usepackage{tikz-feynman}
        Each diagram is wrapped in a `feynman` environment; tikz-feynman's
        automatic layout engine handles vertex placement.

        If `file` is given, also writes the TikZ to disk.
        """
        if not self._tikz_template:
            raise NotImplementedError(
                f"No TikZ template for {self.process_name}/{self.channel}"
            )
        code = self._tikz_template
        if file is not None:
            Path(file).write_text(code)
        return code

    # ---------- Repr ----------

    def __repr__(self) -> str:
        return (f"Diagram(process={self.process_name!r}, "
                f"channel={self.channel!r})")


# ---------------------------------------------------------------------------
# TikZ-Feynman templates (idiomatic, no manual coordinates)
# ---------------------------------------------------------------------------

TIKZ_PREAMBLE_NOTE = (
    "% Requires:  \\usepackage{tikz}\n"
    "%            \\usepackage{tikz-feynman}\n"
    "% Compile with lualatex (recommended for tikz-feynman auto-layout).\n"
)


TIKZ_COMPTON_S = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [right=2cm of a] (b);
    \vertex [above left=of a] (i1) {\(e^{-}(p)\)};
    \vertex [below left=of a] (i2) {\(\gamma(k)\)};
    \vertex [above right=of b] (f1) {\(e^{-}(p')\)};
    \vertex [below right=of b] (f2) {\(\gamma(k')\)};
    \diagram* {
      (i1) -- [fermion] (a),
      (i2) -- [photon]  (a),
      (a)  -- [fermion, edge label=\(e^{*}\)] (b),
      (b)  -- [fermion] (f1),
      (b)  -- [photon]  (f2),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_COMPTON_U = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [right=2cm of a] (b);
    \vertex [above left=of a] (i1) {\(e^{-}(p)\)};
    \vertex [below left=of a] (i2) {\(\gamma(k)\)};
    \vertex [above right=of b] (f1) {\(e^{-}(p')\)};
    \vertex [below right=of b] (f2) {\(\gamma(k')\)};
    \diagram* {
      (i1) -- [fermion] (a),
      (a)  -- [fermion, edge label=\(e^{*}\)] (b),
      (b)  -- [fermion] (f1),
      (i2) -- [photon]  (b),
      (a)  -- [photon]  (f2),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_EEMUMU_S = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [above left=of a]  (i1) {\(e^{-}\)};
    \vertex [below left=of a]  (i2) {\(e^{+}\)};
    \vertex [above right=of b] (f1) {\(\mu^{-}\)};
    \vertex [below right=of b] (f2) {\(\mu^{+}\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (i2),
      (a)  -- [photon, edge label=\(\gamma^{*}\)] (b),
      (f2) -- [fermion] (b) -- [fermion] (f1),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_MOLLER_T = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [left=of a]  (i1) {\(e^{-}\)};
    \vertex [right=of a] (f1) {\(e^{-}\)};
    \vertex [left=of b]  (i2) {\(e^{-}\)};
    \vertex [right=of b] (f2) {\(e^{-}\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (f1),
      (a)  -- [photon, edge label=\(\gamma^{*}\)] (b),
      (i2) -- [fermion] (b) -- [fermion] (f2),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_MOLLER_U = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [left=of a]  (i1) {\(e^{-}\)};
    \vertex [right=of a] (f1) {\(e^{-}\)};
    \vertex [left=of b]  (i2) {\(e^{-}\)};
    \vertex [right=of b] (f2) {\(e^{-}\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (f2),
      (a)  -- [photon, edge label=\(\gamma^{*}\)] (b),
      (i2) -- [fermion] (b) -- [fermion] (f1),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_BHABHA_S = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [above left=of a]  (i1) {\(e^{-}\)};
    \vertex [below left=of a]  (i2) {\(e^{+}\)};
    \vertex [above right=of b] (f1) {\(e^{-}\)};
    \vertex [below right=of b] (f2) {\(e^{+}\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (i2),
      (a)  -- [photon, edge label=\(\gamma^{*}\)] (b),
      (f2) -- [fermion] (b) -- [fermion] (f1),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_BHABHA_T = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [below=2cm of a] (b);
    \vertex [left=of a]  (i1) {\(e^{-}\)};
    \vertex [right=of a] (f1) {\(e^{-}\)};
    \vertex [left=of b]  (i2) {\(e^{+}\)};
    \vertex [right=of b] (f2) {\(e^{+}\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (f1),
      (a)  -- [photon, edge label=\(\gamma^{*}\)] (b),
      (f2) -- [fermion] (b) -- [fermion] (i2),
    };
  \end{feynman}
\end{tikzpicture}
"""

TIKZ_MUON_DECAY = TIKZ_PREAMBLE_NOTE + r"""
\begin{tikzpicture}
  \begin{feynman}
    \vertex (a);
    \vertex [right=2cm of a] (b);
    \vertex [left=of a]       (i1) {\(\mu^{-}\)};
    \vertex [above right=of a](f1) {\(\nu_{\mu}\)};
    \vertex [above right=of b](f2) {\(e^{-}\)};
    \vertex [below right=of b](f3) {\(\bar{\nu}_{e}\)};
    \diagram* {
      (i1) -- [fermion] (a) -- [fermion] (f1),
      (a)  -- [boson, edge label=\(W^{-}\)] (b),
      (b)  -- [fermion] (f2),
      (f3) -- [fermion] (b),
    };
  \end{feynman}
\end{tikzpicture}
"""
