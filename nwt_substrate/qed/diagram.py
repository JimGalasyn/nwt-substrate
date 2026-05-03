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
        LaTeX string of the substrate-algebra amplitude i M.
    feynman_rules : dict
        Metadata: which propagators / vertices appear in this diagram.
    _render_fn : callable
        Internal: matplotlib renderer that draws the diagram on an Axes.
    _tikz_template : str
        Internal: TikZ-Feynman LaTeX template (idiomatic, lets tikz-feynman
        do auto-layout).
    """

    process_name: str
    channel: str
    expression: str
    feynman_rules: dict = field(default_factory=dict)
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
