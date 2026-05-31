#!/usr/bin/env python3
"""
README hero figure: one QED amplitude, two views.

Left  — the Compton process gamma e- -> gamma e- (s-channel), drawn by the
        QED shim's Feynman renderer with each line colored to match its
        chunk of the substrate amplitude i M printed underneath.
Right — the substrate those gamma-matrices descend from: the complete graph
        K_7 in its Heffter triangular embedding on the torus (7 vertices,
        21 edges = the so(7) generators), with the Cl(0,7) -> Cl(1,3)
        reduction spelled out.

Everything is produced by the library's own `diagrams` machinery:
  * nwt_substrate.qed.compton.diagrams.s_channel.render_color_mapped()
  * nwt_substrate.diagrams.draw_torus / draw_K7_heffter

Usage:
    python3 diagrams/readme_hero.py [OUTPUT.png]

Default output: docs/assets/readme_hero.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nwt_substrate.qed as qed
from nwt_substrate.diagrams import draw_K7_heffter, draw_torus


def build_hero(figsize=(15.0, 6.6)):
    """Build and return the two-view hero Figure."""
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.30, 1.0],
        height_ratios=[7.0, 1.35],
        wspace=0.02, hspace=0.04,
        left=0.012, right=0.988, top=0.88, bottom=0.02,
    )
    ax_diag = fig.add_subplot(gs[0, 0])
    ax_expr = fig.add_subplot(gs[1, 0])
    ax_torus = fig.add_subplot(gs[0, 1], projection="3d")
    ax_text = fig.add_subplot(gs[1, 1])

    # ---- Left: color-mapped Compton (algebra <-> picture) --------------
    qed.compton.diagrams.s_channel.render_color_mapped(
        ax_diagram=ax_diag, ax_expression=ax_expr, expression_fontsize=13,
    )
    ax_diag.set_xlim(0.0, 1.0)
    ax_diag.set_ylim(0.12, 0.92)          # crop dead space above/below
    ax_diag.set_title(
        r"Compton  $\gamma\,e^{-}\!\to\gamma\,e^{-}$  (s-channel)",
        fontsize=13.5, pad=4,
    )

    # ---- Right: the K_7 / Cl(0,7) substrate ----------------------------
    draw_torus(ax_torus, alpha=0.10)
    draw_K7_heffter(
        ax_torus, edge_linewidth=1.35, vertex_size=95, label_fontsize=9.0,
    )
    ax_torus.set_axis_off()
    ax_torus.set_box_aspect((1, 1, 0.46))
    ax_torus.view_init(elev=32, azim=-58)
    # Torus extent is +/-(R+r)=+/-2.05 in x,y and +/-r=0.55 in z; vertex
    # labels are pushed out to ~2.4 (xy) and ~0.9 (z).  Crop close to that
    # so the surface fills the panel instead of floating.
    ax_torus.set_xlim(-2.15, 2.15)
    ax_torus.set_ylim(-2.15, 2.15)
    ax_torus.set_zlim(-0.95, 0.95)
    ax_torus.set_title(
        r"$K_7$ on the torus  $\Rightarrow$  Cl(0,7)",
        fontsize=13.5, pad=-4,
    )

    # ---- Right caption: the Cl(0,7) -> Cl(1,3) descent -----------------
    ax_text.axis("off")
    ax_text.text(
        0.5, 0.95,
        "7 imaginary generators  =  |V($K_7$)|\n"
        r"pick 4  $\to$  Dirac $\gamma^{\mu}$ (8$\times$8)"
        "      remaining 3  $\\to$  internal SU(2)\n"
        "21 edges of $K_7$  =  so(7) generators",
        transform=ax_text.transAxes, ha="center", va="top",
        fontsize=11.0, color="#333333", linespacing=1.45,
    )

    fig.suptitle(
        r"One amplitude, two views — the QED process (left) and the "
        r"Cl(0,7)/$K_7$ substrate its $\gamma$-matrices come from (right)",
        fontsize=15, y=0.985,
    )
    return fig


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        Path(__file__).resolve().parent.parent / "docs" / "assets"
        / "readme_hero.png"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig = build_hero()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
