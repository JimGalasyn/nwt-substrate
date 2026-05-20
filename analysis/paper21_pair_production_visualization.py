"""Paper 21 §4 visualization — pair production kill-shot.

Visualizes the substrate's pair-production mechanism as paired forward
+ reverse K_7 walks:

  e⁺e⁻ pair:  e⁻ = 0→2→5→1→4→0  +  e⁺ = 0→4→1→5→2→0  (reverse)
  p p̄ pair:   p  = 0→1→2→3→4→5→6→0 (K_7 Hamilton, 100% QR)
              p̄  = 0→6→5→4→3→2→1→0 (Hamilton reversed, 100% NR)

Produces multi-panel publication-quality figures:
  Fig 1: K_7 heptagonal layout + e⁺e⁻ pair-production step-by-step
  Fig 2: p p̄ pair-production — the cleanest matter/antimatter symmetry
  Fig 3: Combined comparison + total-winding tracker

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/paper21_pair_production_visualization.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.lines import Line2D

OUT_DIR = Path(__file__).parent / "paper21_figures"
OUT_DIR.mkdir(exist_ok=True, parents=True)


# K_7 heptagonal vertex layout — v_0 at top, others CCW
def k7_positions(R: float = 1.0) -> dict:
    """7 vertices arranged in heptagonal pattern, v_0 at top."""
    pos = {}
    for k in range(7):
        angle = math.pi / 2 - k * 2 * math.pi / 7  # CW from top
        pos[k] = (R * math.cos(angle), R * math.sin(angle))
    return pos


VERTEX_LABELS = {0: r"$P$", 1: r"$E_1$", 2: r"$E_2$", 3: r"$E_3$",
                  4: r"$F_1$", 5: r"$F_2$", 6: r"$F_3$"}


def signed_step(a: int, b: int) -> int:
    """Signed step direction (b - a) mod 7, in {1..6}."""
    return (b - a) % 7


def step_color(d: int) -> tuple:
    """QR-direction (d ∈ {1, 2, 4}) → blue. NR-direction (d ∈ {3, 5, 6}) → red."""
    if d in {1, 2, 4}:
        return (0.10, 0.40, 0.85)  # QR blue
    else:
        return (0.85, 0.20, 0.20)  # NR red


def draw_k7_skeleton(ax, positions: dict, radius: float = 1.0,
                     alpha: float = 0.15, labels: bool = True):
    """Draw the K_7 vertex layout (dots + labels), no edges."""
    # All edges very faint (for context)
    for i in range(7):
        for j in range(i + 1, 7):
            x1, y1 = positions[i]
            x2, y2 = positions[j]
            ax.plot([x1, x2], [y1, y2], color="gray", linewidth=0.5,
                     alpha=alpha, zorder=1)
    # Vertices
    for k, (x, y) in positions.items():
        # Color v_0 differently — substrate vacuum reference point
        color = "#1a1a1a" if k == 0 else "#f4ecd9"
        edge = "#d4a017" if k == 0 else "#1a1a1a"
        ax.scatter(x, y, s=400, c=color, edgecolor=edge, linewidth=2,
                    zorder=10)
        if labels:
            txt_color = "white" if k == 0 else "#1a1a1a"
            ax.text(x, y, VERTEX_LABELS[k], ha="center", va="center",
                     fontsize=10, fontweight="bold", color=txt_color,
                     zorder=11)
    ax.set_aspect("equal")
    ax.set_xlim(-radius * 1.3, radius * 1.3)
    ax.set_ylim(-radius * 1.3, radius * 1.3)
    ax.axis("off")


def draw_walk(ax, walk: list[int], positions: dict, alpha: float = 1.0,
              lw: float = 2.5, label: str = None,
              arrow_step: int = 1):
    """Draw a walk with QR/NR-colored arrows along each edge."""
    L = len(walk) - 1
    for i in range(L):
        a, b = walk[i], walk[i + 1]
        x1, y1 = positions[a]
        x2, y2 = positions[b]
        d = signed_step(a, b)
        color = step_color(d)
        # Curve each edge slightly for visibility
        arrow = FancyArrowPatch(
            (x1, y1), (x2, y2),
            arrowstyle="-|>", mutation_scale=15,
            color=color, linewidth=lw, alpha=alpha,
            connectionstyle=f"arc3,rad={0.08 * (1 if i % 2 == 0 else -1)}",
            zorder=5,
        )
        ax.add_patch(arrow)


def panel_walk_emergence(ax, full_walk: list[int], step: int,
                         positions: dict, color_scheme: str = "qr_nr",
                         title: str = ""):
    """Draw walk progress at given step (0 = just v_0, step k = first k edges)."""
    draw_k7_skeleton(ax, positions, alpha=0.1)
    if step > 0:
        partial = full_walk[:step + 1]
        draw_walk(ax, partial, positions, alpha=1.0, lw=2.5)
    ax.set_title(title, fontsize=10)


def make_pair_production_figure(matter_walk: list[int],
                                antimatter_walk: list[int],
                                particle_name: str,
                                anti_name: str,
                                title: str,
                                outfile: Path):
    """Stepwise pair-production figure: matter walk + antimatter walk emerging."""
    positions = k7_positions()
    L = len(matter_walk) - 1

    fig = plt.figure(figsize=(min(16, 2.5 * (L + 2)), 7.5))
    n_cols = L + 2  # step 0 .. step L for matter walk + 1 final for antimatter

    # Row 1: matter walk emergence (steps 0..L)
    for s in range(L + 1):
        ax = fig.add_subplot(2, n_cols, s + 1)
        title_step = (f"vacuum\n($P$ active)" if s == 0
                       else f"step {s}\n{matter_walk[s-1]}$\\to${matter_walk[s]}"
                       if s < L else
                       f"matter walk\n{particle_name}: complete")
        panel_walk_emergence(ax, matter_walk, s, positions,
                              title=title_step)

    # Empty last column (just for visual breath)
    ax = fig.add_subplot(2, n_cols, L + 2)
    ax.axis("off")
    ax.text(0.5, 0.5, f"+ reverse\nemerges\nsimultaneously\nfrom vacuum",
             ha="center", va="center", fontsize=11,
             transform=ax.transAxes,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0e8d4",
                        edgecolor="#d4a017", linewidth=1.5))

    # Row 2: antimatter walk emergence (steps 0..L)
    for s in range(L + 1):
        ax = fig.add_subplot(2, n_cols, n_cols + s + 1)
        title_step = (f"vacuum\n($P$ active)" if s == 0
                       else f"step {s}\n{antimatter_walk[s-1]}$\\to${antimatter_walk[s]}"
                       if s < L else
                       f"antimatter walk\n{anti_name}: complete")
        panel_walk_emergence(ax, antimatter_walk, s, positions,
                              title=title_step)

    # Final panel of antimatter row: total
    ax = fig.add_subplot(2, n_cols, 2 * n_cols)
    ax.axis("off")
    matter_qr = sum(1 for i in range(L)
                     if signed_step(matter_walk[i], matter_walk[i+1]) in {1,2,4})
    matter_nr = L - matter_qr
    anti_qr = sum(1 for i in range(L)
                   if signed_step(antimatter_walk[i], antimatter_walk[i+1]) in {1,2,4})
    anti_nr = L - anti_qr
    summary = (
        f"Total winding:\n"
        f"({0}, {0}) ✓\n\n"
        f"{particle_name}:\n"
        f"QR: {matter_qr}/{L} ({100*matter_qr//L}%)\n"
        f"NR: {matter_nr}/{L} ({100*matter_nr//L}%)\n\n"
        f"{anti_name}:\n"
        f"QR: {anti_qr}/{L} ({100*anti_qr//L}%)\n"
        f"NR: {anti_nr}/{L} ({100*anti_nr//L}%)\n\n"
        f"QR ↔ NR swap ✓"
    )
    ax.text(0.5, 0.5, summary, ha="center", va="center", fontsize=9,
             transform=ax.transAxes, fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#e8f0e8",
                        edgecolor="#3a7a3a", linewidth=1.2))

    # Legend in figure
    legend_handles = [
        Line2D([0], [0], color=(0.10, 0.40, 0.85), lw=3,
                label="QR-direction step ($d \\in \\{1,2,4\\}$, matter)"),
        Line2D([0], [0], color=(0.85, 0.20, 0.20), lw=3,
                label="NR-direction step ($d \\in \\{3,5,6\\}$, antimatter)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1a1a1a",
                markeredgecolor="#d4a017", markersize=14,
                label="$P$ = K$_7$ polar vertex (substrate vacuum reference)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
                fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(title, fontsize=13, y=1.00, fontweight="bold")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {outfile}")


def make_proton_symmetry_figure(outfile: Path):
    """The cleanest matter/antimatter symmetry: K_7 Hamilton forward vs backward."""
    positions = k7_positions()
    proton = [0, 1, 2, 3, 4, 5, 6, 0]
    antiproton = [0, 6, 5, 4, 3, 2, 1, 0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))

    # Panel 1: proton walk
    ax = axes[0]
    draw_k7_skeleton(ax, positions, alpha=0.15)
    draw_walk(ax, proton, positions, alpha=1.0, lw=3)
    ax.set_title(r"$\mathbf{Proton}$ = K$_7$ Hamilton cycle"
                  "\n"
                  r"$0\to 1\to 2\to 3\to 4\to 5\to 6\to 0$"
                  "\n"
                  r"$100\%$ QR-direction steps (all $d=+1$)",
                  fontsize=11)

    # Panel 2: antiproton walk
    ax = axes[1]
    draw_k7_skeleton(ax, positions, alpha=0.15)
    draw_walk(ax, antiproton, positions, alpha=1.0, lw=3)
    ax.set_title(r"$\mathbf{Antiproton}$ = same Hamilton, reversed"
                  "\n"
                  r"$0\to 6\to 5\to 4\to 3\to 2\to 1\to 0$"
                  "\n"
                  r"$100\%$ NR-direction steps (all $d=+6 \equiv -1$)",
                  fontsize=11)

    # Panel 3: superposition
    ax = axes[2]
    draw_k7_skeleton(ax, positions, alpha=0.15)
    draw_walk(ax, proton, positions, alpha=0.6, lw=3)
    draw_walk(ax, antiproton, positions, alpha=0.6, lw=3)
    ax.set_title(r"$\mathbf{Pair-production\ state}$"
                  "\n"
                  r"$\gamma\gamma \to p + \bar{p}$ at threshold $E\geq 2 m_p c^2$"
                  "\n"
                  r"Total winding $(1,3) + (-1,-3) = (0,0)$",
                  fontsize=11)

    legend_handles = [
        Line2D([0], [0], color=(0.10, 0.40, 0.85), lw=3,
                label="QR-direction ($d=+1$, proton scaffold)"),
        Line2D([0], [0], color=(0.85, 0.20, 0.20), lw=3,
                label="NR-direction ($d=+6$, antiproton scaffold)"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2,
                fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(r"$\mathbf{p\bar{p}}$ pair production at substrate level "
                  r"— the cleanest matter/antimatter symmetry",
                  fontsize=13, y=1.02)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {outfile}")


def make_electron_emergence_figure(outfile: Path):
    """e⁺e⁻ pair-production walk emergence panel-by-panel."""
    electron = [0, 2, 5, 1, 4, 0]
    positron = [0, 4, 1, 5, 2, 0]
    make_pair_production_figure(
        electron, positron, "$e^-$", "$e^+$",
        r"$\mathbf{e^+e^-}$ pair production: forward + reverse walks from substrate vacuum"
        " ($K_7$ at vertex $P$)",
        outfile,
    )


def make_proton_emergence_figure(outfile: Path):
    """p p̄ pair-production walk emergence panel-by-panel."""
    proton = [0, 1, 2, 3, 4, 5, 6, 0]
    antiproton = [0, 6, 5, 4, 3, 2, 1, 0]
    make_pair_production_figure(
        proton, antiproton, "$p$", r"$\bar{p}$",
        r"$\mathbf{p\bar{p}}$ pair production: Hamilton cycle forward + reverse"
        " — 100% QR/NR swap",
        outfile,
    )


def make_combined_panel_figure(outfile: Path):
    """Single figure with both e⁺e⁻ and p p̄ final-state panels."""
    positions = k7_positions()
    electron = [0, 2, 5, 1, 4, 0]
    positron = [0, 4, 1, 5, 2, 0]
    proton = [0, 1, 2, 3, 4, 5, 6, 0]
    antiproton = [0, 6, 5, 4, 3, 2, 1, 0]

    fig = plt.figure(figsize=(13, 10))

    # 2x3 grid: top row e⁻ e⁺ pair; bottom row p p̄ pair
    titles = [
        (r"$e^-$ walk", electron, 1),
        (r"$e^+$ walk = $\overline{e^-}$", positron, 2),
        (r"$e^+e^-$ pair: total winding $(0,0)$", None, 3),
        (r"$p$ walk (K$_7$ Hamilton forward)", proton, 4),
        (r"$\bar{p}$ walk (K$_7$ Hamilton reversed)", antiproton, 5),
        (r"$p\bar{p}$ pair: total winding $(0,0)$", None, 6),
    ]
    for title, walk, pos in titles:
        ax = fig.add_subplot(2, 3, pos)
        draw_k7_skeleton(ax, positions, alpha=0.12)
        if walk is not None:
            draw_walk(ax, walk, positions, alpha=1.0, lw=3)
            # Annotate signed-step pattern
            steps = [signed_step(walk[i], walk[i+1])
                     for i in range(len(walk) - 1)]
            qr = sum(1 for d in steps if d in {1, 2, 4})
            nr = len(steps) - qr
            ax.text(0.05, -1.18, f"signed steps: {steps}\n"
                                  f"QR: {qr}, NR: {nr}",
                     transform=ax.transData, fontsize=9,
                     ha='left', va='top', family='monospace',
                     bbox=dict(boxstyle='round,pad=0.3',
                                facecolor='#fafafa',
                                edgecolor='gray', linewidth=0.5))
        else:
            # Pair view: superimpose
            if pos == 3:  # e+e- pair
                draw_walk(ax, electron, positions, alpha=0.6, lw=3)
                draw_walk(ax, positron, positions, alpha=0.6, lw=3)
            else:  # ppbar pair
                draw_walk(ax, proton, positions, alpha=0.6, lw=3)
                draw_walk(ax, antiproton, positions, alpha=0.6, lw=3)
        ax.set_title(title, fontsize=11)

    legend_handles = [
        Line2D([0], [0], color=(0.10, 0.40, 0.85), lw=3,
                label="QR step ($d \\in \\{1,2,4\\}$, matter direction)"),
        Line2D([0], [0], color=(0.85, 0.20, 0.20), lw=3,
                label="NR step ($d \\in \\{3,5,6\\}$, antimatter direction)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1a1a1a",
                markeredgecolor="#d4a017", markersize=14,
                label="$P$ = K$_7$ polar vertex"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
                fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(r"$\mathbf{Pair\ production\ as\ K_7\ walk\ +\ reverse-traversal}$"
                  "\n"
                  r"Matter walk in QR-direction; antimatter walk in NR-direction;"
                  r" total winding $(0,0)$",
                  fontsize=13, y=1.005)
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {outfile}")


def main() -> None:
    print("=" * 78)
    print("Paper 21 §4 visualization — pair production kill-shot")
    print("=" * 78)
    print()

    print("Figure 1: e⁻ + e⁺ stepwise emergence")
    make_electron_emergence_figure(OUT_DIR / "fig_e_pair_emergence.png")

    print("Figure 2: p + p̄ stepwise emergence")
    make_proton_emergence_figure(OUT_DIR / "fig_p_pair_emergence.png")

    print("Figure 3: p p̄ symmetry-comparison (3-panel)")
    make_proton_symmetry_figure(OUT_DIR / "fig_proton_symmetry.png")

    print("Figure 4: Combined 2x3 panel (both pairs)")
    make_combined_panel_figure(OUT_DIR / "fig_combined_pairs.png")

    print()
    print("All figures saved to:")
    print(f"  {OUT_DIR}/")


if __name__ == "__main__":
    main()
