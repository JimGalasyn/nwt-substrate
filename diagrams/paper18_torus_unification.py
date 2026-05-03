#!/usr/bin/env python3
"""
Paper 18 §4.4 figure: trefoil + K_7 on the same Heegaard torus.

Two-panel layout illustrating the structural identity
M_Pl^2 / m_e^2 = |A_{K_7}|^{-2}: the (2,3) trefoil BPS background
of Paper 16 and the K_7 Heffter graph of Paper 17 live on the SAME
Heegaard torus T^2 in S^3 / 2I, so they are not two different
substrate objects but two views of one.

Output: paper18_fig_unified_torus.pdf  (vector, paper-grade)
        paper18_fig_unified_torus.png  (raster, for git preview)
"""
from __future__ import annotations
import os
from itertools import combinations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


HERE = os.path.dirname(os.path.abspath(__file__))


PUB_STYLE = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "text.color": "#1a1a1a",
    "axes.labelcolor": "#1a1a1a",
    "axes.edgecolor": "#1a1a1a",
    "xtick.color": "#1a1a1a",
    "ytick.color": "#1a1a1a",
    "grid.color": "#cccccc",
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.10,
}


# Torus geometry
R_MAJOR = 1.5
R_MINOR = 0.55


def torus_xyz(u, v):
    x = (R_MAJOR + R_MINOR * np.cos(v)) * np.cos(u)
    y = (R_MAJOR + R_MINOR * np.cos(v)) * np.sin(u)
    z = R_MINOR * np.sin(v)
    return x, y, z


def draw_torus_surface(ax, alpha=0.15, color="#cfd8e3"):
    """Light translucent torus surface so overlaid lines stay visible."""
    n_u, n_v = 120, 40
    uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, n_u),
                          np.linspace(0, 2 * np.pi, n_v))
    xx, yy, zz = torus_xyz(uu, vv)
    ax.plot_surface(xx, yy, zz, color=color, alpha=alpha,
                     linewidth=0, antialiased=True, shade=True)
    # Light grid lines along the meridians and longitudes
    for u0 in np.linspace(0, 2 * np.pi, 13)[:-1]:
        v_line = np.linspace(0, 2 * np.pi, 80)
        x, y, z = torus_xyz(np.full_like(v_line, u0), v_line)
        ax.plot(x, y, z, color="#a8b3c2", linewidth=0.35, alpha=0.55)
    for v0 in np.linspace(0, 2 * np.pi, 9)[:-1]:
        u_line = np.linspace(0, 2 * np.pi, 200)
        x, y, z = torus_xyz(u_line, np.full_like(u_line, v0))
        ax.plot(x, y, z, color="#a8b3c2", linewidth=0.35, alpha=0.55)


def style_panel(ax, title):
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.2, 2.2); ax.set_zlim(-1.4, 1.4)
    ax.set_box_aspect((1, 1, 0.55))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    # Hide the axis "panels" so the figure reads as a clean 3D scene
    ax.xaxis.set_pane_color((1, 1, 1, 0))
    ax.yaxis.set_pane_color((1, 1, 1, 0))
    ax.zaxis.set_pane_color((1, 1, 1, 0))
    for spine_axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        spine_axis._axinfo["grid"]["linewidth"] = 0
    ax.view_init(elev=30, azim=-58)
    ax.set_title(title, fontsize=12, pad=8)


def fig_unified_torus(out_path_pdf: str, out_path_png: str) -> None:
    plt.rcParams.update(PUB_STYLE)
    fig = plt.figure(figsize=(13.5, 6.3))

    # ----- LEFT PANEL: torus + trefoil (matter / BPS background) -----
    axL = fig.add_subplot(1, 2, 1, projection="3d")
    draw_torus_surface(axL)

    # (2,3) torus knot
    t = np.linspace(0, 2 * np.pi, 800)
    p_w, q_w = 2, 3
    u_k = p_w * t
    v_k = q_w * t
    xk, yk, zk = torus_xyz(u_k, v_k)
    axL.plot(xk, yk, zk, color="#c0392b", linewidth=2.8,
              solid_capstyle="round", zorder=10)

    style_panel(axL,
                 r"(a) BPS background: $(2,3)$ trefoil on $T^2 \subset S^3 / 2I$"
                 "\n"
                 r"matter sector $\;-\;$ Paper 16 vortex $\;-\;$ "
                 r"$(p,q)$ winding gives $(p^2+q^2)$ mass")

    # ----- RIGHT PANEL: torus + K_7 Heffter embedding (gravity / Wilson) -----
    axR = fig.add_subplot(1, 2, 2, projection="3d")
    draw_torus_surface(axR)

    # K_7 Heffter coords on the unit fundamental domain (then * 2pi)
    vert_uv = np.array([[(k % 7) / 7.0, (3 * k % 7) / 7.0] for k in range(7)]) * 2 * np.pi
    vx, vy, vz = torus_xyz(vert_uv[:, 0], vert_uv[:, 1])

    # Edges: shortest geodesic from i to j on the flat torus, then mapped
    edge_color = "#2a4a7a"
    for i, j in combinations(range(7), 2):
        u1, v1 = vert_uv[i]
        u2_raw, v2_raw = vert_uv[j]
        # find closest periodic image of (u2, v2) to (u1, v1)
        best = (u2_raw, v2_raw); best_d = np.inf
        for du in (-2 * np.pi, 0, 2 * np.pi):
            for dv in (-2 * np.pi, 0, 2 * np.pi):
                d = (u2_raw + du - u1) ** 2 + (v2_raw + dv - v1) ** 2
                if d < best_d:
                    best = (u2_raw + du, v2_raw + dv); best_d = d
        u2, v2 = best
        ts = np.linspace(0, 1, 80)
        u_seg = u1 + (u2 - u1) * ts
        v_seg = v1 + (v2 - v1) * ts
        xe, ye, ze = torus_xyz(u_seg, v_seg)
        axR.plot(xe, ye, ze, color=edge_color, linewidth=1.0,
                  alpha=0.78, solid_capstyle="round", zorder=6)

    axR.scatter(vx, vy, vz, s=110, c="#f4ecd9", edgecolors="#1a1a1a",
                 linewidths=1.4, zorder=10, depthshade=False)
    for k, (x, y, z) in enumerate(zip(vx, vy, vz)):
        axR.text(x * 1.18, y * 1.18, z * 1.6, str(k),
                  ha="center", va="center", fontsize=9.5, fontweight="bold",
                  color="#1a1a1a", zorder=11)

    style_panel(axR,
                 r"(b) Wilson graph: $K_7$ Heffter embedding on $T^2 \subset S^3 / 2I$"
                 "\n"
                 r"gravity sector $\;-\;$ Paper 17 amplitude $\;-\;$ "
                 r"$21$ edges $\to A_{K_7} = \mathrm{bracket}(\alpha)\,\alpha^{21/2}$")

    fig.suptitle(
        r"Same Heegaard torus, two sectors:  "
        r"$M_\mathrm{Pl}^2 / m_e^2 = |A_{K_7}|^{-2}$",
        fontsize=14, y=1.00)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path_pdf)
    fig.savefig(out_path_png, dpi=180)
    plt.close(fig)
    print(f"  wrote {out_path_pdf}")
    print(f"  wrote {out_path_png}")


def main() -> None:
    fig_unified_torus(
        os.path.join(HERE, "paper18_fig_unified_torus.pdf"),
        os.path.join(HERE, "paper18_fig_unified_torus.png"),
    )


if __name__ == "__main__":
    main()
