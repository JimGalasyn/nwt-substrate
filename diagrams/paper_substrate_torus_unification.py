#!/usr/bin/env python3
"""
NWT substrate unification figures.

Produces two paper-grade figures:

  paper_torus_unification.png
      The Heegaard torus of S^3 / 2I with BOTH the (2,3) trefoil
      carrier knot (matter sector, local winding) AND the K_7 Heffter
      triangular embedding (gravity sector, global Wilson graph)
      overlaid on a single stereographically-projected torus.
      One surface, two sectors.

  paper_K7_traversals.png
      The same K_7 graph traversed three different ways for three
      different sectors:

        (A) Eulerian circuit (21 edges) -> alpha^(21/2) gravity amplitude
        (B) Triangular faces (14 of them) -> octonion associators,
            matter-sector 3-body content (Sigma|associator|^2 = 112)
        (C) 4-cycles -> substrate-algebra walks, cycle hierarchy
            21 / 112 / 24 (memory:walk-phase-4g)

The point is to make explicit, in one place, that K_7 isn't just a
graph -- it's a substrate object whose different traversals compute
different physical amplitudes.
"""
from __future__ import annotations
import argparse
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

# Torus geometry: major radius R, minor radius r
R_MAJOR = 1.6
R_MINOR = 0.55


def _torus_xyz(u: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Standard torus parameterisation in R^3 (stereographic-like)."""
    x = (R_MAJOR + R_MINOR * np.cos(v)) * np.cos(u)
    y = (R_MAJOR + R_MINOR * np.cos(v)) * np.sin(u)
    z = R_MINOR * np.sin(v)
    return x, y, z


# ---------------------------------------------------------------------------
# Figure 1: unified torus -- matter (trefoil) + gravity (K_7) on one surface
# ---------------------------------------------------------------------------

def fig_unified_torus(out_path: str) -> None:
    plt.rcParams.update(PUB_STYLE)
    fig = plt.figure(figsize=(8.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")

    # ----- Torus surface (light grey wireframe to keep both knot and graph visible) -----
    n_u, n_v = 80, 30
    uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, n_u),
                          np.linspace(0, 2 * np.pi, n_v))
    xx, yy, zz = _torus_xyz(uu, vv)
    ax.plot_wireframe(xx, yy, zz, rstride=3, cstride=4,
                       color="#bbbbbb", linewidth=0.4, alpha=0.5)

    # ----- (2,3) torus knot (the trefoil, matter sector) -----
    t = np.linspace(0, 2 * np.pi, 600)
    p_w, q_w = 2, 3
    u_k = p_w * t
    v_k = q_w * t
    xk, yk, zk = _torus_xyz(u_k, v_k)
    ax.plot(xk, yk, zk, color="#c0392b", linewidth=2.6,
             label=r"$(2,3)$ trefoil — matter sector $(p^2+q^2)$ winding",
             zorder=5)

    # ----- K_7 Heffter embedding (gravity sector) -----
    # Vertex (u, v) coords: Heffter pattern from paper19 figure
    vert_uv = np.array([[(k % 7) / 7.0,
                          (3 * k % 7) / 7.0] for k in range(7)]) * 2 * np.pi
    vx, vy, vz = _torus_xyz(vert_uv[:, 0], vert_uv[:, 1])

    # Edges: shortest toroidal geodesic between each pair
    edge_color = "#2a4a7a"
    for i, j in combinations(range(7), 2):
        # find shortest periodic-image of j relative to i
        u1, v1 = vert_uv[i]
        u2_raw, v2_raw = vert_uv[j]
        best = (u2_raw, v2_raw)
        best_d = np.inf
        for du in (-2 * np.pi, 0, 2 * np.pi):
            for dv in (-2 * np.pi, 0, 2 * np.pi):
                d = (u2_raw + du - u1) ** 2 + (v2_raw + dv - v1) ** 2
                if d < best_d:
                    best = (u2_raw + du, v2_raw + dv)
                    best_d = d
        u2, v2 = best
        # parametric line in (u, v) space, mapped to torus
        ts = np.linspace(0, 1, 60)
        u_seg = u1 + (u2 - u1) * ts
        v_seg = v1 + (v2 - v1) * ts
        # wrap into [0, 2pi) for parameterisation
        u_seg = u_seg % (2 * np.pi)
        v_seg = v_seg % (2 * np.pi)
        xe, ye, ze = _torus_xyz(u_seg, v_seg)
        ax.plot(xe, ye, ze, color=edge_color, linewidth=0.9, alpha=0.7,
                 zorder=4)

    # Vertices
    ax.scatter(vx, vy, vz, s=80, c="#f4ecd9", edgecolors="#1a1a1a",
                linewidths=1.4, zorder=6,
                label=r"$K_7$ vertices — gravity sector, 21 Wilson edges $\to \alpha^{21/2}$")
    for k, (x, y, z) in enumerate(zip(vx, vy, vz)):
        ax.text(x * 1.18, y * 1.18, z * 1.6, str(k),
                ha="center", va="center", fontsize=9, fontweight="bold",
                color="#1a1a1a", zorder=7)

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_box_aspect((1, 1, 0.4))

    # Camera
    ax.view_init(elev=28, azim=-55)

    ax.set_title(
        "Heegaard torus $T^2 \\subset S^3/2I$ : matter and gravity on one surface\n"
        r"local $(2,3)$ trefoil winding $\;+\;$ global $K_7$ Heffter embedding",
        fontsize=12, pad=18)

    leg = ax.legend(loc="upper left", fontsize=9, framealpha=0.85)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: K_7 traversal modes (Eulerian, triangles, 4-cycles)
# ---------------------------------------------------------------------------

def _heffter_coords() -> np.ndarray:
    """Heffter embedding vertex coordinates on the unit fundamental domain."""
    return np.array([[(k % 7) / 7.0, (3 * k % 7) / 7.0] for k in range(7)])


def _shortest_segment(p1: np.ndarray, p2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return endpoints of the shortest toroidal segment from p1 to a periodic image of p2."""
    best = p2.copy()
    best_d = np.inf
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            cand = p2 + np.array([dx, dy])
            d = np.linalg.norm(cand - p1)
            if d < best_d:
                best = cand
                best_d = d
    return p1.copy(), best


def _draw_torus_panel_background(ax, coords: np.ndarray) -> None:
    """Draw fundamental-domain box, K_7 vertices, and identification marks."""
    boundary_color = "#666666"
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0],
             color=boundary_color, linewidth=1.4)
    # arrow markers on boundaries (top<->bottom, left<->right) to indicate identification
    for y in (0, 1):
        ax.annotate("", xy=(0.55, y), xytext=(0.45, y),
                     arrowprops=dict(arrowstyle="->", color=boundary_color, lw=1.2))
    for x in (0, 1):
        ax.annotate("", xy=(x, 0.55), xytext=(x, 0.45),
                     arrowprops=dict(arrowstyle="->", color=boundary_color, lw=1.2))


def _draw_vertices(ax, coords: np.ndarray) -> None:
    for k, (x, y) in enumerate(coords):
        ax.plot(x, y, "o", markersize=11,
                 markerfacecolor="#f4ecd9",
                 markeredgecolor="#1a1a1a",
                 markeredgewidth=1.2, zorder=5)
        ax.text(x, y, str(k), color="#1a1a1a",
                 ha="center", va="center", fontsize=8.5,
                 fontweight="bold", zorder=6)


def _draw_edge(ax, coords, i, j, color, linewidth=0.9, alpha=0.85,
                 zorder=2):
    p1, p2 = _shortest_segment(coords[i], coords[j])
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
             color=color, linewidth=linewidth, alpha=alpha, zorder=zorder)
    out_x = int(p2[0] < 0) - int(p2[0] > 1)
    out_y = int(p2[1] < 0) - int(p2[1] > 1)
    if out_x != 0 or out_y != 0:
        shift = np.array([out_x, out_y])
        ax.plot([p1[0] + shift[0], p2[0] + shift[0]],
                 [p1[1] + shift[1], p2[1] + shift[1]],
                 color=color, linewidth=linewidth, alpha=alpha, zorder=zorder)


def _eulerian_circuit_K7() -> list[tuple[int, int]]:
    """
    Return one Eulerian circuit on K_7 as a list of 21 edges.
    K_7 has all-even degree (6), so Eulerian circuits exist.  Hierholzer's
    algorithm gives one constructively.
    """
    from collections import defaultdict
    adj: dict[int, list[int]] = defaultdict(list)
    for i, j in combinations(range(7), 2):
        adj[i].append(j)
        adj[j].append(i)
    # Hierholzer
    stack = [0]
    circuit: list[int] = []
    while stack:
        v = stack[-1]
        if adj[v]:
            u = adj[v].pop()
            adj[u].remove(v)
            stack.append(u)
        else:
            circuit.append(stack.pop())
    circuit.reverse()
    edges = [(circuit[i], circuit[i + 1]) for i in range(len(circuit) - 1)]
    assert len(edges) == 21
    return edges


def _heffter_triangles() -> list[tuple[int, int, int]]:
    """
    Return the 14 triangular faces of the K_7 Heffter triangular embedding
    on the torus.  K_7 has C(7,3) = 35 possible triangles; only 14 are
    the embedded faces (Euler relation: V - E + F = 7 - 21 + 14 = 0,
    genus 1).

    Picked geometrically: in the toroidal embedding (u, v) = (k/7,
    3k/7 mod 7), the embedded face triangles have the smallest area
    when each side is taken as the shortest periodic image.  The 14
    smallest-area unordered triples are the embedded faces.
    """
    coords = _heffter_coords()
    tri_areas = []
    for i, j, k in combinations(range(7), 3):
        best = np.inf
        for dx_j in (-1, 0, 1):
            for dy_j in (-1, 0, 1):
                for dx_k in (-1, 0, 1):
                    for dy_k in (-1, 0, 1):
                        pj = coords[j] + np.array([dx_j, dy_j])
                        pk = coords[k] + np.array([dx_k, dy_k])
                        v1 = pj - coords[i]
                        v2 = pk - coords[i]
                        a = 0.5 * abs(v1[0] * v2[1] - v1[1] * v2[0])
                        if a < best:
                            best = a
        tri_areas.append((best, (i, j, k)))
    tri_areas.sort()
    return [t[1] for t in tri_areas[:14]]


def _draw_triangle_face(ax, coords, i, j, k, fill_color, edge_color):
    """Fill the triangle (i, j, k) at its smallest periodic image."""
    # Find which periodic-image triple gives the smallest triangle
    best = None
    best_a = np.inf
    for dx_j in (-1, 0, 1):
        for dy_j in (-1, 0, 1):
            for dx_k in (-1, 0, 1):
                for dy_k in (-1, 0, 1):
                    pj = coords[j] + np.array([dx_j, dy_j])
                    pk = coords[k] + np.array([dx_k, dy_k])
                    v1 = pj - coords[i]; v2 = pk - coords[i]
                    a = 0.5 * abs(v1[0] * v2[1] - v1[1] * v2[0])
                    if a < best_a:
                        best_a = a
                        best = (coords[i], pj, pk)
    p_i, p_j, p_k = best
    poly = plt.Polygon([p_i, p_j, p_k], facecolor=fill_color,
                        edgecolor=edge_color, alpha=0.55, linewidth=0.6,
                        zorder=2)
    ax.add_patch(poly)


def fig_K7_traversals(out_path: str) -> None:
    plt.rcParams.update(PUB_STYLE)
    fig, axes = plt.subplots(1, 3, figsize=(18.0, 7.5))

    coords = _heffter_coords()

    # ----- Panel A: Eulerian circuit (gravity sector) -----
    ax = axes[0]
    _draw_torus_panel_background(ax, coords)
    eul = _eulerian_circuit_K7()
    cmap = plt.get_cmap("plasma")
    for step, (i, j) in enumerate(eul):
        col = cmap(step / 20.0)
        _draw_edge(ax, coords, i, j, color=col, linewidth=2.0, alpha=0.95,
                    zorder=3)
    _draw_vertices(ax, coords)
    ax.set_title("(A) Eulerian circuit  →  gravity\n"
                  r"21 edges as one closed walk: $\alpha^{21/2}$ amplitude",
                  fontsize=13, pad=10)
    ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.18, 1.18)
    ax.set_aspect("equal"); ax.axis("off")

    # ----- Panel B: triangular faces (matter sector via associators) -----
    ax = axes[1]
    _draw_torus_panel_background(ax, coords)
    tris = _heffter_triangles()
    palette = ["#c0392b", "#27ae60", "#2980b9", "#f39c12",
                "#9b59b6", "#16a085", "#e67e22"]
    for n, (i, j, k) in enumerate(tris):
        _draw_triangle_face(ax, coords, i, j, k,
                              fill_color=palette[n % len(palette)],
                              edge_color="#1a1a1a")
    # All 21 edges in light black on top
    for i, j in combinations(range(7), 2):
        _draw_edge(ax, coords, i, j, color="#1a1a1a", linewidth=0.6,
                    alpha=0.7, zorder=2)
    _draw_vertices(ax, coords)
    ax.set_title("(B) Triangular faces  →  matter associators\n"
                  r"14 triangles, octonion 3-body $\Sigma\,|[\,,\,,]\,|^{2}=112$",
                  fontsize=13, pad=10)
    ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.18, 1.18)
    ax.set_aspect("equal"); ax.axis("off")

    # ----- Panel C: 4-cycles (substrate algebra walks) -----
    ax = axes[2]
    _draw_torus_panel_background(ax, coords)
    # Light all 21 edges
    for i, j in combinations(range(7), 2):
        _draw_edge(ax, coords, i, j, color="#bbbbbb", linewidth=0.7,
                    alpha=0.8, zorder=1)
    # Highlight one 4-cycle (e.g. 0-1-3-5-0) and another (0-2-4-6-0)
    cycles = [
        ([0, 1, 3, 5], "#c0392b"),
        ([0, 2, 4, 6], "#2980b9"),
    ]
    for cyc, col in cycles:
        for a, b in zip(cyc, cyc[1:] + [cyc[0]]):
            _draw_edge(ax, coords, a, b, color=col, linewidth=2.4,
                        alpha=0.9, zorder=3)
    _draw_vertices(ax, coords)
    ax.set_title("(C) 4-cycles  →  substrate-algebra walks\n"
                  r"cycle hierarchy 21 / 112 / 24 (walk-phase 4f, 4g)",
                  fontsize=13, pad=10)
    ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.18, 1.18)
    ax.set_aspect("equal"); ax.axis("off")

    fig.suptitle(r"$K_{7}$ traversal modes: same graph, different sectors",
                  fontsize=15, y=0.99)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--unified", action="store_true",
                    help="just the unified-torus figure")
    p.add_argument("--traversals", action="store_true",
                    help="just the K_7 traversal-modes figure")
    args = p.parse_args()
    do_all = not (args.unified or args.traversals)

    if args.unified or do_all:
        fig_unified_torus(os.path.join(HERE, "paper_torus_unification.png"))
    if args.traversals or do_all:
        fig_K7_traversals(os.path.join(HERE, "paper_K7_traversals.png"))


if __name__ == "__main__":
    main()
