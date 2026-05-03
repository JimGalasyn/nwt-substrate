#!/usr/bin/env python3
"""
Substrate-primitive figures for Papers 18 / 19.

Generates three publication-quality vector PDFs that illustrate the
core substrate objects:

  paper19_fig_trefoil.pdf
      The (2,3) torus knot (trefoil), the canonical NWT carrier knot
      that supports the BPS background of the Sakharov derivation
      and whose Dehn surgery yields the Brieskorn-Poincare sphere.

  paper19_fig_poincare_sphere.pdf
      The Brieskorn-Poincare homology 3-sphere Sigma(2,3,5) = S^3 / 2I,
      drawn as a regular dodecahedron with opposite faces identified
      via a 36-degree twist.

  paper19_fig_K7_heffter.pdf
      The complete graph K_7 on its unique Heffter triangular embedding
      on the torus, drawn in the square fundamental domain.

Style: light publication theme, vector PDF output, axes off where possible.

Usage:
    python3 diagrams/paper19_substrate_figures.py            # all three
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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


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


# ---------------------------------------------------------------------------
# Figure 1: trefoil knot drawn on (and slightly outside) its supporting torus
# ---------------------------------------------------------------------------

def fig_trefoil(out_path: str) -> None:
    plt.rcParams.update(PUB_STYLE)
    fig = plt.figure(figsize=(6.0, 5.5))
    ax = fig.add_subplot(111, projection="3d")

    # Torus wireframe (sparse) -- 3D matplotlib's plot_surface occludes
    # lines regardless of zorder, so we use plot_wireframe with light
    # styling instead.
    R, r = 2.0, 0.7
    theta = np.linspace(0, 2 * np.pi, 60)
    phi = np.linspace(0, 2 * np.pi, 24)
    THETA, PHI = np.meshgrid(theta, phi)
    Xt = (R + r * np.cos(PHI)) * np.cos(THETA)
    Yt = (R + r * np.cos(PHI)) * np.sin(THETA)
    Zt = r * np.sin(PHI)
    ax.plot_wireframe(Xt, Yt, Zt,
                       color="#a3c4dc", linewidth=0.3, alpha=0.6,
                       rstride=3, cstride=3)

    # Trefoil curve.  Plot a single 3D polyline; matplotlib's 3D
    # depth-sorting is fragile when many short segments mix with surfaces,
    # so a single plot call with one strong colour reads more cleanly.
    t = np.linspace(0, 2 * np.pi, 600)
    rr = r + 0.18
    Xk = (R + rr * np.cos(3 * t)) * np.cos(2 * t)
    Yk = (R + rr * np.cos(3 * t)) * np.sin(2 * t)
    Zk = rr * np.sin(3 * t)
    ax.plot(Xk, Yk, Zk, color="#c0392b", linewidth=2.8, zorder=10)
    # Add a short shadow under the curve at z = -r (projected to torus
    # equator) for additional depth cue.
    ax.plot(Xk, Yk, np.full_like(Zk, -r - 0.01),
             color="#888", linewidth=0.9, alpha=0.35, zorder=1)

    # Title and caption
    ax.text2D(0.50, 0.97,
              r"$(p,q) = (2,3)$ torus knot $\equiv\;3_1$ (trefoil)",
              transform=ax.transAxes, ha="center", fontsize=12)
    ax.text2D(0.50, 0.04,
              r"3 crossings, $\gcd(p,q)=1$; Dehn surgery yields "
              r"$\Sigma(2,3,5)$",
              transform=ax.transAxes, ha="center", fontsize=10,
              color="#444444")

    ax.set_axis_off()
    ax.view_init(elev=30, azim=35)
    ax.set_box_aspect((1, 1, 0.45))
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: regular dodecahedron with opposite-face identification arrows
# ---------------------------------------------------------------------------

def _regular_dodecahedron_verts():
    """Return 20 vertices of a unit regular dodecahedron."""
    phi = (1 + np.sqrt(5)) / 2
    iphi = 1 / phi
    cube = [(sx, sy, sz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    extra = []
    for sa in (-iphi, iphi):
        for sb in (-phi, phi):
            extra.append((0, sa, sb))
            extra.append((sb, 0, sa))
            extra.append((sa, sb, 0))
    verts = np.array(cube + extra, dtype=float)
    return verts / np.linalg.norm(verts[0])


def fig_poincare_sphere(out_path: str) -> None:
    """Use scipy ConvexHull to triangulate the dodecahedron surface."""
    from scipy.spatial import ConvexHull

    plt.rcParams.update(PUB_STYLE)
    fig = plt.figure(figsize=(6.0, 5.6))
    ax = fig.add_subplot(111, projection="3d")

    verts = _regular_dodecahedron_verts()
    hull = ConvexHull(verts)
    # Each ConvexHull simplex is a triangle (i, j, k) of vertex indices.
    # The 12 pentagonal faces of the dodecahedron are each triangulated
    # into 3 triangles, giving 36 simplices total -- this is fine for
    # surface rendering.
    tris = [verts[s] for s in hull.simplices]

    # We also identify the 12 pentagon centroids by grouping coplanar
    # triangles via face normals.  This gives us face centres for the
    # twist-arrow annotation.
    normals = []
    for s in hull.simplices:
        a, b, c = verts[s]
        n = np.cross(b - a, c - a)
        n /= np.linalg.norm(n)
        normals.append(n)
    normals = np.array(normals)

    # Cluster simplices by their face-normal direction.  Two normals
    # belong to the same pentagonal face iff they are within ~1 deg.
    face_groups = []   # list of lists of simplex indices
    for i, n in enumerate(normals):
        placed = False
        for group in face_groups:
            n0 = normals[group[0]]
            if np.dot(n, n0) > 0.999:
                group.append(i)
                placed = True
                break
        if not placed:
            face_groups.append([i])

    face_centres = np.array([
        np.mean([verts[s].mean(axis=0)
                 for s in hull.simplices[group]], axis=0)
        for group in face_groups
    ])

    poly = Poly3DCollection(
        tris, facecolor="#f4ecd9", edgecolor="#222222",
        linewidths=0.7, alpha=0.92,
    )
    ax.add_collection3d(poly)

    # Mark the 30 dodecahedron edges (each shared by 2 pentagonal faces).
    # Easier: extract them from the convex hull edge list.
    edge_set = set()
    for s in hull.simplices:
        for i in range(3):
            a, b = sorted((s[i], s[(i + 1) % 3]))
            edge_set.add((a, b))
    # Filter to "true" dodecahedron edges (length close to the minimum
    # edge length, not the diagonals introduced by triangulation).
    if edge_set:
        lengths = np.array([np.linalg.norm(verts[a] - verts[b])
                            for a, b in edge_set])
        edge_len = np.min(lengths)
        true_edges = [(a, b) for (a, b), L in zip(edge_set, lengths)
                      if abs(L - edge_len) < 1e-6]
        for a, b in true_edges:
            ax.plot(*zip(verts[a], verts[b]), color="#222222",
                     linewidth=1.3, zorder=5)

    # Identify the topmost and bottom-most face centres
    f_top = int(np.argmax(face_centres[:, 2]))
    f_bot = int(np.argmin(face_centres[:, 2]))
    centres = face_centres

    for fidx, sense in [(f_top, +1), (f_bot, -1)]:
        c = centres[fidx]
        n = c / np.linalg.norm(c)
        u = np.array([1.0, 0, 0])
        u = u - (u @ n) * n
        u /= np.linalg.norm(u)
        v = np.cross(n, u)
        rho = 0.42
        # Pull the arrow slightly outward along the face normal so it hovers
        # above the face instead of intersecting it.
        offset = 0.18 * n
        ang = np.linspace(0, sense * 0.65, 30)
        pts = (c[:, None] + offset[:, None] +
               rho * (np.cos(ang)[None, :] * u[:, None]
                      + np.sin(ang)[None, :] * v[:, None]))
        ax.plot(pts[0], pts[1], pts[2], color="#c0392b",
                 linewidth=2.4, zorder=20)
        # Arrow head as a marker at the curve tip
        ax.scatter([pts[0, -1]], [pts[1, -1]], [pts[2, -1]],
                    s=70, c="#c0392b", marker=">", zorder=21)
        # Label
        ax.text((c + 1.3 * n)[0], (c + 1.3 * n)[1], (c + 1.3 * n)[2],
                 r"$+\pi/5$ twist", color="#c0392b", fontsize=11,
                 ha="center")

    ax.text2D(0.50, 0.97,
              "Brieskorn–Poincaré sphere "
              r"$\Sigma(2,3,5) = S^3 / 2I$",
              transform=ax.transAxes, ha="center", fontsize=12)
    ax.text2D(0.50, 0.03,
              "regular dodecahedron, opposite faces identified with "
              r"$\pi/5$ twist; "
              r"$|\pi_1|=|2I|=120$, "
              r"$\lambda_1=168=7{\cdot}24$",
              transform=ax.transAxes, ha="center", fontsize=9,
              color="#444444")

    ax.set_axis_off()
    ax.view_init(elev=22, azim=40)
    L = 1.5
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
    ax.set_box_aspect((1, 1, 1))

    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: K_7 on Heffter triangular embedding (square fundamental domain)
# ---------------------------------------------------------------------------

def fig_K7_heffter(out_path: str) -> None:
    plt.rcParams.update(PUB_STYLE)
    fig, ax = plt.subplots(figsize=(6.0, 6.4))

    # Vertex placement: skew sublattice (k/7, 3k mod 7 / 7), gives the
    # Heffter triangular embedding when edges are drawn as shortest
    # toroidal images.
    coords = np.array([[(k % 7) / 7.0,
                         (3 * k % 7) / 7.0] for k in range(7)])

    # Square fundamental domain with identification arrows
    boundary_color = "#666666"
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0],
             color=boundary_color, linewidth=1.6)
    for x_arrow, y_arrow in [
        ((0.50, 0.55), (0, 0)),     # bottom edge double-arrow group
        ((0.50, 0.55), (1, 1)),     # top
    ]:
        ax.annotate("", xy=(0.55, y_arrow[0]), xytext=(0.45, y_arrow[0]),
                     arrowprops=dict(arrowstyle="->", color=boundary_color,
                                       lw=1.5))
    ax.annotate("", xy=(0, 0.55), xytext=(0, 0.45),
                 arrowprops=dict(arrowstyle="->", color=boundary_color, lw=1.5))
    ax.annotate("", xy=(1, 0.55), xytext=(1, 0.45),
                 arrowprops=dict(arrowstyle="->", color=boundary_color, lw=1.5))

    # Edge drawing helper:\ for each pair (i, j), find the shortest toroidal
    # image of j relative to i, draw the segment.  If the segment crosses
    # the boundary, also draw the periodic-image segment so connectivity
    # is visually clear.
    def shortest(i, j):
        p1 = coords[i]
        best = (None, None)
        best_d = np.inf
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                p2 = coords[j] + np.array([dx, dy])
                d = np.linalg.norm(p2 - p1)
                if d < best_d:
                    best = (p1.copy(), p2.copy())
                    best_d = d
        return best

    edge_color = "#1a1a1a"
    for i, j in combinations(range(7), 2):
        p1, p2 = shortest(i, j)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                 color=edge_color, linewidth=0.9, alpha=0.85)
        # If the wrapped p2 is outside the unit square, also draw the
        # mirror segment in the unit square (other half of the wrap).
        out_x = int(p2[0] < 0) - int(p2[0] > 1)
        out_y = int(p2[1] < 0) - int(p2[1] > 1)
        if out_x != 0 or out_y != 0:
            shift = np.array([out_x, out_y])
            ax.plot([p1[0] + shift[0], p2[0] + shift[0]],
                     [p1[1] + shift[1], p2[1] + shift[1]],
                     color=edge_color, linewidth=0.9, alpha=0.85)

    # Vertices on top
    for k, (x, y) in enumerate(coords):
        ax.plot(x, y, "o", markersize=12,
                 markerfacecolor="#f4ecd9",
                 markeredgecolor="#1a1a1a",
                 markeredgewidth=1.5, zorder=5)
        ax.text(x, y, str(k), color="#1a1a1a",
                 ha="center", va="center", fontsize=9,
                 fontweight="bold", zorder=6)

    ax.text(0.5, 1.10,
             r"$K_7$ Heffter triangular embedding on $T^2$",
             ha="center", fontsize=12, transform=ax.transData)
    ax.text(0.5, -0.10,
             r"$V=7$,  $E=21$,  $F=14$,  $\chi=0$,  genus $=1$",
             ha="center", fontsize=10, color="#444444",
             transform=ax.transData)

    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.18, 1.18)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trefoil", action="store_true")
    p.add_argument("--poincare", action="store_true")
    p.add_argument("--k7", action="store_true")
    p.add_argument("--outdir", default=HERE)
    args = p.parse_args()

    do_all = not (args.trefoil or args.poincare or args.k7)

    if args.trefoil or do_all:
        fig_trefoil(os.path.join(args.outdir, "paper19_fig_trefoil.pdf"))
    if args.poincare or do_all:
        fig_poincare_sphere(os.path.join(args.outdir, "paper19_fig_poincare_sphere.pdf"))
    if args.k7 or do_all:
        fig_K7_heffter(os.path.join(args.outdir, "paper19_fig_K7_heffter.pdf"))


if __name__ == "__main__":
    main()
