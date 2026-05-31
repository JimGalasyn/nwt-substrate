"""
Substrate diagram primitives and figure factories.

A programmatic interface to the canonical NWT substrate diagrams that
have appeared in the gravitational-trilogy papers (15/16/17/18) and
the substrate-monism paper (19).  The module exposes:

  Geometry primitives
    torus_xyz(u, v, R, r)        torus parameterisation
    draw_torus(ax, ...)           render the torus surface
    draw_torus_knot(ax, p, q)    (p,q) torus knot on the surface
    draw_K7_heffter(ax, ...)      K_7 vertices + 21 edges on the torus

  Combinatorial helpers
    K7_eulerian_circuit()         one Eulerian walk through all 21 edges
    K7_heffter_triangles()        14 triangular faces of the embedding

  Figure factories (return matplotlib.figure.Figure)
    figure_torus_knot(p, q, ...)
    figure_K7_on_torus(...)
    figure_unified(...)            both on the same torus, side-by-side
    figure_K7_traversals(mode=...) one of {"eulerian", "triangles", "cycles"}
    figure_K7_all_traversals(...)  3-panel showing all three modes
    figure_paper18_unified(...)    the Paper 18 §4.4 figure

Quick start:

    >>> import nwt_substrate as nwt
    >>> fig = nwt.diagrams.figure_unified()
    >>> fig.savefig("substrate.pdf")

    >>> fig = nwt.diagrams.figure_torus_knot(p=2, q=5)   # cinquefoil
    >>> fig.savefig("cinquefoil.png", dpi=200)

    >>> # raw helpers for custom figures
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure(); ax = fig.add_subplot(111, projection="3d")
    >>> nwt.diagrams.draw_torus(ax)
    >>> nwt.diagrams.draw_torus_knot(ax, p=3, q=2, color="purple")

All factories accept an optional `save_to=path` argument that writes
the figure to disk before returning it.
"""
from __future__ import annotations

from itertools import combinations
from typing import Iterable, Literal, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Default torus geometry (consistent with paper figures)
# ---------------------------------------------------------------------------

DEFAULT_R_MAJOR = 1.5
DEFAULT_R_MINOR = 0.55

# K_7 Heffter embedding: vertex (u, v) coordinates on the unit torus
# (rotation system rho_v cyclic order at v gives a triangular embedding;
# (k/7, 3k/7 mod 7 / 7) is the canonical placement).
HEFFTER_VERT_UV = np.array(
    [[(k % 7) / 7.0, (3 * k % 7) / 7.0] for k in range(7)],
    dtype=float,
)


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------

def torus_xyz(u, v, R: float = DEFAULT_R_MAJOR, r: float = DEFAULT_R_MINOR):
    """
    Parameterise a torus T^2 in R^3.

    Parameters
    ----------
    u, v : float or array-like
        Toroidal and poloidal angles in [0, 2 pi).
    R, r : float
        Major and minor radii.

    Returns
    -------
    (x, y, z) : tuple of arrays
        Cartesian coordinates of the torus surface.
    """
    x = (R + r * np.cos(v)) * np.cos(u)
    y = (R + r * np.cos(v)) * np.sin(u)
    z = r * np.sin(v)
    return x, y, z


def torus_knot_curve(p: int, q: int, *,
                     R: float = DEFAULT_R_MAJOR, r: float = DEFAULT_R_MINOR,
                     n_points: int = 800, closed: bool = False) -> np.ndarray:
    """
    Return the (p, q) torus-knot curve as an array of 3-D points.

    The knot winds ``p`` times toroidally and ``q`` times poloidally on the
    torus of radii (R, r): the curve is ``torus_xyz(p t, q t)`` for
    t in [0, 2 pi).  (p, q) = (2, 3) is the trefoil, (2, 5) the cinquefoil,
    (1, q)/(p, 1) unknots.

    Parameters
    ----------
    p, q : int
        Toroidal and poloidal winding numbers.
    R, r : float
        Major and minor radii of the embedding torus.
    n_points : int
        Number of samples along the curve.
    closed : bool
        If True, repeat the first point at the end so the returned polyline
        is explicitly closed (useful for plotting); if False (default), the
        ``n_points`` samples cover [0, 2 pi) without duplication (useful as a
        point cloud, e.g. for nearest-distance queries).

    Returns
    -------
    np.ndarray
        Array of shape ``(n_points, 3)`` (or ``(n_points + 1, 3)`` if
        ``closed``).
    """
    endpoint = bool(closed)
    t = np.linspace(0.0, 2 * np.pi, n_points, endpoint=False)
    if endpoint:
        t = np.append(t, 0.0)
    x, y, z = torus_xyz(p * t, q * t, R, r)
    return np.stack([x, y, z], axis=1)


def hopf_link_curves(*, R: float = DEFAULT_R_MAJOR, r: float = DEFAULT_R_MINOR,
                     n_points: int = 800) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the two component curves of a Hopf link (linking number 1).

    The Hopf link is the meson carrier (Paper 8): two unknots that each pass
    through the other's disc exactly once.  Component A is the core circle of
    radius R in the z = 0 plane; component B is an identical circle rotated
    90 degrees about the x-axis and shifted by R along x, so the two rings
    interlock once.

    Returns
    -------
    (A, B) : tuple of np.ndarray
        Two arrays of shape ``(n_points, 3)``.
    """
    t = np.linspace(0.0, 2 * np.pi, n_points, endpoint=False)
    A = np.stack([R * np.cos(t), R * np.sin(t), np.zeros_like(t)], axis=1)
    B = np.stack([R * np.cos(t) + R, np.zeros_like(t), R * np.sin(t)], axis=1)
    return A, B


def draw_torus(ax, *, R: float = DEFAULT_R_MAJOR, r: float = DEFAULT_R_MINOR,
                color: str = "#cfd8e3", alpha: float = 0.15,
                grid_color: str = "#a8b3c2", grid_alpha: float = 0.55,
                n_u: int = 120, n_v: int = 40,
                n_meridians: int = 13, n_longitudes: int = 9) -> None:
    """
    Render a translucent torus surface with light meridian/longitude grid.

    Designed to keep overlaid knots/graphs visually clean.
    """
    uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, n_u),
                          np.linspace(0, 2 * np.pi, n_v))
    xx, yy, zz = torus_xyz(uu, vv, R, r)
    ax.plot_surface(xx, yy, zz, color=color, alpha=alpha,
                     linewidth=0, antialiased=True, shade=True)
    for u0 in np.linspace(0, 2 * np.pi, n_meridians)[:-1]:
        v_line = np.linspace(0, 2 * np.pi, 80)
        x, y, z = torus_xyz(np.full_like(v_line, u0), v_line, R, r)
        ax.plot(x, y, z, color=grid_color, linewidth=0.35, alpha=grid_alpha)
    for v0 in np.linspace(0, 2 * np.pi, n_longitudes)[:-1]:
        u_line = np.linspace(0, 2 * np.pi, 200)
        x, y, z = torus_xyz(u_line, np.full_like(u_line, v0), R, r)
        ax.plot(x, y, z, color=grid_color, linewidth=0.35, alpha=grid_alpha)


def draw_torus_knot(ax, p: int, q: int, *,
                     R: float = DEFAULT_R_MAJOR, r: float = DEFAULT_R_MINOR,
                     color: str = "#c0392b", linewidth: float = 2.6,
                     n_points: int = 800) -> None:
    """
    Draw a (p, q) torus knot on the torus surface.

    Parameters
    ----------
    p, q : int
        Toroidal and poloidal winding numbers.  (p, q) = (2, 3) is the
        trefoil; (2, 5) is the cinquefoil; (1, q) and (p, 1) are unknots.
    """
    curve = torus_knot_curve(p, q, R=R, r=r, n_points=n_points, closed=True)
    ax.plot(curve[:, 0], curve[:, 1], curve[:, 2], color=color,
             linewidth=linewidth, solid_capstyle="round", zorder=10)


def draw_K7_heffter(ax, *,
                     R: float = DEFAULT_R_MAJOR, r: float = DEFAULT_R_MINOR,
                     edge_color: str = "#2a4a7a", edge_linewidth: float = 1.0,
                     edge_alpha: float = 0.78,
                     vertex_face: str = "#f4ecd9",
                     vertex_edge: str = "#1a1a1a",
                     vertex_size: float = 110.0,
                     show_labels: bool = True,
                     label_fontsize: float = 9.5) -> None:
    """
    Draw the complete graph K_7 in its unique Heffter triangular embedding
    on the torus.

    Each of the 21 edges is rendered as the shortest geodesic segment
    (interpolated linearly in the (u, v) flat-torus coordinates and
    pushed onto the surface).  Vertices are labelled 0..6.
    """
    vert_uv = HEFFTER_VERT_UV * 2 * np.pi
    vx, vy, vz = torus_xyz(vert_uv[:, 0], vert_uv[:, 1], R, r)

    for i, j in combinations(range(7), 2):
        u1, v1 = vert_uv[i]
        u2_raw, v2_raw = vert_uv[j]
        # find shortest periodic image of vertex j relative to i
        best_u, best_v = u2_raw, v2_raw
        best_d = np.inf
        for du in (-2 * np.pi, 0, 2 * np.pi):
            for dv in (-2 * np.pi, 0, 2 * np.pi):
                d = (u2_raw + du - u1) ** 2 + (v2_raw + dv - v1) ** 2
                if d < best_d:
                    best_u, best_v = u2_raw + du, v2_raw + dv
                    best_d = d
        ts = np.linspace(0, 1, 80)
        u_seg = u1 + (best_u - u1) * ts
        v_seg = v1 + (best_v - v1) * ts
        xe, ye, ze = torus_xyz(u_seg, v_seg, R, r)
        ax.plot(xe, ye, ze, color=edge_color, linewidth=edge_linewidth,
                 alpha=edge_alpha, solid_capstyle="round", zorder=6)

    ax.scatter(vx, vy, vz, s=vertex_size, c=vertex_face,
                edgecolors=vertex_edge, linewidths=1.4, zorder=10,
                depthshade=False)
    if show_labels:
        for k, (x, y, z) in enumerate(zip(vx, vy, vz)):
            ax.text(x * 1.18, y * 1.18, z * 1.6, str(k),
                     ha="center", va="center",
                     fontsize=label_fontsize, fontweight="bold",
                     color="#1a1a1a", zorder=11)


# ---------------------------------------------------------------------------
# K_7 combinatorial helpers
# ---------------------------------------------------------------------------

def K7_eulerian_circuit(start: int = 0) -> list[tuple[int, int]]:
    """
    Return one Eulerian circuit on K_7 as a list of 21 directed edges.

    K_7 has all-even vertex degree (6), so Eulerian circuits exist
    (constructed via Hierholzer's algorithm).  The traversal supplies
    the carrier of the alpha^(21/2) gravity-sector amplitude in
    Paper 17 / Paper 18 §4.4.
    """
    from collections import defaultdict
    adj: dict[int, list[int]] = defaultdict(list)
    for i, j in combinations(range(7), 2):
        adj[i].append(j)
        adj[j].append(i)
    stack = [start]
    walk: list[int] = []
    while stack:
        v = stack[-1]
        if adj[v]:
            u = adj[v].pop()
            adj[u].remove(v)
            stack.append(u)
        else:
            walk.append(stack.pop())
    walk.reverse()
    edges = [(walk[i], walk[i + 1]) for i in range(len(walk) - 1)]
    if len(edges) != 21:
        raise RuntimeError(f"Eulerian circuit produced {len(edges)} edges, "
                            "expected 21")
    return edges


def K7_heffter_triangles() -> list[tuple[int, int, int]]:
    """
    Return the 14 triangular faces of the K_7 Heffter triangular embedding
    on the torus.  Selected geometrically as the 14 unordered triples of
    K_7 vertices with the smallest area in the toroidal embedding.

    Euler relation V - E + F = 0 (genus 1) with V = 7, E = 21 gives
    F = 14 triangular faces; these correspond to the octonion-associator
    structure  Sum_face |[e_i, e_j, e_k]|^2 = 112  in Paper 17 §6.11.
    """
    coords = HEFFTER_VERT_UV
    tri_areas: list[tuple[float, tuple[int, int, int]]] = []
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


# ---------------------------------------------------------------------------
# Figure factories
# ---------------------------------------------------------------------------

_PUB_STYLE = {
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
}


def _style_3d_panel(ax, title: str = "") -> None:
    ax.set_xlim(-2.2, 2.2); ax.set_ylim(-2.2, 2.2); ax.set_zlim(-1.4, 1.4)
    ax.set_box_aspect((1, 1, 0.55))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.xaxis.set_pane_color((1, 1, 1, 0))
    ax.yaxis.set_pane_color((1, 1, 1, 0))
    ax.zaxis.set_pane_color((1, 1, 1, 0))
    ax.view_init(elev=30, azim=-58)
    if title:
        ax.set_title(title, fontsize=12, pad=8)


def _save(fig, path: Optional[str]):
    if path is not None:
        fig.savefig(path, bbox_inches="tight")
    return fig


def figure_torus_knot(p: int = 2, q: int = 3, *,
                       title: Optional[str] = None,
                       save_to: Optional[str] = None,
                       figsize: tuple[float, float] = (6.5, 6.0)):
    """
    Single-panel figure: a (p, q) torus knot on the standard torus.

    Defaults to the trefoil (p=2, q=3).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with plt.rc_context(_PUB_STYLE):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="3d")
        draw_torus(ax)
        draw_torus_knot(ax, p, q)
        _style_3d_panel(ax, title or fr"$({p},{q})$ torus knot on $T^{{2}}$")
        fig.tight_layout()
    return _save(fig, save_to)


def figure_K7_on_torus(*, title: Optional[str] = None,
                         save_to: Optional[str] = None,
                         figsize: tuple[float, float] = (6.5, 6.0)):
    """Single-panel figure: K_7 Heffter embedding on the torus."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with plt.rc_context(_PUB_STYLE):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="3d")
        draw_torus(ax)
        draw_K7_heffter(ax)
        _style_3d_panel(ax, title or r"$K_{7}$ Heffter embedding on $T^{2}$")
        fig.tight_layout()
    return _save(fig, save_to)


def figure_unified(*, p: int = 2, q: int = 3,
                    layout: Literal["side-by-side", "single"] = "side-by-side",
                    title: Optional[str] = None,
                    save_to: Optional[str] = None,
                    figsize: Optional[tuple[float, float]] = None):
    """
    Matter and gravity sectors on the same Heegaard torus.

    layout="side-by-side"  (default):  two panels, trefoil and K_7 each
                                       on their own torus
    layout="single":                   one panel with both objects
                                       overlaid on a single torus
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with plt.rc_context(_PUB_STYLE):
        if layout == "single":
            fig = plt.figure(figsize=figsize or (8.0, 7.0))
            ax = fig.add_subplot(111, projection="3d")
            draw_torus(ax)
            draw_torus_knot(ax, p, q)
            draw_K7_heffter(ax)
            _style_3d_panel(ax,
                title or fr"Heegaard torus: $({p},{q})$ knot $+$ $K_{{7}}$")
        else:
            fig = plt.figure(figsize=figsize or (13.5, 6.3))
            axL = fig.add_subplot(1, 2, 1, projection="3d")
            draw_torus(axL); draw_torus_knot(axL, p, q)
            _style_3d_panel(axL, fr"(a) $({p},{q})$ trefoil on $T^{{2}}$")
            axR = fig.add_subplot(1, 2, 2, projection="3d")
            draw_torus(axR); draw_K7_heffter(axR)
            _style_3d_panel(axR, r"(b) $K_{7}$ Heffter embedding on $T^{2}$")
            fig.suptitle(title
                or r"Same Heegaard torus, two sectors:  "
                   r"$M_\mathrm{Pl}^{2}/m_e^{2} = |A_{K_{7}}|^{-2}$",
                fontsize=14, y=1.00)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save(fig, save_to)


def figure_K7_traversals(mode: Literal["eulerian", "triangles", "cycles"]
                          = "eulerian", *,
                          title: Optional[str] = None,
                          save_to: Optional[str] = None,
                          figsize: tuple[float, float] = (6.5, 6.0)):
    """
    Single-panel figure: K_7 traversed in one of three sector-specific ways.

    mode="eulerian"   (gravity) full Eulerian circuit, 21 edges colour-graded
    mode="triangles"  (matter)  14 triangular faces filled
    mode="cycles"     (algebra) two highlighted 4-cycles
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coords = HEFFTER_VERT_UV
    with plt.rc_context(_PUB_STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        # fundamental-domain box + identification arrows
        ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color="#666666", lw=1.4)
        for y in (0, 1):
            ax.annotate("", xy=(0.55, y), xytext=(0.45, y),
                         arrowprops=dict(arrowstyle="->", color="#666666",
                                          lw=1.2))
        for x in (0, 1):
            ax.annotate("", xy=(x, 0.55), xytext=(x, 0.45),
                         arrowprops=dict(arrowstyle="->", color="#666666",
                                          lw=1.2))

        def _draw_edge(i, j, color, lw=0.9, alpha=0.85, zorder=2):
            p1 = coords[i]; p2 = coords[j]
            best = p2.copy(); best_d = np.inf
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    cand = p2 + np.array([dx, dy])
                    d = np.linalg.norm(cand - p1)
                    if d < best_d:
                        best = cand; best_d = d
            ax.plot([p1[0], best[0]], [p1[1], best[1]],
                     color=color, lw=lw, alpha=alpha, zorder=zorder)
            out_x = int(best[0] < 0) - int(best[0] > 1)
            out_y = int(best[1] < 0) - int(best[1] > 1)
            if out_x or out_y:
                shift = np.array([out_x, out_y])
                ax.plot([p1[0] + shift[0], best[0] + shift[0]],
                         [p1[1] + shift[1], best[1] + shift[1]],
                         color=color, lw=lw, alpha=alpha, zorder=zorder)

        if mode == "eulerian":
            cmap = plt.get_cmap("plasma")
            for step, (i, j) in enumerate(K7_eulerian_circuit()):
                _draw_edge(i, j, color=cmap(step / 20.0), lw=2.0,
                            alpha=0.95, zorder=3)
            default_title = ("Eulerian circuit on $K_{7}$"
                              "  →  21 edges, $\\alpha^{21/2}$ amplitude")
        elif mode == "triangles":
            for i, j in combinations(range(7), 2):
                _draw_edge(i, j, color="#bbbbbb", lw=0.7, alpha=0.7, zorder=1)
            palette = ["#c0392b", "#27ae60", "#2980b9", "#f39c12",
                        "#9b59b6", "#16a085", "#e67e22"]
            for n, (i, j, k) in enumerate(K7_heffter_triangles()):
                # smallest periodic-image triple
                best_a = np.inf; best_pts = None
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
                                    best_pts = (coords[i], pj, pk)
                p_i, p_j, p_k = best_pts
                poly = plt.Polygon([p_i, p_j, p_k],
                                     facecolor=palette[n % len(palette)],
                                     edgecolor="#1a1a1a", alpha=0.55,
                                     linewidth=0.6, zorder=2)
                ax.add_patch(poly)
            default_title = ("$K_{7}$ Heffter triangulation"
                              "  →  14 faces, $\\sum |[\\,,\\,,]\\,|^{2} = 112$")
        elif mode == "cycles":
            for i, j in combinations(range(7), 2):
                _draw_edge(i, j, color="#bbbbbb", lw=0.7, alpha=0.8, zorder=1)
            for cyc, col in [([0, 1, 3, 5], "#c0392b"),
                              ([0, 2, 4, 6], "#2980b9")]:
                for a, b in zip(cyc, cyc[1:] + [cyc[0]]):
                    _draw_edge(a, b, color=col, lw=2.4, alpha=0.9, zorder=3)
            default_title = ("$K_{7}$ 4-cycles"
                              "  →  cycle hierarchy $21 / 112 / 24$")
        else:
            raise ValueError(f"unknown mode {mode!r}; "
                              "use 'eulerian', 'triangles', or 'cycles'")

        # vertices on top
        for k, (x, y) in enumerate(coords):
            ax.plot(x, y, "o", markersize=11, markerfacecolor="#f4ecd9",
                     markeredgecolor="#1a1a1a", markeredgewidth=1.2,
                     zorder=5)
            ax.text(x, y, str(k), ha="center", va="center",
                     fontsize=8.5, fontweight="bold", zorder=6)

        ax.set_title(title or default_title, fontsize=12, pad=10)
        ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.18, 1.18)
        ax.set_aspect("equal"); ax.axis("off")
        fig.tight_layout()
    return _save(fig, save_to)


def figure_K7_all_traversals(*, save_to: Optional[str] = None,
                               figsize: tuple[float, float] = (18.0, 7.5)):
    """
    Three-panel figure showing all three K_7 traversal modes side-by-side.

    Useful for the "K_7 isn't just a graph -- it's a substrate object whose
    traversals compute different physics" framing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Render each mode separately, then composite
    with plt.rc_context(_PUB_STYLE):
        fig = plt.figure(figsize=figsize)
        for idx, mode in enumerate(("eulerian", "triangles", "cycles")):
            ax = fig.add_subplot(1, 3, idx + 1)
            # Re-use figure_K7_traversals's drawing logic by calling it
            # to a temporary figure, then copying axis content.  Cleaner:
            # inline the per-mode drawing here.  For brevity we reproduce
            # by building each panel directly.
            sub = figure_K7_traversals(mode=mode)
            # Move drawn artists from the sub-figure's axes to ours.
            sub_ax = sub.axes[0]
            for line in sub_ax.lines:
                ax.plot(line.get_xdata(), line.get_ydata(),
                         color=line.get_color(),
                         linewidth=line.get_linewidth(),
                         alpha=line.get_alpha(), zorder=line.zorder,
                         marker=line.get_marker(),
                         markersize=line.get_markersize(),
                         markerfacecolor=line.get_markerfacecolor(),
                         markeredgecolor=line.get_markeredgecolor(),
                         markeredgewidth=line.get_markeredgewidth())
            for patch in sub_ax.patches:
                from matplotlib.patches import Polygon
                if isinstance(patch, Polygon):
                    ax.add_patch(Polygon(patch.get_xy(),
                                            facecolor=patch.get_facecolor(),
                                            edgecolor=patch.get_edgecolor(),
                                            alpha=patch.get_alpha(),
                                            linewidth=patch.get_linewidth(),
                                            zorder=patch.zorder))
            for txt in sub_ax.texts:
                ax.text(txt.get_position()[0], txt.get_position()[1],
                         txt.get_text(),
                         fontsize=txt.get_fontsize(),
                         fontweight=txt.get_fontweight(),
                         ha=txt.get_ha(), va=txt.get_va(), zorder=txt.zorder)
            ax.set_title(sub_ax.get_title(), fontsize=12, pad=10)
            ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.18, 1.18)
            ax.set_aspect("equal"); ax.axis("off")
            import matplotlib.pyplot as plt2
            plt2.close(sub)

        fig.suptitle(r"$K_{7}$ traversal modes: same graph, different sectors",
                      fontsize=14, y=0.99)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save(fig, save_to)


def figure_paper18_unified(*, save_to: Optional[str] = None,
                             figsize: tuple[float, float] = (13.5, 6.3)):
    """
    The Paper 18 §4.4 figure: trefoil + K_7 on the same Heegaard torus,
    side-by-side.  Equivalent to figure_unified(p=2, q=3, layout='side-by-side').
    """
    return figure_unified(p=2, q=3, layout="side-by-side",
                            save_to=save_to, figsize=figsize)
