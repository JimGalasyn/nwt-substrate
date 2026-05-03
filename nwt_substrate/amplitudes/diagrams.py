"""
Feynman-diagram renderer using matplotlib.

Renders standard tree-level diagrams alongside their substrate-algebra
expressions.  Designed to produce publication-quality figures for Paper 18
and pedagogical materials.

Line styles:
  fermion : straight arrow (->)
  photon  : sinusoidal wavy line
  W/Z     : sinusoidal wavy line, dashed/double for distinction
  gluon   : helical / spring-like coil
  scalar  : dashed straight line

Each render_<process>() function returns a matplotlib Figure with the
diagram drawn; standalone scripts can save these to PNG/PDF.

Usage:

    fig = diagrams.render_compton()
    fig.savefig("compton.png", dpi=150)
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


# =============================================================================
# Primitive line drawers
# =============================================================================

def draw_fermion(ax, start, end, color="black", lw=1.6, arrow=True, label=None):
    """Straight line with arrowhead at midpoint."""
    x0, y0 = start
    x1, y1 = end
    ax.plot([x0, x1], [y0, y1], "-", color=color, lw=lw, zorder=1)
    if arrow:
        # Arrowhead at midpoint
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx = (x1 - x0) * 0.08
        dy = (y1 - y0) * 0.08
        ax.annotate("", xy=(mx + dx, my + dy), xytext=(mx - dx, my - dy),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw))
    if label:
        # Label slightly off the midpoint
        nx, ny = -(y1 - y0), (x1 - x0)
        n = np.hypot(nx, ny)
        if n > 0:
            nx, ny = nx / n * 0.10, ny / n * 0.10
        ax.text((x0 + x1) / 2 + nx, (y0 + y1) / 2 + ny, label,
                fontsize=11, ha="center", va="center")


def draw_photon(ax, start, end, color="darkblue", lw=1.6, n_waves=6,
                amplitude=0.06, label=None):
    """Sinusoidal wavy line for photon."""
    x0, y0 = start
    x1, y1 = end
    L = np.hypot(x1 - x0, y1 - y0)
    if L < 1e-9:
        return
    # Parametrize along line, with perpendicular wave
    t = np.linspace(0, 1, 240)
    # Tangent
    tx, ty = (x1 - x0), (y1 - y0)
    # Normal
    nx, ny = -ty / L, tx / L
    wave = amplitude * np.sin(2 * np.pi * n_waves * t)
    xs = x0 + t * tx + wave * nx
    ys = y0 + t * ty + wave * ny
    ax.plot(xs, ys, "-", color=color, lw=lw, zorder=1)
    if label:
        nx2, ny2 = -ty, tx
        n2 = np.hypot(nx2, ny2)
        if n2 > 0:
            nx2, ny2 = nx2 / n2 * 0.15, ny2 / n2 * 0.15
        ax.text((x0 + x1) / 2 + nx2, (y0 + y1) / 2 + ny2, label,
                fontsize=11, ha="center", va="center", color=color)


def draw_gluon(ax, start, end, color="forestgreen", lw=1.6, n_loops=8,
               radius=0.06, label=None):
    """Helical coil for gluon."""
    x0, y0 = start
    x1, y1 = end
    L = np.hypot(x1 - x0, y1 - y0)
    if L < 1e-9:
        return
    t = np.linspace(0, 1, 400)
    tx, ty = (x1 - x0), (y1 - y0)
    nx, ny = -ty / L, tx / L
    # Loop at frequency n_loops; r * cos along axis, r * sin off-axis
    phase = 2 * np.pi * n_loops * t
    along = (radius * 0.5) * (1 - np.cos(phase))
    off = radius * np.sin(phase)
    xs = x0 + t * tx + along * (tx / L) + off * nx
    ys = y0 + t * ty + along * (ty / L) + off * ny
    ax.plot(xs, ys, "-", color=color, lw=lw, zorder=1)
    if label:
        nx2, ny2 = -ty, tx
        n2 = np.hypot(nx2, ny2)
        if n2 > 0:
            nx2, ny2 = nx2 / n2 * 0.15, ny2 / n2 * 0.15
        ax.text((x0 + x1) / 2 + nx2, (y0 + y1) / 2 + ny2, label,
                fontsize=11, ha="center", va="center", color=color)


def draw_W(ax, start, end, color="darkred", lw=1.7, n_waves=6,
           amplitude=0.06, label=None):
    """Wavy line for W boson, drawn dashed to distinguish from photon."""
    x0, y0 = start
    x1, y1 = end
    L = np.hypot(x1 - x0, y1 - y0)
    if L < 1e-9:
        return
    t = np.linspace(0, 1, 240)
    tx, ty = (x1 - x0), (y1 - y0)
    nx, ny = -ty / L, tx / L
    wave = amplitude * np.sin(2 * np.pi * n_waves * t)
    xs = x0 + t * tx + wave * nx
    ys = y0 + t * ty + wave * ny
    ax.plot(xs, ys, "--", color=color, lw=lw, zorder=1)
    if label:
        nx2, ny2 = -ty, tx
        n2 = np.hypot(nx2, ny2)
        if n2 > 0:
            nx2, ny2 = nx2 / n2 * 0.15, ny2 / n2 * 0.15
        ax.text((x0 + x1) / 2 + nx2, (y0 + y1) / 2 + ny2, label,
                fontsize=11, ha="center", va="center", color=color)


def draw_vertex(ax, xy, color="black", size=12):
    """Small filled circle at a vertex."""
    ax.plot(xy[0], xy[1], "o", color=color, markersize=size * 0.4, zorder=2)


# =============================================================================
# Specific process diagrams (s-channel + u-channel where applicable)
# =============================================================================

def render_compton(ax=None, channel="both"):
    """
    Compton scattering e- + gamma -> e- + gamma.

    channel = "s" : s-channel only
    channel = "u" : u-channel only
    channel = "both" : both side-by-side
    """
    if channel == "both":
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        render_compton(axes[0], "s")
        render_compton(axes[1], "u")
        for a, t in zip(axes, ["s-channel", "u-channel"]):
            a.set_title(t, fontsize=11)
        fig.suptitle(r"$e^{-}\gamma \to e^{-}\gamma$  (Compton)", fontsize=12)
        fig.tight_layout()
        return fig

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))

    if channel == "s":
        # s-channel:  e- + gamma -> internal e -> e- + gamma
        # incoming e- bottom-left, incoming gamma top-left
        # outgoing e- bottom-right, outgoing gamma top-right
        # vertices on horizontal axis
        v1, v2 = (0.35, 0.5), (0.65, 0.5)
        in_e = (0.05, 0.2)
        in_g = (0.05, 0.8)
        out_e = (0.95, 0.2)
        out_g = (0.95, 0.8)
        draw_fermion(ax, in_e, v1, label=r"$e^{-}(p)$")
        draw_photon(ax, in_g, v1, label=r"$\gamma(k)$")
        draw_fermion(ax, v1, v2, label=r"$e^{*}(p+k)$")
        draw_fermion(ax, v2, out_e, label=r"$e^{-}(p')$")
        draw_photon(ax, v2, out_g, label=r"$\gamma(k')$")
        draw_vertex(ax, v1)
        draw_vertex(ax, v2)
    elif channel == "u":
        # u-channel: photons crossed
        v1, v2 = (0.35, 0.5), (0.65, 0.5)
        in_e = (0.05, 0.2)
        in_g = (0.05, 0.8)
        out_e = (0.95, 0.2)
        out_g = (0.95, 0.8)
        # Incoming photon -> top vertex; outgoing photon connects from same top vertex
        # but with inputs/outputs swapped
        draw_fermion(ax, in_e, v1, label=r"$e^{-}(p)$")
        draw_photon(ax, in_g, v2, label=r"$\gamma(k)$")
        draw_fermion(ax, v1, v2, label=r"$e^{*}(p-k')$")
        draw_fermion(ax, v2, out_e, label=r"$e^{-}(p')$")
        draw_photon(ax, v1, out_g, label=r"$\gamma(k')$")
        draw_vertex(ax, v1)
        draw_vertex(ax, v2)
    else:
        raise ValueError(f"Unknown channel {channel!r}")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    if standalone:
        return fig
    return ax


def render_eemumu(ax=None):
    """e+e- -> mu+mu- via s-channel virtual photon."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))

    v1, v2 = (0.5, 0.65), (0.5, 0.35)
    in_e_minus = (0.05, 0.85)
    in_e_plus  = (0.05, 0.15)
    out_mu_minus = (0.95, 0.85)
    out_mu_plus  = (0.95, 0.15)

    draw_fermion(ax, in_e_minus, v1, label=r"$e^{-}$")
    draw_fermion(ax, v1, in_e_plus, label=r"$e^{+}$")    # arrow reversed = positron
    draw_fermion(ax, v2, out_mu_minus, label=r"$\mu^{-}$")
    draw_fermion(ax, out_mu_plus, v2, label=r"$\mu^{+}$")
    draw_photon(ax, v1, v2, label=r"$\gamma^{*}$", n_waves=4)
    draw_vertex(ax, v1)
    draw_vertex(ax, v2)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    if standalone:
        return fig
    return ax


def render_muon_decay(ax=None):
    """mu- -> e- + nu_mu + nu_e_bar (V-A weak)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))

    v1, v2 = (0.4, 0.5), (0.7, 0.5)
    in_mu = (0.05, 0.5)
    out_numu = (0.6, 0.9)
    out_e = (0.95, 0.7)
    out_nuebar = (0.95, 0.3)

    draw_fermion(ax, in_mu, v1, label=r"$\mu^{-}$")
    draw_fermion(ax, v1, out_numu, label=r"$\nu_{\mu}$")
    draw_W(ax, v1, v2, label=r"$W^{-}$", n_waves=4)
    draw_fermion(ax, v2, out_e, label=r"$e^{-}$")
    draw_fermion(ax, out_nuebar, v2, label=r"$\bar{\nu}_{e}$")
    draw_vertex(ax, v1)
    draw_vertex(ax, v2)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    if standalone:
        return fig
    return ax


def render_neutron_decay(ax=None):
    """n -> p + e- + nu_e_bar (V-A weak, quark level: d -> u + W- -> u + e- nu_e_bar)."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))

    v1, v2 = (0.4, 0.5), (0.7, 0.5)
    in_d = (0.05, 0.5)
    out_u = (0.6, 0.9)
    out_e = (0.95, 0.7)
    out_nuebar = (0.95, 0.3)

    draw_fermion(ax, in_d, v1, label=r"$d$")
    draw_fermion(ax, v1, out_u, label=r"$u$")
    draw_W(ax, v1, v2, label=r"$W^{-}$", n_waves=4)
    draw_fermion(ax, v2, out_e, label=r"$e^{-}$")
    draw_fermion(ax, out_nuebar, v2, label=r"$\bar{\nu}_{e}$")
    draw_vertex(ax, v1)
    draw_vertex(ax, v2)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    if standalone:
        return fig
    return ax


def render_moller(ax=None, channel="both"):
    """Møller scattering e- + e- -> e- + e- (t-channel + u-channel)."""
    if channel == "both":
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        render_moller(axes[0], "t")
        render_moller(axes[1], "u")
        for a, t in zip(axes, ["t-channel", "u-channel"]):
            a.set_title(t, fontsize=11)
        fig.suptitle(r"$e^{-}e^{-} \to e^{-}e^{-}$  (Møller)", fontsize=12)
        fig.tight_layout()
        return fig

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5, 4))

    if channel == "t":
        v1, v2 = (0.4, 0.7), (0.4, 0.3)
        out_v1, out_v2 = (0.95, 0.7), (0.95, 0.3)
        in_v1, in_v2 = (0.05, 0.7), (0.05, 0.3)
        draw_fermion(ax, in_v1, v1, label=r"$e^{-}$")
        draw_fermion(ax, v1, out_v1, label=r"$e^{-}$")
        draw_fermion(ax, in_v2, v2, label=r"$e^{-}$")
        draw_fermion(ax, v2, out_v2, label=r"$e^{-}$")
        draw_photon(ax, v1, v2, label=r"$\gamma^{*}$", n_waves=4)
    else:  # u-channel: cross outgoing
        v1, v2 = (0.4, 0.7), (0.4, 0.3)
        in_v1, in_v2 = (0.05, 0.7), (0.05, 0.3)
        # cross output
        out_v1 = (0.95, 0.3)   # was (0.95, 0.7)
        out_v2 = (0.95, 0.7)
        draw_fermion(ax, in_v1, v1, label=r"$e^{-}$")
        draw_fermion(ax, v1, out_v1, label=r"$e^{-}$")
        draw_fermion(ax, in_v2, v2, label=r"$e^{-}$")
        draw_fermion(ax, v2, out_v2, label=r"$e^{-}$")
        draw_photon(ax, v1, v2, label=r"$\gamma^{*}$", n_waves=4)
    draw_vertex(ax, v1)
    draw_vertex(ax, v2)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    if standalone:
        return fig
    return ax


# =============================================================================
# Substrate-algebra expression LaTeX strings
# =============================================================================

EXPRESSIONS = {
    "compton": (
        r"$i\mathcal{M} = -ie^2\,\bar{u}(p')\,[\,"
        r"\gamma^{\nu}\epsilon'^{*}_{\nu}\, S_F(p+k)\, \gamma^{\mu}\epsilon_{\mu}"
        r"\;+\;\gamma^{\mu}\epsilon_{\mu}\, S_F(p-k')\, \gamma^{\nu}\epsilon'^{*}_{\nu}"
        r"\,]\,u(p)$"
    ),
    "eemumu": (
        r"$i\mathcal{M} = \frac{i\,e^2}{s}\,\bar{v}(p_{e^+})\gamma^{\mu}u(p_{e^-})\,"
        r"\bar{u}(p_{\mu^-})\gamma_{\mu}v(p_{\mu^+})$"
    ),
    "muon_decay": (
        r"$i\mathcal{M} = \frac{-i\,g^2}{8 M_W^2}\,"
        r"\bar{u}(p_{\nu_\mu})\gamma^{\mu}(1-\gamma^5)u(p_{\mu})\,"
        r"\bar{u}(p_e)\gamma_{\mu}(1-\gamma^5)v(p_{\bar{\nu}_e})$"
    ),
    "neutron_decay": (
        r"$i\mathcal{M} = \frac{-i\,G_F}{\sqrt{2}}\,V_{ud}\,"
        r"\bar{u}_p\gamma^{\mu}(g_V - g_A\gamma^5)u_n\,"
        r"\bar{u}_e\gamma_{\mu}(1-\gamma^5)v_{\bar{\nu}_e}$"
    ),
    "moller": (
        r"$i\mathcal{M} = \frac{i\,e^2}{t}\,[\bar{u}_3\gamma^{\mu}u_1][\bar{u}_4\gamma_{\mu}u_2]"
        r"\;-\;\frac{i\,e^2}{u}\,[\bar{u}_4\gamma^{\mu}u_1][\bar{u}_3\gamma_{\mu}u_2]$"
    ),
}


def render_with_expression(name, render_fn, **kwargs):
    """
    Render a process diagram with its substrate-algebra expression
    underneath as a caption.

    Returns matplotlib Figure.
    """
    fig = render_fn(**kwargs)
    if name in EXPRESSIONS:
        fig.text(0.5, 0.02, EXPRESSIONS[name], ha="center", fontsize=10)
        fig.subplots_adjust(bottom=0.18)
    return fig
