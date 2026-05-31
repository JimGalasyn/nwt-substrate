"""Particle portraits — render an NWT particle from its field configuration.

A particle's carrier knot (Paper 6 / Paper 8) is realised as a BPS vortex
defect of the abelian-Higgs condensate and rendered as a glowing,
phase-coloured filament:

    >>> import nwt_substrate.portraits as pp
    >>> fig = pp.portrait(2, 3, title="proton (trefoil)")   # doctest: +SKIP
    >>> fig = pp.portrait(hopf=True, title="meson")          # doctest: +SKIP
    >>> fig = pp.gallery(save_to="carriers.png")             # doctest: +SKIP

The geometry comes from :func:`nwt_substrate.diagrams.torus_knot_curve`, the
scalar cross-section from :func:`nwt_substrate.condensate.solve_bps_vortex`,
and the phase from the abelian-Higgs winding; see
:mod:`nwt_substrate.portraits.field`.
"""
from __future__ import annotations

from collections import OrderedDict

from .field import ParticleField, build_particle_field
from .render import render_field, save_portrait

__all__ = [
    "ParticleField",
    "build_particle_field",
    "render_field",
    "save_portrait",
    "CARRIERS",
    "n_q_to_knot",
    "portrait",
    "carrier_portrait",
    "gallery",
]


# Curated carrier zoo (the toroidal carriers + the Hopf-link meson).  Each
# value is the kwargs handed to :func:`build_particle_field`.
CARRIERS: "OrderedDict[str, dict]" = OrderedDict([
    ("lepton — unknot (2,1)", dict(p=2, q=1)),
    ("meson — Hopf link", dict(hopf=True)),
    ("baryon — trefoil 3₁ (2,3)", dict(p=2, q=3)),
    ("nucleon — cinquefoil 5₁ (2,5)", dict(p=2, q=5)),
    ("septafoil 7₁ (2,7)", dict(p=2, q=7)),
])


def n_q_to_knot(n_q: int) -> dict:
    """Map a carrier crossing number ``n_q`` to field-builder kwargs.

    Supports the toroidal carrier family T(2, n_q) for odd ``n_q`` plus the
    unknot (n_q 0/1) and the Hopf link (n_q 2).  The figure-eight (4₁) and
    6-crossing carriers are genuinely non-toroidal and are not yet supported.
    """
    if n_q in (0, 1):
        return dict(p=2, q=1)
    if n_q == 2:
        return dict(hopf=True)
    if n_q >= 3 and n_q % 2 == 1:
        return dict(p=2, q=n_q)
    raise ValueError(
        f"n_q={n_q} is a non-toroidal carrier (e.g. figure-eight 4₁); "
        f"portrait support is currently limited to the T(2,q) family and "
        f"the Hopf link."
    )


def portrait(p: int = 2, q: int = 3, *, hopf: bool = False,
             title: str | None = None, save_to=None, N: int = 190,
             cmap: str = "hsv", k_opacity: float = 7.0, gain: float = 1.7,
             gamma: float = 1.0, axis: int = 2, **field_kwargs):
    """Build and render one particle portrait; return the matplotlib Figure.

    Geometry/field keywords (``box``, ``R``, ``r``, ``xi``, ``glow_width``,
    ``tilt`` ...) are forwarded to :func:`build_particle_field`; rendering
    keywords (``cmap``, ``k_opacity``, ``gain``, ``gamma``, ``axis``) to
    :func:`render_field`.  If ``save_to`` is given the figure is also written
    to disk.
    """
    pf = build_particle_field(p, q, hopf=hopf, N=N, **field_kwargs)
    fig = save_portrait(
        pf, save_to, title=title,
        cmap=cmap, k_opacity=k_opacity, gain=gain, gamma=gamma, axis=axis,
    )
    return fig


def carrier_portrait(n_q: int, *, title: str | None = None, **kwargs):
    """Render the canonical carrier for crossing number ``n_q``."""
    spec = n_q_to_knot(n_q)
    return portrait(title=title, **spec, **kwargs)


def gallery(specs: "OrderedDict[str, dict] | None" = None, *, ncols: int = 3,
            N: int = 160, save_to=None, cmap: str = "hsv",
            facecolor: str = "black", **render_kwargs):
    """Render a contact sheet of carrier portraits; return the Figure.

    ``specs`` maps titles to :func:`build_particle_field` kwargs (defaults to
    :data:`CARRIERS`).
    """
    import matplotlib.pyplot as plt
    from .render import render_field, _orient

    specs = specs if specs is not None else CARRIERS
    items = list(specs.items())
    nrows = (len(items) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 4.2 * nrows))
    axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for ax, (title, spec) in zip(axes, items):
        pf = build_particle_field(N=N, **spec)
        ax.imshow(_orient(render_field(pf, cmap=cmap, **render_kwargs)),
                  origin="upper")
        ax.set_title(title, color="#e8e8e8", fontsize=11, pad=6)
    for ax in axes[len(items):]:
        ax.set_visible(False)
    for ax in axes:
        ax.axis("off")
    fig.patch.set_facecolor(facecolor)
    fig.tight_layout()
    if save_to is not None:
        fig.savefig(save_to, dpi=150, bbox_inches="tight",
                    facecolor=facecolor, pad_inches=0.1)
    return fig
