"""Render a :class:`ParticleField` as an emissive volume projection.

Each voxel emits light coloured by the local condensate phase (a cyclic
colour map) with an opacity set by the emission field; we composite
front-to-back along the view axis (emission-absorption integral), so the
defect reads as a glowing knotted filament with correct over/under-crossing
occlusion.  The result is a plain RGB image array; no physics depends on the
choice of colour map or tone curve.
"""
from __future__ import annotations

import numpy as np
import matplotlib

from .field import ParticleField

__all__ = ["render_field", "save_portrait"]


def render_field(
    pf: ParticleField,
    *,
    cmap: str = "hsv",
    k_opacity: float = 7.0,
    gain: float = 1.7,
    gamma: float = 1.0,
    axis: int = 2,
) -> np.ndarray:
    """Composite a particle field into an ``(H, W, 3)`` RGB image in [0, 1].

    Parameters
    ----------
    pf : ParticleField
        The field to render.
    cmap : str
        Cyclic matplotlib colour map applied to the phase (``hsv``,
        ``twilight``, ``twilight_shifted`` are good choices).
    k_opacity : float
        Opacity gain: per-voxel alpha = 1 - exp(-k_opacity * glow).
    gain : float
        Brightness gain in the final tone map ``1 - exp(-gain * I)``.
    gamma : float
        Optional gamma applied after tone mapping (1.0 = none).
    axis : int
        View/projection axis (0, 1, or 2).
    """
    phase01 = (pf.phase % (2 * np.pi)) / (2 * np.pi)
    rgb = matplotlib.colormaps[cmap](phase01)[..., :3]
    alpha = 1.0 - np.exp(-k_opacity * pf.glow)

    # Move the view axis to the front, then composite front-to-back.
    a = np.moveaxis(alpha, axis, 0)
    c = np.moveaxis(rgb, axis, 0)
    ones = np.ones_like(a[:1])
    transmittance = np.cumprod(
        np.concatenate([ones, 1.0 - a[:-1]], axis=0), axis=0
    )
    img = (transmittance[..., None] * a[..., None] * c).sum(axis=0)

    img = 1.0 - np.exp(-gain * img)
    if gamma != 1.0:
        img = np.clip(img, 0, 1) ** gamma
    return np.clip(img, 0.0, 1.0)


def _orient(img: np.ndarray) -> np.ndarray:
    """Orient a composited image for natural display (image-row = -y)."""
    return np.transpose(img, (1, 0, 2))[::-1]


def save_portrait(
    pf: ParticleField,
    path=None,
    *,
    title: str | None = None,
    dpi: int = 150,
    facecolor: str = "black",
    title_color: str = "#e8e8e8",
    **render_kwargs,
):
    """Render ``pf`` into a square portrait Figure; save it if ``path`` given.

    Extra keyword arguments are passed to :func:`render_field`.  Returns the
    matplotlib Figure (whether or not it was written to disk).
    """
    import matplotlib.pyplot as plt

    img = _orient(render_field(pf, **render_kwargs))
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(img, origin="upper")
    ax.axis("off")
    fig.patch.set_facecolor(facecolor)
    if title:
        ax.set_title(title, color=title_color, fontsize=12, pad=8)
    if path is not None:
        fig.savefig(path, dpi=dpi, bbox_inches="tight",
                    facecolor=facecolor, pad_inches=0.05)
    return fig
