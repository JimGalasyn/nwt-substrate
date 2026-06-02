"""Periodic cubic coordinate grid for lattice soliton fields.

A tiny, rendering-free, numpy-only mirror of the JAX `Grid` used by the soliton
relaxer engine, so the engine and this reference library agree on coordinate and
wavenumber conventions exactly (centered, endpoint-excluded sample points;
angular wavenumbers from `fftfreq`).  Shared across the soliton sectors
(`faddeev` now, `skyrme` later).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BoxGrid:
    """Periodic cubic grid of N**3 points on a box of side L (in healing-length
    units).  Coordinates are centered and exclude the right endpoint, so a field
    sampled here is periodic with period L."""

    N: int
    L: float

    @property
    def dx(self) -> float:
        return self.L / self.N

    def coords(self):
        """(X, Y, Z) meshgrids, centered on the box, ij-indexed."""
        c = np.linspace(-self.L / 2, self.L / 2, self.N, endpoint=False)
        return np.meshgrid(c, c, c, indexing="ij")

    def k_vectors(self):
        """(KX, KY, KZ, K2) angular wavenumbers (2*pi/L spacing), ij-indexed."""
        k1 = 2.0 * np.pi * np.fft.fftfreq(self.N, d=self.dx)
        KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
        return KX, KY, KZ, KX**2 + KY**2 + KZ**2
