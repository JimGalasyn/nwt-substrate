"""Smooth, compact, in-basin hopfion initial conditions (rational maps).

`rational_hopfion` builds the CP^1 spinor directly as a composition of smooth
maps -- the rational-map construction (Battye-Sutcliffe / Hietarinta-Salo) --
avoiding the two defects of the analytic stereographic seed that push it OUTSIDE
the soliton basin:

  (1) the single-patch lift Z1 = (n1 - i n2)/sqrt(2(1-n3)) is 0/0 on the core
      ring (n3 -> +1), a spinor singularity sitting on the soliton core;
  (2) the field reaches vacuum only as 1/r^2, so on a periodic box it carries a
      few-percent seam mismatch.

Construction (charge Q_H = n*m):

  base map (inverse stereographic R^3 -> S^3), smooth angles
      phi1 = atan2(y, x)                      azimuthal, about the z-axis
      phi2 = atan2(2 R z, r^2 - R^2)          meridional, around the core ring
  rational (degree-(n, m)) map on the Hopf fibre, smooth because each amplitude
  vanishes where its phase is ill-defined (like z^k at z=0):
      Z1 = cos(lam) e^{i n phi1},   Z2 = sin(lam) e^{i m phi2}
  compact profile in the minor radius d = sqrt((rho - R)^2 + z^2):
      lam(d) = (pi/2) * smootherstep(d / w),   lam = pi/2 (vacuum) for d >= w,

so the field is EXACTLY the south-pole vacuum beyond d = w (periodic-safe, with
a vacuum margin to the box wall).  Regularity: on the z-axis d >= R >= w so
cos(lam)=0 kills the ill-defined phi1; on the core ring d=0 so sin(lam)=0 kills
the ill-defined phi2.  Hence w <= R is required.  The core ring has major radius
R in the z=0 plane and tube radius w.
"""
from __future__ import annotations

import numpy as np


def rational_hopfion(grid, R=3.5, w=None, n=1, m=1, center=(0.0, 0.0, 0.0)):
    """Smooth rational-map hopfion spinor (Z1, Z2) of charge Q_H = n*m.

    grid is a `BoxGrid`.  R is the core-ring major radius, w the tube radius
    (default 0.85 R; must be <= R).  Returns complex arrays (Z1, Z2) with
    |Z1|^2 + |Z2|^2 = 1."""
    if w is None:
        w = 0.85 * R
    if w > R:
        raise ValueError(f"need w <= R for axis regularity (w={w}, R={R})")
    X, Y, Z = (np.asarray(c) for c in grid.coords())
    x0, y0, z0 = center
    x, y, z = X - x0, Y - y0, Z - z0
    rho = np.sqrt(x**2 + y**2)
    r2 = x**2 + y**2 + z**2
    phi1 = np.arctan2(y, x)
    phi2 = np.arctan2(2.0 * R * z, r2 - R**2)        # oriented so Q_H = +n*m
    d = np.sqrt((rho - R) ** 2 + z**2)
    t = np.clip(d / w, 0.0, 1.0)
    s = t**3 * (10.0 - 15.0 * t + 6.0 * t**2)        # smootherstep (C^2)
    lam = 0.5 * np.pi * s
    Z1 = np.cos(lam) * np.exp(1j * n * phi1)
    Z2 = np.sin(lam) * np.exp(1j * m * phi2)
    nrm = np.sqrt(np.abs(Z1) ** 2 + np.abs(Z2) ** 2)
    return (Z1 / nrm).astype(np.complex128), (Z2 / nrm).astype(np.complex128)
