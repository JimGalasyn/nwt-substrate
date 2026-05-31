"""Build the 3-D condensate field of a particle from its carrier knot.

An NWT particle is a vortex defect of the abelian-Higgs condensate whose core
follows the particle's carrier-knot curve (Paper 6 / Paper 8).  Here we realise
that field on a Cartesian grid:

  * the scalar amplitude ``|psi|`` is the exact BPS vortex profile ``f(rho)``
    (:mod:`nwt_substrate.condensate.vortex_profile`) evaluated at the distance
    to the nearest point of the carrier curve, in healing-length units;
  * the phase is the abelian-Higgs winding ``p*phi + q*theta`` in the
    embedding-torus angles, so the field carries the correct topological
    circulation;
  * an emission field ``glow`` (a thin Gaussian about the core) marks where the
    defect "lives" for rendering.

This is the physically-seeded configuration: the BPS cross-section is exact and
the topology is exact; relaxing it under the Gross-Pitaevskii flow would only
settle the bends, not change the class.
"""
from __future__ import annotations

from dataclasses import dataclass, field as _dc_field

import numpy as np
from scipy.spatial import cKDTree

from ..diagrams import torus_knot_curve, hopf_link_curves
from ..condensate.vortex_profile import solve_bps_vortex

__all__ = ["ParticleField", "build_particle_field"]


@dataclass
class ParticleField:
    """A particle's condensate field sampled on an ``N^3`` grid.

    Attributes
    ----------
    phase : np.ndarray
        Phase arg(psi) on the grid, shape ``(N, N, N)``.
    glow : np.ndarray
        Emission intensity (peaks on the vortex core), shape ``(N, N, N)``.
    amp : np.ndarray
        Scalar amplitude |psi| (0 on core, 1 in bulk), shape ``(N, N, N)``.
    box : float
        Half-width of the cubic domain (units of major radius R).
    meta : dict
        Descriptive metadata (p, q, carrier kind, winding, ...).
    """

    phase: np.ndarray
    glow: np.ndarray
    amp: np.ndarray
    box: float
    meta: dict = _dc_field(default_factory=dict)

    @property
    def N(self) -> int:
        return self.phase.shape[0]


def _rot3(ax: float, ay: float) -> np.ndarray:
    """Rotation: about x by ``ax`` then about y by ``ay`` (for a 3/4 view)."""
    cx, sx = np.cos(ax), np.sin(ax)
    cy, sy = np.cos(ay), np.sin(ay)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    return Ry @ Rx


def build_particle_field(
    p: int = 2,
    q: int = 3,
    *,
    hopf: bool = False,
    N: int = 190,
    box: float = 1.85,
    R: float = 1.0,
    r: float = 0.42,
    xi: float = 0.05,
    glow_width: float = 0.55,
    tilt: tuple[float, float] = (0.62, 0.5),
    curve_points: int = 4000,
) -> ParticleField:
    """Construct the condensate field for a carrier knot.

    Parameters
    ----------
    p, q : int
        Torus-knot winding numbers (ignored if ``hopf``).  (2, 1) unknot
        (lepton), (2, 3) trefoil (baryon), (2, 5) cinquefoil (nucleon), ...
    hopf : bool
        If True, build the two-component Hopf link (meson carrier) instead of
        a torus knot.
    N : int
        Grid points per dimension.
    box : float
        Half-width of the cubic domain in units of R.
    R, r : float
        Major and minor radii of the embedding torus.
    xi : float
        Vortex healing length (core width) in domain units.
    glow_width : float
        Emission width as a multiple of ``xi``.
    tilt : (float, float)
        (x, y) rotation of the scene, for a three-quarter view.
    curve_points : int
        Samples along each carrier curve (denser = smoother core distance).
    """
    g = np.linspace(-box, box, N)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1) @ _rot3(*tilt).T

    if hopf:
        curves = list(hopf_link_curves(R=R, r=r, n_points=curve_points))
        kind, wind = "Hopf link", (1, 1)
    else:
        curves = [torus_knot_curve(p, q, R=R, r=r, n_points=curve_points)]
        kind, wind = "torus knot", (p, q)

    # Distance from each grid point to the nearest carrier-curve point.
    dist = np.full(len(pts), np.inf)
    for c in curves:
        dd, _ = cKDTree(c).query(pts, k=1)
        dist = np.minimum(dist, dd)
    dist = dist.reshape(N, N, N)

    xr = pts[:, 0].reshape(N, N, N)
    yr = pts[:, 1].reshape(N, N, N)
    zr = pts[:, 2].reshape(N, N, N)
    phi = np.arctan2(yr, xr)
    rho_cyl = np.hypot(xr, yr)
    theta = np.arctan2(zr, rho_cyl - R)
    phase = wind[0] * phi + wind[1] * theta

    # Exact BPS scalar amplitude as a function of core distance (in xi units).
    profile = solve_bps_vortex(n=1)
    amp = profile.f_at(dist / xi)
    glow = np.exp(-((dist / (glow_width * xi)) ** 2))

    return ParticleField(
        phase=phase, glow=glow, amp=amp, box=box,
        meta={"p": wind[0], "q": wind[1], "kind": kind, "hopf": hopf,
              "xi": xi, "R": R, "r": r},
    )
