"""Nielsen-Olesen BPS vortex profile — the radial cross-section of a defect.

The abelian-Higgs condensate at its BPS critical point (lambda = e^2/2, see
:mod:`nwt_substrate.condensate.abelian_higgs`) admits a vortex of winding ``n``
whose scalar amplitude ``f(rho)`` and gauge field ``a(rho)`` satisfy the
first-order Bogomolny equations (in units where the healing length xi = 1):

    f'(rho) = (n - a) f / rho
    a'(rho) = rho (1 - f^2) / 2

with boundary conditions f(0) = 0, f(inf) = 1, a(0) = 0, a(inf) = n.  Near the
core the scalar rises as ``f ~ c1 rho^n``; the single shooting parameter is that
leading coefficient ``c1`` (which equals f'(0+) only for n = 1, where it
converges to ~0.6033), fixed by requiring f(rho_max) = 1.  At the BPS point the
scalar and gauge masses coincide (m_s = m_v = 1/xi), so the scalar approaches the
bulk exponentially at rate 1/xi: 1 - f ~ K0(rho) in these units (Nielsen &
Olesen 1973; Bogomolny 1976).

This is the exact profile that the abelian-Higgs scalar takes across a vortex
core; bend it along a carrier-knot curve and you have the field configuration
of an NWT particle (see :mod:`nwt_substrate.portraits`).

Only the unit vortex (n = 1) is solved here.  The first-order shooting is
stiff: for n >= 2 the trajectory overshoots and diverges past f = 1 before
reaching ``rho_max``, so a higher winding would need a relaxation / collocation
boundary-value solver (e.g. ``scipy.integrate.solve_bvp``).  The equations and
the small-rho initial condition below are written for general n so that
extension is a localised change.

References: Bogomolny (1976); Nielsen & Olesen (1973); Paper 16 (L_2 BPS
sector); Paper 6 (mass scale = line tension).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq

__all__ = [
    "BPSVortexProfile",
    "GLVortexProfile",
    "bps_rhs",
    "gl_vortex_rhs",
    "solve_bps_vortex",
    "solve_gl_vortex",
]


def bps_rhs(rho: float, y, n: int = 1):
    """Right-hand side of the first-order BPS vortex equations."""
    f, a = y
    return [(n - a) * f / rho, rho * (1.0 - f**2) / 2.0]


def _shoot(c1: float, rho_start: float, rho_end: float, n: int) -> float:
    """Integrate from rho_start with f ~ c1 rho^n; return f(rho_end) - 1."""
    y0 = [c1 * rho_start**n, rho_start**2 / 4.0]
    sol = solve_ivp(
        bps_rhs, [rho_start, rho_end], y0, args=(n,),
        method="RK45", rtol=1e-11, atol=1e-13, max_step=0.02,
    )
    if not sol.success:
        return 1e10
    return sol.y[0, -1] - 1.0


@dataclass(frozen=True)
class BPSVortexProfile:
    """A solved BPS vortex profile.

    Attributes
    ----------
    rho : np.ndarray
        Radial grid in units of the healing length.
    f : np.ndarray
        Scalar amplitude |psi|/v; rises from 0 at the core to 1 in the bulk.
    a : np.ndarray
        Gauge profile; rises from 0 to the winding n.
    c1 : float
        Leading small-rho coefficient (f ~ c1 rho^n); equals f'(0+) for n = 1.
    n : int
        Vortex winding number.
    """

    rho: np.ndarray
    f: np.ndarray
    a: np.ndarray
    c1: float
    n: int

    def f_at(self, r):
        """Scalar amplitude at radius ``r`` (clamped to bulk value 1 beyond
        the solved grid).  Accepts scalars or arrays."""
        return np.interp(r, self.rho, self.f, left=0.0, right=1.0)

    def a_at(self, r):
        """Gauge profile at radius ``r`` (clamped to n beyond the grid)."""
        return np.interp(r, self.rho, self.a, left=0.0, right=float(self.n))


@lru_cache(maxsize=16)
def _solve_cached(n: int, rho_max: float, n_points: int) -> BPSVortexProfile:
    rho_start = 1e-4
    # The leading coefficient c1 (f ~ c1 rho^n) shrinks with n, so scan from a
    # low floor for the single sign change of the shooting residual (negative
    # = undershoot, large positive = overshoot), then refine with brentq.
    grid = np.linspace(0.01, 3.0, 150)
    vals = [_shoot(c, rho_start, rho_max, n) for c in grid]
    sign_change = None
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            sign_change = (grid[i], grid[i + 1])
            break
    if sign_change is None:
        raise RuntimeError(
            f"Could not bracket BPS shooting parameter for n={n} in "
            f"c1 in [0.01, 3.0]."
        )
    c1 = brentq(
        _shoot, sign_change[0], sign_change[1], args=(rho_start, rho_max, n),
        xtol=1e-10, rtol=1e-10,
    )
    rho_grid = np.linspace(rho_start, rho_max, n_points)
    sol = solve_ivp(
        bps_rhs, [rho_start, rho_max],
        [c1 * rho_start**n, rho_start**2 / 4.0],
        args=(n,), t_eval=rho_grid, method="RK45",
        rtol=1e-12, atol=1e-14, max_step=0.01,
    )
    return BPSVortexProfile(rho=sol.t, f=sol.y[0], a=sol.y[1], c1=float(c1), n=n)


def solve_bps_vortex(
    n: int = 1, rho_max: float = 15.0, n_points: int = 3000,
) -> BPSVortexProfile:
    """Solve the BPS unit vortex (results are cached).

    Parameters
    ----------
    n : int
        Vortex winding number.  Only ``n = 1`` is currently supported: the
        first-order shooting is numerically unstable for higher windings
        (the trajectory blows up past f = 1 before ``rho_max``), which needs a
        stiff boundary-value solver rather than shooting — see the module note.
    rho_max : float
        Outer radius of the integration domain, in healing lengths.
    n_points : int
        Number of grid points returned.
    """
    if n < 1:
        raise ValueError("BPS vortex winding n must be >= 1.")
    if n != 1:
        raise NotImplementedError(
            "Only the unit vortex (n=1) is currently solved; multi-winding "
            "BPS profiles need a stiff BVP solver (shooting is unstable for "
            "n>=2) and are future work.  Use solve_gl_vortex(kappa=1.0, n=n), "
            "whose boundary-value collocation handles n>=2."
        )
    return _solve_cached(int(n), float(rho_max), int(n_points))


# ---------------------------------------------------------------------------
# General gauged Ginzburg-Landau vortex (any GL parameter kappa = lambda/xi)
# ---------------------------------------------------------------------------
# Away from the BPS critical point the first-order Bogomolny equations no
# longer hold; the static vortex obeys the full second-order Ginzburg-Landau /
# de Vega-Schaposnik system.  With length in *penetration-depth* (lambda)
# units and the internal coefficient kc = kappa / sqrt(2):
#
#     f'' + f'/rho - (n - a)^2 f / rho^2 + kc^2 (1 - f^2) f = 0       (scalar)
#     a'' - a'/rho - f^2 (n - a)                            = 0       (gauge)
#
# (m_s = sqrt(2) kc = kappa, m_v = 1 here, so kappa = m_s/m_v = lambda/xi.)
# kappa = 1 is the BPS self-dual point (m_s = m_v) and reduces to the
# first-order equations above; kappa < 1 is type-I, kappa > 1 type-II.  The
# first-order shooting in solve_bps_vortex is replaced by stiff boundary-value
# collocation (scipy.solve_bvp), which is also stable for higher winding n.
# Profiles are returned in *healing-length* (xi) units (rho_xi = kappa * rho_lambda),
# so solve_gl_vortex(kappa=1.0) matches solve_bps_vortex point-for-point.


def gl_vortex_rhs(rho, y, kappa: float = 1.0, n: int = 1):
    """Right-hand side of the second-order GL vortex system (lambda units).

    ``y = [f, f', a, a']``; ``kappa`` is the GL parameter ``lambda/xi``
    (BPS self-dual at ``kappa = 1``).  Accepts the vectorised ``rho`` that
    :func:`scipy.integrate.solve_bvp` passes, returning a ``(4, len(rho))``
    array.
    """
    kc = kappa / np.sqrt(2.0)
    f, fp, a, ap = y
    return np.vstack([
        fp,
        -fp / rho + (n - a) ** 2 * f / rho**2 - kc**2 * (1.0 - f**2) * f,
        ap,
        ap / rho - f**2 * (n - a),
    ])


@dataclass(frozen=True)
class GLVortexProfile:
    """A solved Ginzburg-Landau vortex profile at GL parameter ``kappa``.

    Attributes
    ----------
    rho : np.ndarray
        Radial grid in units of the healing length xi.
    f : np.ndarray
        Scalar amplitude |psi|/v; rises from 0 at the core to 1 in the bulk
        (over ~1 healing length).
    a : np.ndarray
        Gauge profile; rises from 0 to the winding n over ~kappa healing
        lengths (the penetration depth).
    kappa : float
        GL parameter lambda/xi.  kappa = 1 is the BPS self-dual point;
        kappa < 1 type-I, kappa > 1 type-II.
    n : int
        Vortex winding number.
    """

    rho: np.ndarray
    f: np.ndarray
    a: np.ndarray
    kappa: float
    n: int

    def f_at(self, r):
        """Scalar amplitude at radius ``r`` (0 at the core, clamped to 1
        beyond the solved grid).  Accepts scalars or arrays."""
        return np.interp(r, self.rho, self.f, left=0.0, right=1.0)

    def a_at(self, r):
        """Gauge profile at radius ``r`` (clamped to n beyond the grid)."""
        return np.interp(r, self.rho, self.a, left=0.0, right=float(self.n))

    def supercurrent(self, r=None):
        """Gauge-invariant azimuthal supercurrent ``f^2 (n - a) / rho``.

        This is the screening (London) current of the vortex; its radial
        sheath sits near the core at the BPS point and migrates outward
        (toward ~kappa healing lengths) as kappa increases into the type-II
        regime.  Evaluated on the solved grid if ``r`` is None, else at ``r``.
        """
        if r is None:
            rho, f, a = self.rho, self.f, self.a
        else:
            r = np.asarray(r, dtype=float)
            rho, f, a = r, self.f_at(r), self.a_at(r)
        return f**2 * (self.n - a) / np.maximum(rho, 1e-12)


@lru_cache(maxsize=32)
def _solve_gl_cached(kappa: float, n: int, rho_max: float, n_points: int) -> GLVortexProfile:
    # Solve in penetration-depth (lambda) units, then rescale to healing
    # lengths.  Float the outer radius up so the gauge field saturates even
    # for large kappa (penetration ~ a few lambda).
    L_lam = max(rho_max / kappa, 12.0)
    rho0 = 1e-3
    r = np.linspace(rho0, L_lam, max(n_points, 1200))
    # Initial guess: scalar heals over ~1/kappa (lambda units), gauge over ~1.
    y0 = np.zeros((4, r.size))
    y0[0] = np.tanh(kappa * r)
    y0[1] = kappa / np.cosh(kappa * r) ** 2
    y0[2] = n * (1.0 - np.exp(-(r**2)))
    y0[3] = 2.0 * n * r * np.exp(-(r**2))

    def bc(ya, yb):
        return np.array([ya[0], ya[2], yb[0] - 1.0, yb[2] - n])

    sol = solve_bvp(
        lambda rr, yy: gl_vortex_rhs(rr, yy, kappa, n),
        bc, r, y0, max_nodes=200000, tol=1e-7,
    )
    if not sol.success:
        raise RuntimeError(
            f"Ginzburg-Landau vortex BVP failed to converge for kappa={kappa}, "
            f"n={n} ({sol.message}).  Very large kappa may need continuation "
            "(solve ascending kappa, seeding each from the previous solution)."
        )
    rho_xi = kappa * sol.x
    return GLVortexProfile(rho=rho_xi, f=sol.y[0], a=sol.y[2], kappa=float(kappa), n=int(n))


def solve_gl_vortex(
    kappa: float = 1.0, n: int = 1, rho_max: float = 20.0, n_points: int = 2000,
) -> GLVortexProfile:
    """Solve the gauged Ginzburg-Landau vortex at GL parameter ``kappa`` (cached).

    Unlike :func:`solve_bps_vortex` (locked to the BPS self-dual point via
    first-order shooting), this solves the full second-order GL / de Vega-
    Schaposnik system by boundary-value collocation, for any ``kappa`` and any
    winding ``n``.

    Parameters
    ----------
    kappa : float
        GL parameter lambda/xi.  ``kappa = 1`` is the BPS self-dual point
        (m_s = m_v) and reproduces :func:`solve_bps_vortex`; ``kappa < 1`` is
        type-I, ``kappa > 1`` type-II (the screening-current sheath sits
        further from the core).
    n : int
        Vortex winding number (>= 1).  The collocation solver is stable for
        n >= 2, unlike the BPS shooting.
    rho_max : float
        Outer radius of the returned grid, in healing lengths xi (floored
        internally so the gauge field saturates even at large kappa).
    n_points : int
        Number of grid points used by the solver.

    Returns
    -------
    GLVortexProfile
        Radial profile with ``rho`` in healing lengths.
    """
    if kappa <= 0:
        raise ValueError("GL parameter kappa must be > 0.")
    if n < 1:
        raise ValueError("Vortex winding n must be >= 1.")
    return _solve_gl_cached(float(kappa), int(n), float(rho_max), int(n_points))
