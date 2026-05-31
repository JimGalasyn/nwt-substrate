"""Tests for nwt_substrate.condensate.vortex_profile (BPS + general GL vortex)."""

import numpy as np
import pytest

from nwt_substrate.condensate import (
    solve_bps_vortex,
    BPSVortexProfile,
    solve_gl_vortex,
    GLVortexProfile,
)


def test_n1_converges_to_known_slope():
    """The n=1 BPS shooting parameter is c1 ~ 0.6033."""
    pr = solve_bps_vortex(n=1)
    assert isinstance(pr, BPSVortexProfile)
    assert pr.c1 == pytest.approx(0.6033, abs=2e-3)


def test_boundary_conditions():
    """f: 0 -> 1, a: 0 -> n across the domain."""
    pr = solve_bps_vortex(n=1)
    assert pr.f[0] == pytest.approx(0.0, abs=1e-3)
    assert pr.f[-1] == pytest.approx(1.0, abs=1e-3)
    assert pr.a[0] == pytest.approx(0.0, abs=1e-3)
    assert pr.a[-1] == pytest.approx(1.0, abs=2e-3)


@pytest.mark.parametrize("n", [2, 3])
def test_higher_winding_not_supported(n):
    """Multi-winding is documented as unsupported and raises cleanly
    (shooting is unstable for n>=2; needs a stiff BVP solver)."""
    with pytest.raises(NotImplementedError):
        solve_bps_vortex(n=n)


def test_scalar_decays_at_unit_bps_rate():
    """At BPS m_s = m_v = 1/xi, so 1 - f decays at rate ~1 (not sqrt(2))."""
    import numpy as np
    pr = solve_bps_vortex(n=1, rho_max=25.0, n_points=6000)
    m = (pr.rho > 4) & (pr.rho < 8) & (pr.f < 0.9999) & (pr.f > 0)
    # log(1-f) + 0.5 log(rho) has slope -rate for 1-f ~ K0(rate*rho)
    rate = -np.polyfit(pr.rho[m], np.log(1 - pr.f[m]) + 0.5 * np.log(pr.rho[m]), 1)[0]
    assert rate == pytest.approx(1.0, abs=0.05)


def test_scalar_is_monotone_increasing():
    """f rises essentially monotonically from core to bulk."""
    pr = solve_bps_vortex(n=1)
    # allow tiny numerical noise at the saturated tail
    assert np.min(np.diff(pr.f)) > -1e-4
    assert pr.f[-1] > pr.f[len(pr.f) // 2] > pr.f[0]


def test_f_at_clamps_outside_grid():
    """f_at is 0 at the core and clamps to 1 beyond the solved grid."""
    pr = solve_bps_vortex(n=1)
    assert pr.f_at(0.0) == pytest.approx(0.0, abs=1e-3)
    assert pr.f_at(1e6) == pytest.approx(1.0)
    vals = pr.f_at(np.array([0.0, 1.0, 3.0, 1e6]))
    assert vals.shape == (4,)
    assert np.all(np.diff(vals) >= 0)


def test_solver_is_cached():
    """Repeated calls with identical args return the same cached object."""
    a = solve_bps_vortex(n=1)
    b = solve_bps_vortex(n=1)
    assert a is b


def test_invalid_winding_rejected():
    with pytest.raises(ValueError):
        solve_bps_vortex(n=0)


# ---------------------------------------------------------------------------
# General gauged Ginzburg-Landau vortex (any kappa = lambda/xi)
# ---------------------------------------------------------------------------
def test_gl_kappa1_reproduces_bps():
    """At the BPS self-dual point (kappa=1, m_s=m_v) the second-order GL
    solution must coincide with the first-order BPS profile."""
    bps = solve_bps_vortex(n=1)
    gl = solve_gl_vortex(kappa=1.0, n=1)
    assert isinstance(gl, GLVortexProfile)
    mask = (bps.rho > 0.2) & (bps.rho < 8.0)
    f_gl = gl.f_at(bps.rho[mask])
    a_gl = gl.a_at(bps.rho[mask])
    assert np.max(np.abs(f_gl - bps.f[mask])) < 5e-3
    assert np.max(np.abs(a_gl - bps.a[mask])) < 5e-3


def test_gl_boundary_conditions():
    """f: 0 -> 1, a: 0 -> n for a type-II vortex."""
    gl = solve_gl_vortex(kappa=1.5, n=1)
    assert gl.f[0] == pytest.approx(0.0, abs=2e-3)
    assert gl.f[-1] == pytest.approx(1.0, abs=2e-3)
    assert gl.a[0] == pytest.approx(0.0, abs=2e-3)
    assert gl.a[-1] == pytest.approx(1.0, abs=3e-3)


@pytest.mark.parametrize("n", [2, 3])
def test_gl_solves_higher_winding(n):
    """The collocation solver handles n>=2 (where BPS shooting is unstable)."""
    gl = solve_gl_vortex(kappa=1.0, n=n)
    assert gl.f[-1] == pytest.approx(1.0, abs=3e-3)
    assert gl.a[-1] == pytest.approx(float(n), abs=5e-3)


def test_gl_supercurrent_sheath_moves_outward_with_kappa():
    """The screening-current sheath f^2(n-a)/rho sits near the core at BPS and
    migrates outward as kappa increases into the type-II regime."""
    peaks = {}
    for kappa in (1.0, 1.5, 2.0):
        gl = solve_gl_vortex(kappa=kappa, n=1)
        J = gl.supercurrent()
        peaks[kappa] = gl.rho[np.argmax(np.clip(J, 0, None))]
    assert peaks[1.0] < peaks[1.5] < peaks[2.0]


def test_gl_invalid_params_rejected():
    with pytest.raises(ValueError):
        solve_gl_vortex(kappa=0.0)
    with pytest.raises(ValueError):
        solve_gl_vortex(kappa=1.0, n=0)


def test_gl_solver_is_cached():
    a = solve_gl_vortex(kappa=1.3, n=1)
    b = solve_gl_vortex(kappa=1.3, n=1)
    assert a is b


def test_gl_supercurrent_at_explicit_radii():
    """supercurrent(r=...) evaluates the screening current off-grid (0 at the
    core, decaying in the bulk) — covers the non-None branch."""
    gl = solve_gl_vortex(kappa=1.2, n=1)
    r = np.array([0.0, 1.0, 3.0, 1e6])
    J = gl.supercurrent(r)
    assert J.shape == r.shape
    assert J[0] == pytest.approx(0.0, abs=1e-6)    # f^2 -> 0 at the core
    assert J[-1] == pytest.approx(0.0, abs=1e-6)   # (n - a) -> 0 in the bulk
    assert J[1] > 0 and J[2] > 0                    # nonzero sheath in between


def test_gl_rhs_gauge_equation_sign():
    """gl_vortex_rhs returns a'' = a'/rho - f^2(n-a), i.e. the homogeneous gauge
    equation is a'' - a'/rho + f^2(n-a) = 0 (the documented + sign on the source
    term, not -).  Pins the sign the PR review caught in the comment block."""
    from nwt_substrate.condensate import gl_vortex_rhs
    rho = np.array([2.0])
    y = np.array([[0.5], [0.1], [0.3], [0.2]])      # f, f', a, a'
    out = gl_vortex_rhs(rho, y, kappa=1.5, n=1)
    f, fp, a, ap, n = 0.5, 0.1, 0.3, 0.2, 1
    assert out[3, 0] == pytest.approx(ap / 2.0 - f**2 * (n - a))     # a'' = a'/rho - f^2(n-a)
    assert out[3, 0] != pytest.approx(ap / 2.0 + f**2 * (n - a))     # NOT the wrong sign
    kc = 1.5 / np.sqrt(2.0)                          # scalar-eq sign too
    assert out[1, 0] == pytest.approx(-fp / 2.0 + (n - a) ** 2 * f / 4.0 - kc**2 * (1 - f**2) * f)
