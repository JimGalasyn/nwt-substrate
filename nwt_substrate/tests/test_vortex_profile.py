"""Tests for nwt_substrate.condensate.vortex_profile (BPS vortex)."""

import numpy as np
import pytest

from nwt_substrate.condensate import solve_bps_vortex, BPSVortexProfile


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
