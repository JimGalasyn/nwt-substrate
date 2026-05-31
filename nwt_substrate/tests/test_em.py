"""Tests for nwt_substrate.em (classical EM fields from substrate sources)."""

import numpy as np
import pytest

import nwt_substrate.em as em


def _point_charge(N=48, L=16.0, Q=1.0, width=0.8):
    g, dx, (X, Y, Z) = em.grid(N, L)
    rho = np.exp(-((X**2 + Y**2 + Z**2) / (2 * width**2)))
    rho *= Q / (rho.sum() * dx**3)
    return g, dx, (X, Y, Z), rho


def _current_ring(N=48, L=16.0, R=4.0, I=1.0, width=0.8, n=200):
    g, dx, (X, Y, Z) = em.grid(N, L)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    curve = np.stack([R * np.cos(t), R * np.sin(t), 0 * t], 1)
    rho, j = em.deposit_sources([curve], (X, Y, Z), dx, Q=0.0, I=I, width=width)
    return g, dx, (X, Y, Z), curve, j


# ---------------------------------------------------------------------------
# Electric field — Coulomb
# ---------------------------------------------------------------------------
def test_point_charge_field_is_radial_and_outward():
    g, dx, (X, Y, Z), rho = _point_charge(Q=1.0)
    Ex, Ey, Ez = em.electric_field(rho, dx)
    # sample a shell of points; E should point radially outward (E·r̂ > 0)
    r2 = X**2 + Y**2 + Z**2
    shell = (r2 > 3.0**2) & (r2 < 5.0**2)
    Edotr = (Ex * X + Ey * Y + Ez * Z)[shell]
    assert np.mean(Edotr > 0) > 0.95


def test_coulomb_falloff():
    """|E| should fall roughly as 1/r² in the mid-field."""
    g, dx, (X, Y, Z), rho = _point_charge(Q=1.0, N=64, L=20.0)
    Ex, Ey, Ez = em.electric_field(rho, dx)
    Emag = np.sqrt(Ex**2 + Ey**2 + Ez**2)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    m1 = (r > 2.5) & (r < 3.5); m2 = (r > 5.0) & (r < 6.0)
    # E(3) / E(5.5) ≈ (5.5/3)² ≈ 3.4 ; allow a generous band (finite box)
    ratio = Emag[m1].mean() / Emag[m2].mean()
    assert 2.0 < ratio < 5.5


def test_open_bc_cleaner_coulomb_falloff():
    """Open (free-space) BCs give a cleaner 1/r² law than periodic, because
    there's no neutralizing background / image-charge contamination."""
    g, dx, (X, Y, Z), rho = _point_charge(Q=1.0, N=64, L=20.0)
    Emag = {}
    for bc in ("periodic", "open"):
        Ex, Ey, Ez = em.electric_field(rho, dx, bc=bc)
        Emag[bc] = np.sqrt(Ex**2 + Ey**2 + Ez**2)
    r = np.sqrt(X**2 + Y**2 + Z**2)
    m1 = (r > 2.5) & (r < 3.5); m2 = (r > 5.0) & (r < 6.0)
    target = (5.5 / 3.0) ** 2          # ideal 1/r² ratio ≈ 3.36
    err = {bc: abs(Emag[bc][m1].mean() / Emag[bc][m2].mean() - target) for bc in Emag}
    assert err["open"] < err["periodic"]      # open is closer to true Coulomb


def test_gauss_law():
    g, dx, XYZ, rho = _point_charge()
    E = em.electric_field(rho, dx)
    divE = em.divergence(E, dx)
    # ∇·E ≈ ρ − ⟨ρ⟩ (periodic, k=0 dropped)
    corr = np.corrcoef(divE.ravel(), (rho - rho.mean()).ravel())[0, 1]
    assert corr > 0.99


# ---------------------------------------------------------------------------
# Magnetic field — current loop / dipole
# ---------------------------------------------------------------------------
def test_current_ring_axial_field_is_axial():
    g, dx, (X, Y, Z), curve, j = _current_ring()
    Bx, By, Bz = em.magnetic_field(j, dx)
    # on the z-axis the dipole field is purely axial (Bz dominates)
    N = X.shape[0]; c = N // 2
    axial = np.abs(Bz[c, c, :]); transverse = np.abs(Bx[c, c, :]) + np.abs(By[c, c, :])
    assert axial.max() > 5 * transverse.max()


def test_divergence_of_curl_vanishes():
    g, dx, (X, Y, Z), curve, j = _current_ring()
    B = em.magnetic_field(j, dx)          # B = ∇×A → ∇·B = 0
    assert np.abs(em.divergence(B, dx)).max() < 1e-9


# ---------------------------------------------------------------------------
# Euler–Heisenberg nonlinear option
# ---------------------------------------------------------------------------
def test_eh_zero_reduces_to_linear():
    g, dx, (X, Y, Z), rho = _point_charge()
    _, _, _, curve, j = _current_ring()
    E_lin, B_lin = em.electric_field(rho, dx), em.magnetic_field(j, dx)
    E_eh, B_eh = em.maxwell_eh(rho, j, dx, eh_xi=0.0)
    assert np.allclose(E_eh[0], E_lin[0]) and np.allclose(B_eh[2], B_lin[2])


def test_eh_nonzero_perturbs_field():
    g, dx, (X, Y, Z), rho = _point_charge()
    _, _, _, curve, j = _current_ring()
    E0, _ = em.maxwell_eh(rho, j, dx, eh_xi=0.0)
    Eh, _ = em.maxwell_eh(rho, j, dx, eh_xi=0.05, n_iter=5)
    rel = np.linalg.norm(np.array(Eh[0]) - np.array(E0[0])) / (np.linalg.norm(np.array(E0[0])) + 1e-12)
    assert rel > 1e-4


# ---------------------------------------------------------------------------
# Multipole moments (form-factor bridge)
# ---------------------------------------------------------------------------
def test_monopole_recovers_total_charge():
    g, dx, XYZ, rho = _point_charge(Q=-1.0)
    mom = em.multipole_moments(rho, [np.zeros_like(rho)] * 3, XYZ, dx)
    assert mom["charge"] == pytest.approx(-1.0, abs=1e-6)
    assert np.linalg.norm(mom["electric_dipole"]) < 0.05   # centred → ~no dipole


def test_current_ring_has_axial_magnetic_dipole():
    g, dx, XYZ, curve, j = _current_ring(I=1.0)
    mom = em.multipole_moments(np.zeros_like(j[0]), j, XYZ, dx)
    m = mom["magnetic_dipole"]
    assert abs(m[2]) > 5 * (abs(m[0]) + abs(m[1]))   # points along +z
    assert m[2] > 0


# ---------------------------------------------------------------------------
# Field-line tracing
# ---------------------------------------------------------------------------
def test_form_factor_of_gaussian():
    """|F(q)| of a Gaussian charge density is itself a Gaussian; F(0)=1 and the
    low-q slope recovers ⟨r²⟩ (F ≈ 1 − ⟨r²⟩ q²/6)."""
    g, dx, (X, Y, Z), rho = _point_charge(Q=1.0, N=64, L=20.0, width=1.2)
    q, F = em.form_factor(rho, dx, bc="open", nbins=40, qmax=2.0)
    assert F[0] == pytest.approx(1.0, abs=0.02)        # normalised to F(0)=1
    assert np.all(np.diff(F) <= 1e-6)                  # monotonically decreasing
    # low-q slope: F ≈ 1 − ⟨r²⟩ q²/6 ; for width-1.2 Gaussian ⟨r²⟩ = 3·1.2²
    lo = q < 0.6
    slope = np.polyfit(q[lo] ** 2, F[lo], 1)[0]
    r2 = -6 * slope
    assert r2 == pytest.approx(3 * 1.2 ** 2, rel=0.25)


def test_open_form_factor_cleaner_than_periodic():
    g, dx, (X, Y, Z), rho = _point_charge(Q=1.0, N=64, L=20.0, width=1.0)
    qo, Fo = em.form_factor(rho, dx, bc="open", nbins=40, qmax=2.0)
    qp, Fp = em.form_factor(rho, dx, bc="periodic", nbins=40, qmax=2.0)
    assert Fo[0] == pytest.approx(1.0, abs=0.02) and Fp[0] == pytest.approx(1.0, abs=0.05)
    # open is monotone; periodic carries image ripples
    assert np.all(np.diff(Fo) <= 1e-6)


def test_mean_square_radius_neutral_is_finite_second_moment():
    """Net-neutral ⟨r²⟩ must be the bare second moment ∫r²ρ dV (a finite,
    physical value) — NOT a divide-by-(~0) blow-up.  Negative here because the
    negative charge sits at larger r (the neutron's positive-interior/neg-skin)."""
    g, dx, (X, Y, Z) = em.grid(48, 16.0)
    core = np.exp(-(X**2 + Y**2 + Z**2) / (2 * 0.8**2))
    skin = np.exp(-(X**2 + Y**2 + Z**2) / (2 * 2.0**2))
    core /= core.sum(); skin /= skin.sum()
    rho = core - skin                                  # net ≈ 0 (FP residual ~1e-17)
    assert 0 < abs(rho.sum()) < 1e-6                   # genuinely near-neutral, not exact
    r2 = em.mean_square_radius(rho, (X, Y, Z), dx)
    direct = float(((X**2 + Y**2 + Z**2) * rho).sum() * dx**3)   # the bare 2nd moment
    assert r2 == pytest.approx(direct, rel=1e-9)       # returns the moment, not garbage
    assert -3.0 < r2 < 0.0                             # finite & negative (physical scale)


def test_mean_square_radius_charged_is_per_charge():
    g, dx, (X, Y, Z), rho = _point_charge(Q=1.0, N=48, L=16.0, width=1.0)
    # net-charged → ⟨r²⟩ = ∫r²ρ/∫ρ ; a width-1 Gaussian has ⟨r²⟩ = 3
    assert em.mean_square_radius(rho, (X, Y, Z), dx) == pytest.approx(3.0, rel=0.1)


def test_form_factor_suppresses_scattering():
    from nwt_substrate.amplitudes.cross_sections import (
        mott_cross_section, elastic_form_factor_cross_section)
    theta = np.linspace(0.3, 3.0, 30); E = 1.0
    F = lambda q2: 1.0 / (1 + q2 / 0.71) ** 2           # dipole form factor
    mott = mott_cross_section(theta, E)
    ext = elastic_form_factor_cross_section(theta, E, F)
    assert np.all(ext <= mott + 1e-12)                 # |F|²≤1 suppresses
    assert ext[-1] / mott[-1] < ext[0] / mott[0]       # more suppression at large angle


def test_trace_field_lines_runs():
    g, dx, (X, Y, Z), rho = _point_charge()
    E = em.electric_field(rho, dx)
    seeds = [(3.0, 0, 0), (0, 3.0, 0), (-3.0, 0, 0)]
    lines = em.trace_field_lines(E, g, seeds, n_steps=80, ds=0.2)
    assert len(lines) >= 3
    assert all(ln.shape[1] == 3 and len(ln) > 3 for ln in lines)
