"""Tests for nwt_substrate.solitons.faddeev.

Analytic / closed-form assertions only -- no JAX, no GPU, no engine import (the
dependency arrow is one-way: engine -> library).  Covers the CP^1 map, the VOS
solid-angle area form, the rational-map hopfion seed (charge, compactness,
smoothness, Q_H = n*m), the Faddeev energy identity, and the Whitehead charge.
The numeric primitives here also serve as the reference oracle the JAX engine
is cross-checked against (that cross-test lives engine-side, where JAX exists).
"""
from __future__ import annotations

import numpy as np
import pytest

from nwt_substrate.solitons import BoxGrid
from nwt_substrate.solitons.faddeev import (
    n_from_Z,
    solid_angle_triangle,
    plaquette_area_form,
    rational_hopfion,
    faddeev_energy,
    energy_density,
    soliton_size,
    whitehead_hopf_charge,
)


# ---------------------------------------------------------------------------
# solid_angle_triangle (Van Oosterom-Strackee)
# ---------------------------------------------------------------------------
def test_solid_angle_octant_is_half_pi():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    c = np.array([0.0, 0.0, 1.0])
    assert solid_angle_triangle(a, b, c) == pytest.approx(np.pi / 2)


def test_solid_angle_degenerate_is_zero():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    # collinear (a, b, a): the triangle is degenerate -> zero area
    assert solid_angle_triangle(a, b, a) == pytest.approx(0.0)


def test_solid_angle_orientation_flips_sign():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    c = np.array([0.0, 0.0, 1.0])
    assert solid_angle_triangle(a, c, b) == pytest.approx(
        -solid_angle_triangle(a, b, c))


# ---------------------------------------------------------------------------
# n_from_Z
# ---------------------------------------------------------------------------
def test_n_from_Z_is_unit():
    rng = np.random.default_rng(0)
    Z1 = rng.normal(size=(4, 4, 4)) + 1j * rng.normal(size=(4, 4, 4))
    Z2 = rng.normal(size=(4, 4, 4)) + 1j * rng.normal(size=(4, 4, 4))
    n1, n2, n3 = n_from_Z(Z1, Z2)
    assert np.allclose(n1**2 + n2**2 + n3**2, 1.0)


def test_n_from_Z_vacuum_spinor():
    # (Z1, Z2) = (0, 1) is the south pole n = (0, 0, -1)
    n1, n2, n3 = n_from_Z(np.array([0.0 + 0j]), np.array([1.0 + 0j]))
    assert (n1[0], n2[0], n3[0]) == pytest.approx((0.0, 0.0, -1.0))


# ---------------------------------------------------------------------------
# rational_hopfion seed
# ---------------------------------------------------------------------------
def test_rational_hopfion_is_unit_spinor():
    grid = BoxGrid(N=48, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)
    assert np.allclose(np.abs(Z1) ** 2 + np.abs(Z2) ** 2, 1.0)


def test_rational_hopfion_is_compact_vacuum_at_corner():
    grid = BoxGrid(N=48, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)
    n1, n2, n3 = n_from_Z(Z1, Z2)
    # box corner is deep vacuum: south pole n3 = -1 (periodic-safe)
    assert n3[0, 0, 0] == pytest.approx(-1.0, abs=1e-9)


def test_rational_hopfion_charge_is_one():
    grid = BoxGrid(N=64, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)
    n1, n2, n3 = n_from_Z(Z1, Z2)
    assert whitehead_hopf_charge(n1, n2, n3, grid) == pytest.approx(1.0, abs=3e-3)


def test_rational_hopfion_charge_nm():
    # Q_H = n * m
    grid = BoxGrid(N=64, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5, n=1, m=2)
    n1, n2, n3 = n_from_Z(Z1, Z2)
    assert whitehead_hopf_charge(n1, n2, n3, grid) == pytest.approx(2.0, abs=0.1)


def test_rational_hopfion_is_smooth():
    # a singular seed would show a per-step jump near 2 (antipodal); a smooth
    # field jumps by O(dx) -- well under 0.5 at this resolution.
    grid = BoxGrid(N=64, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)
    for c in n_from_Z(Z1, Z2):
        jump = max(np.abs(np.roll(c, -1, ax) - c).max() for ax in range(3))
        assert jump < 0.5


def test_rational_hopfion_requires_w_le_R():
    grid = BoxGrid(N=32, L=16.0)
    with pytest.raises(ValueError):
        rational_hopfion(grid, R=3.0, w=4.0)


# ---------------------------------------------------------------------------
# energy
# ---------------------------------------------------------------------------
def test_vacuum_energy_is_zero():
    grid = BoxGrid(N=16, L=16.0)
    shape = (grid.N, grid.N, grid.N)
    n1 = np.zeros(shape)
    n2 = np.zeros(shape)
    n3 = -np.ones(shape)          # uniform south-pole vacuum
    assert faddeev_energy(n1, n2, n3, grid.dx, c4=4.0) == pytest.approx(0.0)


def test_energy_density_sums_to_energy():
    grid = BoxGrid(N=48, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)
    n1, n2, n3 = n_from_Z(Z1, Z2)
    dens = energy_density(n1, n2, n3, grid.dx, c4=4.0)
    total = faddeev_energy(n1, n2, n3, grid.dx, c4=4.0)
    assert np.sum(dens) * grid.dx**3 == pytest.approx(total)


def test_soliton_size_split_matches_energy():
    grid = BoxGrid(N=48, L=16.0)
    Z1, Z2 = rational_hopfion(grid, R=3.5)
    n1, n2, n3 = n_from_Z(Z1, Z2)
    sz = soliton_size(n1, n2, n3, grid, c4=4.0)
    assert sz["E2"] + sz["E4"] == pytest.approx(
        faddeev_energy(n1, n2, n3, grid.dx, c4=4.0))
    assert sz["R_E"] > 0.0


# ---------------------------------------------------------------------------
# whitehead_hopf_charge
# ---------------------------------------------------------------------------
def test_vacuum_charge_is_zero():
    grid = BoxGrid(N=16, L=16.0)
    shape = (grid.N, grid.N, grid.N)
    n1 = np.zeros(shape)
    n2 = np.zeros(shape)
    n3 = -np.ones(shape)
    assert whitehead_hopf_charge(n1, n2, n3, grid) == pytest.approx(0.0, abs=1e-9)
