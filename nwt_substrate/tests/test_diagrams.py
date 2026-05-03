"""Tests for nwt_substrate.diagrams."""

import matplotlib
matplotlib.use("Agg")  # headless

import pytest
import numpy as np
from matplotlib.figure import Figure

import nwt_substrate as nwt
from nwt_substrate import diagrams as dg


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------

def test_torus_xyz_returns_correct_shapes():
    """Scalar input -> 3 scalars; array input -> 3 arrays of same shape."""
    x, y, z = dg.torus_xyz(0.0, 0.0)
    assert isinstance(x, float) and isinstance(y, float) and isinstance(z, float)

    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, 2 * np.pi, 50)
    x, y, z = dg.torus_xyz(u, v)
    assert x.shape == y.shape == z.shape == (50,)


def test_torus_xyz_geometry_invariant():
    """Distance from torus center axis = R + r cos(v)."""
    R, r = 1.5, 0.55
    for v in np.linspace(0, 2 * np.pi, 13)[:-1]:
        x, y, z = dg.torus_xyz(0.0, v, R=R, r=r)
        rho = np.sqrt(x**2 + y**2)
        assert abs(rho - (R + r * np.cos(v))) < 1e-12


# ---------------------------------------------------------------------------
# K_7 combinatorial helpers
# ---------------------------------------------------------------------------

def test_K7_eulerian_circuit_has_21_edges():
    edges = dg.K7_eulerian_circuit()
    assert len(edges) == 21


def test_K7_eulerian_circuit_uses_all_edges_once():
    edges = dg.K7_eulerian_circuit()
    seen = set(frozenset(e) for e in edges)
    assert len(seen) == 21  # all 21 unique edges traversed exactly once


def test_K7_eulerian_circuit_is_closed():
    edges = dg.K7_eulerian_circuit(start=0)
    # walk closes back to the start vertex
    assert edges[-1][1] == edges[0][0]


def test_K7_heffter_triangles_count():
    """Heffter embedding has F = 14 triangular faces (V - E + F = 0)."""
    tris = dg.K7_heffter_triangles()
    assert len(tris) == 14
    for t in tris:
        assert len(set(t)) == 3
        for v in t:
            assert 0 <= v < 7


# ---------------------------------------------------------------------------
# Figure factories return matplotlib Figures
# ---------------------------------------------------------------------------

def test_figure_torus_knot_returns_figure():
    fig = dg.figure_torus_knot(p=2, q=3)
    assert isinstance(fig, Figure)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_figure_torus_knot_accepts_arbitrary_pq():
    """Cinquefoil (2,5), Wilson (3,2), unknot (1,7) all render."""
    import matplotlib.pyplot as plt
    for p, q in [(2, 5), (3, 2), (1, 7), (4, 3)]:
        fig = dg.figure_torus_knot(p=p, q=q)
        assert isinstance(fig, Figure)
        plt.close(fig)


def test_figure_K7_on_torus_returns_figure():
    fig = dg.figure_K7_on_torus()
    assert isinstance(fig, Figure)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_figure_unified_side_by_side():
    fig = dg.figure_unified(layout="side-by-side")
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_figure_unified_single_panel():
    fig = dg.figure_unified(layout="single")
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 1
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_figure_K7_traversals_modes():
    import matplotlib.pyplot as plt
    for mode in ("eulerian", "triangles", "cycles"):
        fig = dg.figure_K7_traversals(mode=mode)
        assert isinstance(fig, Figure)
        plt.close(fig)


def test_figure_K7_traversals_rejects_unknown_mode():
    with pytest.raises(ValueError):
        dg.figure_K7_traversals(mode="kruskal")


def test_figure_paper18_unified_returns_figure():
    """The exact Paper 18 §4.4 figure factory."""
    fig = dg.figure_paper18_unified()
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 2
    import matplotlib.pyplot as plt
    plt.close(fig)


# ---------------------------------------------------------------------------
# Top-level access via nwt.diagrams
# ---------------------------------------------------------------------------

def test_diagrams_exposed_at_top_level():
    """nwt.diagrams should be the same module."""
    assert nwt.diagrams is dg


def test_save_to_writes_file(tmp_path):
    """save_to=path argument writes the file before returning."""
    out = tmp_path / "torus.png"
    fig = dg.figure_torus_knot(p=2, q=3, save_to=str(out))
    assert out.exists() and out.stat().st_size > 0
    import matplotlib.pyplot as plt
    plt.close(fig)
