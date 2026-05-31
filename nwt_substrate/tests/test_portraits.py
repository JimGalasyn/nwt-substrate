"""Tests for nwt_substrate.portraits (particle portrait engine)."""

import matplotlib
matplotlib.use("Agg")  # headless

import numpy as np
import pytest
from matplotlib.figure import Figure

import nwt_substrate.portraits as pp


# ---------------------------------------------------------------------------
# Field construction
# ---------------------------------------------------------------------------

def test_field_shapes_and_amplitude_range():
    pf = pp.build_particle_field(2, 3, N=48)
    assert pf.N == 48
    assert pf.phase.shape == pf.glow.shape == pf.amp.shape == (48, 48, 48)
    # |psi| is a probability-like amplitude in [0, 1]
    assert pf.amp.min() >= -1e-9
    assert pf.amp.max() <= 1.0 + 1e-9


def test_amplitude_vanishes_on_core_and_saturates_in_bulk():
    """|psi| -> 0 somewhere (the core) and -> 1 in the far corners (bulk)."""
    pf = pp.build_particle_field(2, 3, N=64)
    assert pf.amp.min() < 0.2          # core is resolved
    # the 8 grid corners are far from the knot: bulk condensate
    corners = pf.amp[[0, -1]][:, [0, -1]][:, :, [0, -1]]
    assert corners.min() > 0.9


def test_glow_peaks_where_amplitude_is_low():
    """Emission is concentrated on the vortex core (low |psi|)."""
    pf = pp.build_particle_field(2, 3, N=64)
    brightest = np.unravel_index(np.argmax(pf.glow), pf.glow.shape)
    assert pf.amp[brightest] < 0.3


def test_hopf_field_builds():
    pf = pp.build_particle_field(hopf=True, N=48)
    assert pf.meta["kind"] == "Hopf link"
    assert pf.amp.shape == (48, 48, 48)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def test_render_returns_rgb_image_in_unit_range():
    pf = pp.build_particle_field(2, 1, N=48)
    img = pp.render_field(pf)
    assert img.shape == (48, 48, 3)
    assert img.min() >= 0.0 and img.max() <= 1.0


def test_portrait_returns_figure():
    fig = pp.portrait(2, 3, N=40, title="trefoil")
    assert isinstance(fig, Figure)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_portrait_save_to_writes_file(tmp_path):
    out = tmp_path / "p.png"
    fig = pp.portrait(2, 1, N=40, save_to=str(out))
    assert out.exists() and out.stat().st_size > 0
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_gallery_returns_figure():
    specs = type(pp.CARRIERS)([
        ("a", dict(p=2, q=1)),
        ("b", dict(p=2, q=3)),
    ])
    fig = pp.gallery(specs, N=36, ncols=2)
    assert isinstance(fig, Figure)
    import matplotlib.pyplot as plt
    plt.close(fig)


# ---------------------------------------------------------------------------
# Carrier mapping
# ---------------------------------------------------------------------------

def test_n_q_to_knot_mapping():
    assert pp.n_q_to_knot(1) == {"p": 2, "q": 1}
    assert pp.n_q_to_knot(2) == {"hopf": True}
    assert pp.n_q_to_knot(3) == {"p": 2, "q": 3}
    assert pp.n_q_to_knot(5) == {"p": 2, "q": 5}


def test_n_q_to_knot_rejects_nontoroidal():
    with pytest.raises(ValueError):
        pp.n_q_to_knot(4)  # figure-eight is non-toroidal
