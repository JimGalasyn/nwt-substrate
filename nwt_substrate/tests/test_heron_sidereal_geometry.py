"""Extra coverage for heron.sidereal_geometry beyond the asymmetric-jet tests.

Covers:
  - Observatory dataclass
  - default_k7_lab_positions (pure math)
  - directional_match_score (pure numpy)
  - lab_to_icrs_matrix / lab_to_icrs (orthogonality, single-time consistency)
  - lst_hours / next_lst_match_unix / schedule_triplet_at_lst (astropy-gated)

Astropy-gated tests skip cleanly when astropy is unavailable, matching the
pattern in test_heron_sidereal.py.
"""

from __future__ import annotations

import numpy as np
import pytest

import nwt_substrate.heron as heron
from nwt_substrate.heron.sidereal_geometry import (
    EHNINGEN,
    YORKTOWN,
    Observatory,
    default_k7_lab_positions,
    directional_match_score,
)


# ---------------------------------------------------------------
# Observatory + module-level singletons
# ---------------------------------------------------------------

def test_observatory_is_frozen_dataclass():
    obs = Observatory(name="test", lat_deg=10.0, lon_deg=20.0, height_m=5.0)
    with pytest.raises(Exception):
        obs.name = "mutated"  # frozen → AttributeError or FrozenInstanceError


def test_yorktown_singleton_present_with_negative_longitude():
    assert YORKTOWN.name.startswith("IBM Yorktown")
    assert YORKTOWN.lat_deg == pytest.approx(41.2089, abs=1e-3)
    assert YORKTOWN.lon_deg < 0  # west = negative


def test_ehningen_singleton_present_with_positive_longitude():
    assert EHNINGEN.name == "IBM Ehningen"
    assert EHNINGEN.lon_deg > 0  # east = positive


# ---------------------------------------------------------------
# default_k7_lab_positions (pure)
# ---------------------------------------------------------------

def test_default_k7_lab_positions_shape_and_unit_norm():
    pts = default_k7_lab_positions()
    assert pts.shape == (7, 3)
    norms = np.linalg.norm(pts, axis=1)
    np.testing.assert_allclose(norms, 1.0, atol=1e-15)


def test_default_k7_lab_positions_center_near_zenith():
    """Vertex 0 (center qubit) is placed near the zenith (z-component close to 1)."""
    pts = default_k7_lab_positions()
    assert pts[0, 2] > 0.999


def test_default_k7_lab_positions_hexagon_symmetry():
    """Outer 6 qubits sum to ~0 in chip plane (regular-hexagon symmetry)."""
    pts = default_k7_lab_positions()
    outer_xy = pts[1:, :2]
    np.testing.assert_allclose(outer_xy.sum(axis=0), 0.0, atol=1e-12)


def test_default_k7_lab_positions_radius_argument_independence_of_z():
    """Different radii leave the chip plane normal direction intact."""
    p1 = default_k7_lab_positions(radius_m=0.001)
    p2 = default_k7_lab_positions(radius_m=0.01)
    # All vertices remain unit vectors (normalised post-construction)
    np.testing.assert_allclose(np.linalg.norm(p1, axis=1), 1.0, atol=1e-15)
    np.testing.assert_allclose(np.linalg.norm(p2, axis=1), 1.0, atol=1e-15)


# ---------------------------------------------------------------
# directional_match_score (pure numpy)
# ---------------------------------------------------------------

def test_directional_match_perfect_self_match():
    """obs == pred → pearson_r = 1, cos_sim = 1, chi2 = 0."""
    pred = np.array([1.0, -0.5, 0.3, 0.0, 0.2, -0.1, 0.7])
    out = directional_match_score(pred, pred)
    assert out["pearson_r"] == pytest.approx(1.0, abs=1e-12)
    assert out["cosine_similarity"] == pytest.approx(1.0, abs=1e-12)
    assert out["best_amplitude_chi2"] == pytest.approx(0.0, abs=1e-12)
    assert out["dof"] == 6


def test_directional_match_antiphase_negative_pearson_negative_cos():
    pred = np.array([1.0, -0.5, 0.3, 0.0, 0.2, -0.1, 0.7])
    obs = -pred
    out = directional_match_score(obs, pred)
    assert out["pearson_r"] == pytest.approx(-1.0, abs=1e-12)
    assert out["cosine_similarity"] == pytest.approx(-1.0, abs=1e-12)
    assert out["best_amplitude"] < 0


def test_directional_match_orthogonal_zero_cos():
    obs = np.array([1.0, 0.0, 0.0])
    pred = np.array([0.0, 1.0, 0.0])
    out = directional_match_score(obs, pred)
    assert out["cosine_similarity"] == pytest.approx(0.0, abs=1e-15)


def test_directional_match_shape_mismatch_raises():
    with pytest.raises(ValueError, match="same shape"):
        directional_match_score([1.0, 2.0], [1.0, 2.0, 3.0])


def test_directional_match_with_weights():
    """observed_std arg weighs the chi² fit."""
    pred = np.array([1.0, 2.0, 3.0, 4.0])
    obs = pred * 2.0
    out = directional_match_score(obs, pred, observed_std=[0.1, 0.1, 0.1, 0.1])
    # Perfect rescaling: A_hat = 2.0, chi² = 0
    assert out["best_amplitude"] == pytest.approx(2.0, rel=1e-12)
    assert out["best_amplitude_chi2"] == pytest.approx(0.0, abs=1e-12)


def test_directional_match_zero_predicted_safe():
    """All-zero predicted → no division by zero, cos_sim = 0."""
    obs = np.array([1.0, 2.0, 3.0])
    pred = np.zeros(3)
    out = directional_match_score(obs, pred)
    assert out["cosine_similarity"] == 0.0
    assert out["pearson_r"] == 0.0  # zero variance in pred


# ---------------------------------------------------------------
# Astropy-gated functions
# ---------------------------------------------------------------

T_A = 1.7e9
T_B = T_A + 86164.0905 / 2.0
T_C = T_A + 86164.0905


@pytest.mark.skipif(
    not heron.HAS_ASTROPY,
    reason="requires astropy",
)
class TestAstropyDependentSidereal:
    """Exercise the astropy-gated geometry surface end-to-end."""

    def test_lab_to_icrs_matrix_is_near_orthogonal(self):
        """Lab basis → ICRS basis is a near-rotation (small aberration/nutation
        residuals); R^T R ≈ I to ~1e-3, det ≈ ±1."""
        from nwt_substrate.heron.sidereal_geometry import lab_to_icrs_matrix
        R = lab_to_icrs_matrix(T_A, YORKTOWN)
        np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-3)
        assert abs(abs(np.linalg.det(R)) - 1.0) < 1e-3

    def test_lab_to_icrs_unit_vector_stays_unit(self):
        from nwt_substrate.heron.sidereal_geometry import lab_to_icrs
        out = lab_to_icrs(np.array([1.0, 0.0, 0.0]), T_A, YORKTOWN)
        assert np.linalg.norm(out) == pytest.approx(1.0, rel=1e-10)

    def test_predicted_sigma_pattern_default_shapes(self):
        """Symmetric primary returns 7-vector and bookkeeping keys."""
        from nwt_substrate.heron.sidereal_geometry import predicted_sigma_pattern
        out = predicted_sigma_pattern(T_A, T_B, T_C, YORKTOWN)
        assert out["sigma_v_pred"].shape == (7,)
        assert out["C_A"].shape == (7,)
        assert "predicted_axis_name" in out
        assert "predicted_axis_radec" in out

    def test_predicted_sigma_pattern_wrong_shape_raises(self):
        """Custom positions must be (7,3) — anything else rejected."""
        from nwt_substrate.heron.sidereal_geometry import predicted_sigma_pattern
        with pytest.raises(ValueError, match="qubit_lab_positions must be shape"):
            predicted_sigma_pattern(
                T_A, T_B, T_C, YORKTOWN,
                qubit_lab_positions=np.ones((5, 3)),
            )

    def test_predicted_sigma_pattern_asymmetric_wrong_shape_raises(self):
        from nwt_substrate.heron.sidereal_geometry import predicted_sigma_pattern_asymmetric
        with pytest.raises(ValueError, match="qubit_lab_positions must be shape"):
            predicted_sigma_pattern_asymmetric(
                T_A, T_B, T_C, YORKTOWN,
                qubit_lab_positions=np.ones((4, 3)),
                epsilon=0.1,
            )

    def test_lst_hours_in_valid_range(self):
        from nwt_substrate.heron.sidereal_geometry import lst_hours
        h = lst_hours(T_A, YORKTOWN)
        assert 0.0 <= h < 24.0

    def test_lst_hours_advances_with_solar_time(self):
        """LST advances by ~12 sidereal hours over half a sidereal day."""
        from nwt_substrate.heron.sidereal_geometry import lst_hours
        h0 = lst_hours(T_A, YORKTOWN)
        h_half = lst_hours(T_A + 86164.0905 / 2.0, YORKTOWN)
        delta = (h_half - h0) % 24.0
        assert abs(delta - 12.0) < 0.05  # within ~3 sidereal-minutes of half day

    def test_next_lst_match_unix_lands_at_target(self):
        from nwt_substrate.heron.sidereal_geometry import lst_hours, next_lst_match_unix
        target = 6.0  # 6 sidereal hours
        t_match = next_lst_match_unix(target, YORKTOWN, after_unix=T_A,
                                       tol_seconds=2.0)
        h = lst_hours(t_match, YORKTOWN)
        residual = (target - h + 12.0) % 24.0 - 12.0
        # Within 5 sidereal seconds of target after ≤4 iterations
        assert abs(residual) * 3600.0 < 5.0
        assert t_match >= T_A

    def test_schedule_triplet_at_lst_returns_three_slots(self):
        from nwt_substrate.heron.sidereal_geometry import schedule_triplet_at_lst
        plan = schedule_triplet_at_lst(12.0, YORKTOWN, after_unix=T_A)
        assert "A_unix" in plan and "B_unix" in plan and "C_unix" in plan
        assert plan["observatory"] == YORKTOWN.name
        assert plan["target_lst_hours"] == 12.0
        # B is ~half sidereal day after A; C is ~1 sidereal day after A
        assert plan["B_unix"] - plan["A_unix"] == pytest.approx(86164.0905 / 2, rel=1e-9)
        assert plan["C_unix"] - plan["A_unix"] == pytest.approx(86164.0905, rel=1e-9)


# ---------------------------------------------------------------
# Import-error path (no-astropy code branch)
# ---------------------------------------------------------------

def test_lab_to_icrs_matrix_raises_without_astropy(monkeypatch):
    """The ImportError branch is exercised by mocking the import statement."""
    import builtins
    import nwt_substrate.heron.sidereal_geometry as mod

    real_import = builtins.__import__

    def _block_astropy(name, *args, **kwargs):
        if name.startswith("astropy"):
            raise ImportError("blocked")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _block_astropy)
    with pytest.raises(ImportError, match="requires astropy"):
        mod.lab_to_icrs_matrix(T_A, YORKTOWN)
    with pytest.raises(ImportError, match="requires astropy"):
        mod.lst_hours(T_A, YORKTOWN)
