"""Tests for the heron sidereal-geometry layer.

Focused on the Exp 11 prediction functions, especially the secondary
asymmetric-jet model added 2026-05-11 evening (memory entry
`vortex-vision-hpa-structural-bh-system.md`).
"""
import numpy as np
import pytest

import nwt_substrate.heron as heron


# Three slots covering one sidereal day for ibm_marrakesh's host.
T_A = 1.7e9
T_B = T_A + 86164.0905 / 2.0   # half sidereal day later
T_C = T_A + 86164.0905         # one sidereal day later


@pytest.mark.skipif(
    not heron.HAS_ASTROPY,
    reason="sidereal-geometry tests require astropy (`pip install astropy`)",
)
class TestAsymmetricJetPrediction:
    """The vortex-vision 2026-05-11 structural model adds a signed
    dipolar correction to the symmetric primary prediction. ε=0 must
    recover the primary exactly; ε≠0 must produce a sign-sensitive
    deviation; |ε|>1 must raise to keep C_v non-negative."""

    def test_epsilon_zero_recovers_symmetric(self):
        """epsilon=0 must reproduce predicted_sigma_pattern bit-for-bit.
        This is the pre-registration safety guard: the asymmetric
        function must not silently change the primary prediction."""
        sym = heron.predicted_sigma_pattern(T_A, T_B, T_C, heron.YORKTOWN)
        asym = heron.predicted_sigma_pattern_asymmetric(
            T_A, T_B, T_C, heron.YORKTOWN, epsilon=0.0
        )
        assert np.allclose(sym["sigma_v_pred"], asym["sigma_v_pred"],
                            atol=1e-15)
        assert np.allclose(sym["C_A"], asym["C_A"], atol=1e-15)
        assert np.allclose(sym["C_B"], asym["C_B"], atol=1e-15)
        assert np.allclose(sym["C_C"], asym["C_C"], atol=1e-15)

    def test_epsilon_positive_differs_from_symmetric(self):
        """Nonzero epsilon must produce a different prediction —
        otherwise the asymmetric model is degenerate."""
        sym = heron.predicted_sigma_pattern(T_A, T_B, T_C, heron.YORKTOWN)
        asym = heron.predicted_sigma_pattern_asymmetric(
            T_A, T_B, T_C, heron.YORKTOWN, epsilon=+0.3
        )
        assert not np.allclose(sym["sigma_v_pred"], asym["sigma_v_pred"])

    def test_sign_sensitivity(self):
        """epsilon=+0.3 and epsilon=-0.3 must produce DIFFERENT
        predictions (not just same-magnitude). This is what makes the
        function a useful test of the one-sided-jet hypothesis."""
        asym_pos = heron.predicted_sigma_pattern_asymmetric(
            T_A, T_B, T_C, heron.YORKTOWN, epsilon=+0.3
        )
        asym_neg = heron.predicted_sigma_pattern_asymmetric(
            T_A, T_B, T_C, heron.YORKTOWN, epsilon=-0.3
        )
        assert not np.allclose(asym_pos["sigma_v_pred"],
                                asym_neg["sigma_v_pred"])

    def test_output_shape(self):
        """sigma_v_pred is a 7-vector (one per K_7 vertex)."""
        result = heron.predicted_sigma_pattern_asymmetric(
            T_A, T_B, T_C, heron.YORKTOWN, epsilon=0.2
        )
        assert result["sigma_v_pred"].shape == (7,)
        assert result["C_A"].shape == (7,)
        assert result["C_B"].shape == (7,)
        assert result["C_C"].shape == (7,)

    def test_metadata_echoed(self):
        """Result dict carries epsilon and coupling_model for
        downstream model-comparison bookkeeping."""
        result = heron.predicted_sigma_pattern_asymmetric(
            T_A, T_B, T_C, heron.YORKTOWN, epsilon=0.42
        )
        assert result["epsilon"] == 0.42
        assert result["coupling_model"] == "asymmetric_jet"

    def test_epsilon_too_large_rejected(self):
        """|epsilon|>1 would make C_v = (1+ε·cos) cos² negative for
        some cos values, breaking the coupling interpretation."""
        with pytest.raises(ValueError, match="epsilon must lie in"):
            heron.predicted_sigma_pattern_asymmetric(
                T_A, T_B, T_C, heron.YORKTOWN, epsilon=1.5
            )
        with pytest.raises(ValueError, match="epsilon must lie in"):
            heron.predicted_sigma_pattern_asymmetric(
                T_A, T_B, T_C, heron.YORKTOWN, epsilon=-1.5
            )

    def test_epsilon_boundary_accepted(self):
        """epsilon = ±1 (the edge of validity) must still work."""
        for eps in (-1.0, +1.0):
            result = heron.predicted_sigma_pattern_asymmetric(
                T_A, T_B, T_C, heron.YORKTOWN, epsilon=eps
            )
            assert result["sigma_v_pred"].shape == (7,)
            assert np.all(np.isfinite(result["sigma_v_pred"]))
