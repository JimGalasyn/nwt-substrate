"""Tests for the dynamical-framing Chern-Simons / CME migration.

Covers:
  * `condensate.chiral_magnetic.cme_growth_rate` -- the closed-form helical-mode
    growth rate ``eta(-k^2 +/- 2 mu5 k)`` (the validated chiral-MHD result),
  * `condensate.chiral_magnetic.cme_current` -- the lattice CME term ``2 mu5 B``,
  * `cosmology.eta_B_sign.eta_B_sign` -- the J_parent -> eta_B sign map,
  * `topology.linking_invariants.net_linking` -- the pseudoscalar net linking.
"""
from __future__ import annotations

import numpy as np
import pytest

from nwt_substrate.condensate.chiral_magnetic import cme_growth_rate, cme_current
from nwt_substrate.cosmology.eta_B_sign import eta_B_sign
from nwt_substrate.diagrams import hopf_link_curves
from nwt_substrate.topology import net_linking, borromean_rings


# ---------------------------------------------------------------------------
# cme_growth_rate -- lambda(k) = eta(-k^2 + helicity * 2 mu5 k)
# ---------------------------------------------------------------------------
class TestCMEGrowthRate:
    def test_positive_helicity_grows_when_mu5_above_half_k(self):
        # mu5 > k/2  =>  unstable (rate > 0) for +helicity
        k, mu5 = 1.0, 0.8
        assert mu5 > k / 2
        assert cme_growth_rate(k, mu5, helicity=+1) > 0

    def test_positive_helicity_decays_when_mu5_below_half_k(self):
        # mu5 < k/2  =>  stable (rate < 0) for +helicity
        k, mu5 = 1.0, 0.3
        assert mu5 < k / 2
        assert cme_growth_rate(k, mu5, helicity=+1) < 0

    def test_negative_helicity_is_the_mirror(self):
        # -helicity sector is the parity mirror: it decays where +helicity grows
        k, mu5 = 1.0, 0.8
        assert cme_growth_rate(k, mu5, helicity=-1) < 0
        # explicit mirror: rate(k, mu5, -1) == rate(k, -mu5, +1)
        assert cme_growth_rate(k, mu5, helicity=-1) == pytest.approx(
            cme_growth_rate(k, -mu5, helicity=+1)
        )

    def test_oddness_in_mu5_about_diffusion_floor(self):
        # the CME (curl) part is odd in mu5; the diffusion part -eta k^2 is the
        # even baseline. So rate(+mu5) + rate(-mu5) = 2 * (-eta k^2).
        k, mu5, eta = 1.3, 0.7, 1.0
        s = cme_growth_rate(k, +mu5, eta=eta, helicity=+1) + cme_growth_rate(
            k, -mu5, eta=eta, helicity=+1
        )
        assert s == pytest.approx(2.0 * (-eta * k ** 2))

    def test_mu5_zero_is_pure_diffusion(self):
        k, eta = 1.7, 1.3
        assert cme_growth_rate(k, 0.0, eta=eta, helicity=+1) == pytest.approx(
            -eta * k ** 2
        )
        # helicity-independent at mu5 = 0
        assert cme_growth_rate(k, 0.0, eta=eta, helicity=-1) == pytest.approx(
            -eta * k ** 2
        )

    def test_matches_closed_form(self):
        k, mu5, eta, h = 0.9, 1.1, 0.7, +1
        assert cme_growth_rate(k, mu5, eta=eta, helicity=h) == pytest.approx(
            eta * (-(k ** 2) + h * 2.0 * mu5 * k)
        )


# ---------------------------------------------------------------------------
# cme_current -- j = 2 mu5 B
# ---------------------------------------------------------------------------
class TestCMECurrent:
    def test_scalar(self):
        assert cme_current(3.0, 0.5) == pytest.approx(3.0)
        assert cme_current(2.0, -1.0) == pytest.approx(-4.0)

    def test_array(self):
        B = np.array([1.0, -2.0, 0.5])
        np.testing.assert_allclose(cme_current(B, 0.5), B)

    def test_zero_mu5_zero_current(self):
        B = np.array([1.0, 2.0, 3.0])
        np.testing.assert_allclose(cme_current(B, 0.0), 0.0)


# ---------------------------------------------------------------------------
# eta_B_sign -- sgn(eta_B) = sgn(J_parent)
# ---------------------------------------------------------------------------
class TestEtaBSign:
    def test_matter(self):
        assert eta_B_sign(+1) == +1

    def test_antimatter(self):
        assert eta_B_sign(-1) == -1

    def test_default_is_matter(self):
        assert eta_B_sign() == +1

    def test_unspun_parent_no_asymmetry(self):
        assert eta_B_sign(0) == 0


# ---------------------------------------------------------------------------
# net_linking -- the pseudoscalar total linking
# ---------------------------------------------------------------------------
def _circle(radius=1.0, n=240, center=(0, 0, 0), plane="xy"):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    zero = np.zeros_like(t)
    if plane == "xy":
        pts = np.stack([radius * np.cos(t), radius * np.sin(t), zero], axis=1)
    else:  # "xz"
        pts = np.stack([radius * np.cos(t), zero, radius * np.sin(t)], axis=1)
    return pts + np.asarray(center, dtype=float)


class TestNetLinking:
    def test_two_far_apart_rings_zero(self):
        a = _circle(center=(0, 0, 0))
        b = _circle(center=(8, 0, 0))
        assert net_linking([a, b]) == pytest.approx(0.0, abs=1e-3)

    def test_hopf_link_is_unit(self):
        a, b = hopf_link_curves()
        assert abs(net_linking([a, b])) == pytest.approx(1.0, abs=0.02)

    def test_borromean_net_is_zero(self):
        # all three pairwise linkings vanish -> net linking 0 (achiral reference)
        assert net_linking(list(borromean_rings())) == pytest.approx(0.0, abs=1e-2)

    def test_single_curve_is_zero(self):
        assert net_linking([_circle()]) == 0.0

    def test_equals_half_matrix_sum(self):
        from nwt_substrate.topology import linking_matrix

        curves = list(hopf_link_curves())
        assert net_linking(curves) == pytest.approx(
            linking_matrix(curves).sum() / 2.0
        )
