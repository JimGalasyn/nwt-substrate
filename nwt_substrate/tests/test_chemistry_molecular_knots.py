"""Tests for substrate molecular-knot accessibility (Tier B.6).

Resolves chemistry Tier-B.6 ([[molecular-knot-accessibility-resolution]]):

  Form A — strong: all 3 signals (P1, P2, P3) cleanly confirmed — FAILS
  Form B — sweet spot + gap confirmed; 4_1 inversion weakened — WINS
  Form C — qualitative tiers only — DOMINATED by B
  Form D — partial: some signals fail — DOMINATED by B
  Form E — NULL: signal reversed — REJECTED (none reversed)

Substrate carrier-knot table n_q ∈ {0..6} from 2I irrep dimensions
predicts molecular-knot synthesis accessibility:
  - n_q=5 sweet spot (pentafoil 5_1, K(Cl-)≈3.6e10 M^-1)
  - n_q=6 gap (stevedore 6_1, never synthesized as small molecule)
  - n_q=4 below n_q=5 (figure-8 functionally inferior to pentafoil)
"""
from __future__ import annotations

import pytest

from nwt_substrate.chemistry import (
    KNOT_REFERENCE,
    AccessibilityTier,
    HostGuestData,
    MolecularKnotEntry,
    accessibility_for_knot,
    carrier_class_of,
    substrate_predicted_tier,
)
from nwt_substrate.isa.constants import MAX_CROSSING_NUMBER, CARRIER_NAMES
from nwt_substrate.topology.torus_knots import carrier_for_n_q


# ---------------------------------------------------------------------------
# Carrier class mapping
# ---------------------------------------------------------------------------

class TestCarrierClassOf:
    @pytest.mark.parametrize("c,expected", [
        (0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6),
    ])
    def test_in_table(self, c, expected):
        assert carrier_class_of(c) == expected

    @pytest.mark.parametrize("c", [7, 8, 9, 10, 12])
    def test_outside_table_returns_none(self, c):
        assert carrier_class_of(c) is None

    def test_cap_aligns_with_substrate_constant(self):
        """The cap is MAX_CROSSING_NUMBER from substrate constants."""
        assert carrier_class_of(MAX_CROSSING_NUMBER) == MAX_CROSSING_NUMBER
        assert carrier_class_of(MAX_CROSSING_NUMBER + 1) is None

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            carrier_class_of(-1)


# ---------------------------------------------------------------------------
# Substrate-predicted tier — Form B encoding
# ---------------------------------------------------------------------------

class TestSubstratePredictedTier:
    """Lock in the Form-B tier mapping."""

    @pytest.mark.parametrize("c", [0, 1])
    def test_trivial(self, c):
        assert substrate_predicted_tier(c) == AccessibilityTier.TRIVIAL

    @pytest.mark.parametrize("c", [2, 3])
    def test_accessible(self, c):
        assert substrate_predicted_tier(c) == AccessibilityTier.ACCESSIBLE

    def test_4_1_hard(self):
        """n_q=4 figure-8: HARD tier (below sweet spot, P3 weakened)."""
        assert substrate_predicted_tier(4) == AccessibilityTier.HARD

    def test_5_1_sweet_spot(self):
        """n_q=5 pentafoil: SWEET SPOT (P1 confirmed)."""
        assert substrate_predicted_tier(5) == AccessibilityTier.SWEET_SPOT

    def test_6_1_gap(self):
        """n_q=6 stevedore: GAP (P2 confirmed, never synthesized)."""
        assert substrate_predicted_tier(6) == AccessibilityTier.GAP

    @pytest.mark.parametrize("c", [7, 8, 9, 12])
    def test_outside_table(self, c):
        assert substrate_predicted_tier(c) == AccessibilityTier.OUTSIDE_TABLE

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            substrate_predicted_tier(-1)


# ---------------------------------------------------------------------------
# Pentafoil — P1 sweet spot lock-in
# ---------------------------------------------------------------------------

class TestPentafoilSweetSpot:
    """P1 — pentafoil 5_1 is the n_q=5 carrier-knot sweet spot."""

    def test_pentafoil_in_reference(self):
        assert "5_1" in KNOT_REFERENCE

    def test_pentafoil_is_sweet_spot(self):
        r = accessibility_for_knot("5_1")
        assert r.accessibility_tier == AccessibilityTier.SWEET_SPOT

    def test_pentafoil_n_q_is_5(self):
        r = accessibility_for_knot("5_1")
        assert r.n_q_carrier == 5
        assert r.crossing_number == 5

    def test_pentafoil_carrier_name_is_cinquefoil(self):
        r = accessibility_for_knot("5_1")
        assert r.carrier_name == "cinquefoil"
        # Cross-check with the substrate ISA carrier names
        assert CARRIER_NAMES[5] == "cinquefoil"
        assert carrier_for_n_q(5) == "cinquefoil"

    def test_pentafoil_first_synthesis_leigh_2012(self):
        r = accessibility_for_knot("5_1")
        assert r.first_synthesis_year == 2012
        assert "Leigh" in r.first_synthesis_group

    def test_pentafoil_cl_binding_constant(self):
        """K(Cl-) = 3.6e10 M^-1 (Leigh 2015 JACS)."""
        r = accessibility_for_knot("5_1")
        assert r.representative_host_guest is not None
        assert r.representative_host_guest.guest == "Cl-"
        assert r.representative_host_guest.binding_constant_M_inv == pytest.approx(3.6e10)

    def test_pentafoil_serves_p1(self):
        r = accessibility_for_knot("5_1")
        assert "P1_sweet_spot" in r.substrate_signals


# ---------------------------------------------------------------------------
# Stevedore — P2 gap lock-in
# ---------------------------------------------------------------------------

class TestStevedoreGap:
    """P2 — stevedore 6_1 is the n_q=6 gap (never synthesized as small molecule)."""

    def test_stevedore_in_reference(self):
        assert "6_1" in KNOT_REFERENCE

    def test_stevedore_is_gap(self):
        r = accessibility_for_knot("6_1")
        assert r.accessibility_tier == AccessibilityTier.GAP

    def test_stevedore_n_q_is_6(self):
        r = accessibility_for_knot("6_1")
        assert r.n_q_carrier == 6

    def test_stevedore_never_synthesized(self):
        """Lock in: first_synthesis_year is None as of 2026-05-18."""
        r = accessibility_for_knot("6_1")
        assert r.first_synthesis_year is None
        assert r.first_synthesis_group == ""
        assert r.best_yield_pct is None

    def test_stevedore_serves_p2(self):
        r = accessibility_for_knot("6_1")
        assert "P2_gap" in r.substrate_signals


# ---------------------------------------------------------------------------
# Figure-8 — P3 weakened lock-in
# ---------------------------------------------------------------------------

class TestFigureEightWeakened:
    """P3 — 4_1 is HARD (below pentafoil) but no longer 'only via DCL';
    deliberate templates exist (Zhang/Jin 2019/2020/2024)."""

    def test_figure_eight_in_reference(self):
        assert "4_1" in KNOT_REFERENCE

    def test_figure_eight_is_hard(self):
        r = accessibility_for_knot("4_1")
        assert r.accessibility_tier == AccessibilityTier.HARD

    def test_figure_eight_n_q_is_4(self):
        r = accessibility_for_knot("4_1")
        assert r.n_q_carrier == 4
        assert r.carrier_name == "figure-eight"

    def test_figure_eight_has_synthesis_year(self):
        """Sanders 2014 DCL discovery + Zhang/Jin template syntheses since."""
        r = accessibility_for_knot("4_1")
        assert r.first_synthesis_year is not None
        assert "Zhang/Jin" in r.first_synthesis_group or "Sanders" in r.first_synthesis_group

    def test_figure_eight_no_anion_host_data(self):
        """Lock in P3-weakened claim: no anion-binding K reported approaching
        pentafoil's 1e10 — current lit consistent with substrate prediction
        that 4_1 sits below 5_1 in functional terms."""
        r = accessibility_for_knot("4_1")
        # As of 2026-05-18, no documented pentafoil-comparable host data for 4_1
        assert r.representative_host_guest is None

    def test_figure_eight_serves_p3_weakened(self):
        r = accessibility_for_knot("4_1")
        assert "P3_inversion_weakened" in r.substrate_signals

    def test_pentafoil_more_accessible_than_figure_eight(self):
        """P3 qualitative claim preserved: SWEET_SPOT > HARD."""
        tier_5_1 = accessibility_for_knot("5_1").accessibility_tier
        tier_4_1 = accessibility_for_knot("4_1").accessibility_tier
        # SWEET_SPOT corresponds to the substrate-predicted maximum;
        # HARD is below it. The inversion is preserved at the tier level.
        assert tier_5_1 == AccessibilityTier.SWEET_SPOT
        assert tier_4_1 == AccessibilityTier.HARD


# ---------------------------------------------------------------------------
# Outside-table knots
# ---------------------------------------------------------------------------

class TestOutsideTable:
    """8_19, 7_4, 7_1, triskelion: synthesized but outside carrier-knot table."""

    @pytest.mark.parametrize("knot_id", ["7_1", "7_4", "8_19", "3_1#3_1#3_1"])
    def test_outside_table_tier(self, knot_id):
        r = accessibility_for_knot(knot_id)
        assert r.accessibility_tier == AccessibilityTier.OUTSIDE_TABLE
        assert r.n_q_carrier is None

    def test_8_19_synthesized_2017(self):
        r = accessibility_for_knot("8_19")
        assert r.first_synthesis_year == 2017
        assert r.best_yield_pct == 62.0

    def test_8_19_outside_table_consistency(self):
        r = accessibility_for_knot("8_19")
        # Substrate prediction: outside-table OK, since not all synthesis
        # mechanisms go through the carrier-knot vortex closure.
        assert r.crossing_number > MAX_CROSSING_NUMBER


# ---------------------------------------------------------------------------
# Baseline (n_q=2/3) sanity
# ---------------------------------------------------------------------------

class TestBaselineKnots:
    def test_hopf_catenane_n_q_2(self):
        r = accessibility_for_knot("hopf_catenane")
        assert r.crossing_number == 2
        assert r.n_q_carrier == 2
        assert r.carrier_name == "Hopf"
        assert r.accessibility_tier == AccessibilityTier.ACCESSIBLE

    def test_trefoil_n_q_3(self):
        r = accessibility_for_knot("3_1")
        assert r.crossing_number == 3
        assert r.n_q_carrier == 3
        assert r.carrier_name == "trefoil"
        assert r.accessibility_tier == AccessibilityTier.ACCESSIBLE
        assert r.first_synthesis_year == 1989

    def test_5_2_shares_n_q_with_pentafoil(self):
        """5_2 is at n_q=5 too; predicted SWEET_SPOT by tier rule."""
        r = accessibility_for_knot("5_2")
        assert r.n_q_carrier == 5
        assert r.accessibility_tier == AccessibilityTier.SWEET_SPOT


# ---------------------------------------------------------------------------
# Form B verdict lock-in: tier counts in reference table
# ---------------------------------------------------------------------------

class TestFormBVerdict:
    """Lock in that the reference table reflects the Form-B Tier-B.6 verdict."""

    def test_exactly_one_sweet_spot_prime_knot(self):
        """Pentafoil 5_1 is the only prime knot in SWEET_SPOT with
        documented anion-host data."""
        sweet_spots_with_host = [
            r for r in KNOT_REFERENCE.values()
            if r.accessibility_tier == AccessibilityTier.SWEET_SPOT
            and r.representative_host_guest is not None
        ]
        assert len(sweet_spots_with_host) == 1
        assert sweet_spots_with_host[0].alexander_briggs == "5_1"

    def test_exactly_one_gap_entry(self):
        gaps = [r for r in KNOT_REFERENCE.values()
                if r.accessibility_tier == AccessibilityTier.GAP]
        assert len(gaps) == 1
        assert gaps[0].alexander_briggs == "6_1"

    def test_gap_entry_is_unsynthesized(self):
        gap = next(r for r in KNOT_REFERENCE.values()
                   if r.accessibility_tier == AccessibilityTier.GAP)
        assert gap.first_synthesis_year is None

    def test_unknown_knot_raises(self):
        with pytest.raises(KeyError):
            accessibility_for_knot("99_99")


# ---------------------------------------------------------------------------
# Cross-consistency with topology.torus_knots
# ---------------------------------------------------------------------------

class TestCrossConsistency:
    """The carrier-knot names in our entries must match isa.CARRIER_NAMES via
    topology.torus_knots.carrier_for_n_q."""

    @pytest.mark.parametrize("knot_id", ["hopf_catenane", "3_1", "4_1", "5_1", "5_2", "6_1"])
    def test_carrier_name_matches_topology_module(self, knot_id):
        r = accessibility_for_knot(knot_id)
        assert r.n_q_carrier is not None
        assert r.carrier_name == carrier_for_n_q(r.n_q_carrier)

    def test_max_crossing_number_locked(self):
        """The Form-B GAP at n_q=6 IS the substrate's MAX_CROSSING_NUMBER."""
        assert MAX_CROSSING_NUMBER == 6
