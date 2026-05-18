"""Tests for substrate NMR ring-current predictions (Tier C.8).

Resolves chemistry Tier-C.8 ([[nmr-via-hopf-pair-resolution]]):

  Form A EXACT — sign rule via Hopf-pair parity (14/14), but TRIVIAL
                 extension of A.1/A.2/B.4 aromaticity classification
  Form B FAILS — broad magnitude prediction fails rational-density audit
                 (7/14 strict hits at ±0.5 ppm gives p_random ≈ 0.78)
  Form C DEFERRED — Tr(M^k) invariant extension is future work
  Form D WINS (narrow) — sign rule + 2 structurally distinctive hits
                         (coronene K_7 outer, benzene DIM_OCTONION)
  Form E NOT REJECTED for broad pattern, but narrow distinctive hits
                       survive via independent substrate identifications

Honest verdict: this is the WEAKEST result in the chemistry-sector survey.
Library ships minimal sign rule + 2 narrow structurally distinctive hits;
NO broad magnitude API (would over-claim).
"""
from __future__ import annotations

import pytest

from nwt_substrate.chemistry import (
    NICS_REFERENCE,
    NICSReference,
    NICSSign,
    StructurallyDistinctiveHit,
    nics_reference,
    nics_sign_from_hopf_parity,
)
from nwt_substrate.isa.constants import DIM_OCTONION, N_EDGES_K7, RANK_SO7


# ---------------------------------------------------------------------------
# Sign rule (Form A — exact but trivial)
# ---------------------------------------------------------------------------

class TestSignRule:
    """The Hopf-pair parity rule predicts NICS sign exactly.  This is
    TRIVIAL by aromaticity classification — locking in for completeness."""

    @pytest.mark.parametrize("n_pi", [2, 6, 10, 14, 18])
    def test_4n_plus_2_diatropic(self, n_pi):
        assert nics_sign_from_hopf_parity(n_pi) == NICSSign.DIATROPIC

    @pytest.mark.parametrize("n_pi", [4, 8, 12, 16, 20])
    def test_4n_paratropic(self, n_pi):
        assert nics_sign_from_hopf_parity(n_pi) == NICSSign.PARATROPIC

    def test_zero_pi_is_zero(self):
        assert nics_sign_from_hopf_parity(0) == NICSSign.ZERO

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            nics_sign_from_hopf_parity(-2)

    def test_odd_raises(self):
        """Cyclic closed-shell π systems must have even electron count."""
        with pytest.raises(ValueError):
            nics_sign_from_hopf_parity(5)


# ---------------------------------------------------------------------------
# Reference set sign rule — 14/14 exact across the canonical NICS set
# ---------------------------------------------------------------------------

class TestReferenceSetSignRule:
    """Form A — sign rule holds across the canonical 14-molecule NICS set."""

    def test_reference_set_size(self):
        assert len(NICS_REFERENCE) == 14

    def test_all_signs_correct(self):
        """The substrate Hopf-pair-parity sign rule must match the actual
        NICS sign for ALL 14 reference molecules."""
        for name, entry in NICS_REFERENCE.items():
            assert entry.sign_rule_correct, (
                f"Sign rule fails for {name}: predicted={entry.predicted_sign}, "
                f"actual={entry.actual_sign} (NICS={entry.nics_ppm})"
            )

    def test_aromatic_count(self):
        """12 aromatic (diatropic, negative NICS) molecules in the reference."""
        n = sum(1 for e in NICS_REFERENCE.values()
                if e.predicted_sign == NICSSign.DIATROPIC)
        assert n == 12

    def test_anti_aromatic_count(self):
        """2 anti-aromatic (paratropic, positive NICS) molecules."""
        n = sum(1 for e in NICS_REFERENCE.values()
                if e.predicted_sign == NICSSign.PARATROPIC)
        assert n == 2


# ---------------------------------------------------------------------------
# Structurally distinctive Form-D hits
# ---------------------------------------------------------------------------

class TestCoroneneK7Hit:
    """The coronene-outer NICS ≈ -18 identification is structurally
    distinctive: same integer 18 = N_EDGES_K7 − RANK_SO7 appears in
    periodic-table A.3, C.7 18e rule, and so7-ISA probe."""

    def test_coronene_outer_in_reference(self):
        assert "coronene_outer" in NICS_REFERENCE

    def test_coronene_outer_distinctive_flag(self):
        e = nics_reference("coronene_outer")
        assert e.structurally_distinctive_hit == StructurallyDistinctiveHit.CORONENE_K7_OUTER

    def test_coronene_outer_nics_close_to_target(self):
        """NICS = -18.7 ppm, within 1 ppm of -(N_EDGES_K7 - RANK_SO7) = -18."""
        e = nics_reference("coronene_outer")
        target = -(N_EDGES_K7 - RANK_SO7)   # = -18
        assert target == -18
        assert abs(e.nics_ppm - target) < 1.0   # 0.7 ppm

    def test_coronene_outer_notes_mention_k7(self):
        e = nics_reference("coronene_outer")
        assert "K_7" in e.notes or "K_7-hub" in e.notes


class TestBenzeneDimOctonionHit:
    """Benzene NICS = -8 = -DIM_OCTONION exactly. Structurally distinctive
    because DIM_OCTONION grounds particle physics throughout NWT."""

    def test_benzene_in_reference(self):
        assert "benzene" in NICS_REFERENCE

    def test_benzene_distinctive_flag(self):
        e = nics_reference("benzene")
        assert e.structurally_distinctive_hit == StructurallyDistinctiveHit.BENZENE_DIM_OCTONION

    def test_benzene_nics_exact_minus_dim_octonion(self):
        e = nics_reference("benzene")
        assert e.nics_ppm == -float(DIM_OCTONION)   # exact -8.0


class TestStructurallyDistinctiveCount:
    """Exactly 2 structurally distinctive hits in the reference table —
    no others over-claimed."""

    def test_exactly_two_distinctive_hits(self):
        n = sum(1 for e in NICS_REFERENCE.values()
                if e.structurally_distinctive_hit is not None)
        assert n == 2

    def test_other_molecules_have_no_distinctive_flag(self):
        """Naphthalene, anthracene, Cp-, etc. are NOT flagged as distinctive
        even though they happen to land near canonical integers — broad
        magnitude prediction failed rational-density audit."""
        for name in ["naphthalene_center", "anthracene_center", "Cp- anion",
                     "phenanthrene_center", "pyrene_center", "coronene_hub",
                     "pyrrole", "furan", "thiophene"]:
            e = nics_reference(name)
            assert e.structurally_distinctive_hit is None


# ---------------------------------------------------------------------------
# Form B/E — broad magnitude API is NOT exposed
# ---------------------------------------------------------------------------

class TestNoMagnitudeAPI:
    """The audit found that broad magnitude prediction fails rational-
    density (p_random ≈ 0.78 for 7/14 hits). The library deliberately does
    NOT expose a broad magnitude API — that would over-claim. This test
    locks in the honest-scope decision."""

    def test_no_predict_magnitude_function(self):
        """No public function should predict NICS magnitudes for arbitrary
        molecules from substrate algebra alone."""
        from nwt_substrate.chemistry import nmr
        # The module should NOT export a function like
        # `nics_magnitude_predict` or `nics_substrate_magnitude`.
        forbidden = [
            "nics_magnitude_predict",
            "nics_substrate_magnitude",
            "predict_nics_value",
            "substrate_nics_prediction",
        ]
        for name in forbidden:
            assert not hasattr(nmr, name), (
                f"NMR module exposes {name!r}, which would over-claim the "
                f"broad magnitude prediction that FAILS rational-density audit "
                f"(p_random ≈ 0.78). Form D requires keeping the scope narrow."
            )


# ---------------------------------------------------------------------------
# Cross-arc with other modules
# ---------------------------------------------------------------------------

class TestCrossArcWithC7AndA3:
    """The integer 18 appears in three independent substrate contexts:
    A.3 periodic-table shell-18, C.7 18-electron rule, and C.8 coronene
    K_7 outer NICS. Cross-consistency lock-in."""

    def test_18_eq_n_edges_k7_minus_rank_so7(self):
        from nwt_substrate.isa.constants import N_EDGES_K7, RANK_SO7
        assert N_EDGES_K7 - RANK_SO7 == 18

    def test_coronene_outer_uses_same_18(self):
        from nwt_substrate.chemistry import nics_reference
        from nwt_substrate.isa.constants import N_EDGES_K7, RANK_SO7
        e = nics_reference("coronene_outer")
        target = -(N_EDGES_K7 - RANK_SO7)
        assert abs(e.nics_ppm - target) < 1.0

    def test_c7_substrate_form_uses_same_18(self):
        from nwt_substrate.chemistry import substrate_canonical_form
        # The transition-metal 18-electron rule uses the same
        # N_EDGES_K7 − RANK_SO7 identification.
        form_18 = substrate_canonical_form(18)
        assert form_18 is not None
        assert "N_EDGES_K7" in form_18 and "RANK_SO7" in form_18


# ---------------------------------------------------------------------------
# Unknown molecule error
# ---------------------------------------------------------------------------

class TestUnknownMolecule:
    def test_unknown_raises(self):
        with pytest.raises(KeyError):
            nics_reference("xyz_unknown_molecule")
