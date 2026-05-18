"""Tests for substrate transition-metal electron-count rules (Tier C.7).

Resolves chemistry Tier-C.7 ([[transition-metal-18e-rule-resolution]]):

  Form A WEAK — ladder N_EDGES_K7 − k for k ∈ {3, 5, 7, 9} reproduces
                {18, 16, 14, 12}, but 16/14 have alternative substrate
                forms (uniqueness fails).
  Form B partial — per-geometry table populated.
  Form C FAILS — substrate silent on geometry → count mapping.
  Form D WINS — partial substrate identification, ladder + alternatives.
  Form E REJECTED — {3, 5, 7, 9} ladder is SAME Spin(7) rep-class
                    triple that surfaced in Wade's rules B.5; not coincidence.

  Bonus cross-arc: Fe(CO)_5 (18e, trig bipyramid) and B_5 H_5^{2-}
  (Wade's n=5 closo) share IDENTICAL substrate V/E/F identifications.
"""
from __future__ import annotations

import pytest

from nwt_substrate.chemistry import (
    ElectronCountClass,
    OrganometallicEntry,
    TRANSITION_METAL_REFERENCE,
    closo_borane_substrate_canonical,
    deltahedron_edge_count,
    deltahedron_face_count,
    electron_count_class,
    is_substrate_predicted_stable,
    ladder_k_for_count,
    substrate_canonical_form,
    transition_metal_entry,
)
from nwt_substrate.isa.constants import (
    H_COXETER_SO7,
    H_V_SO7,
    N_EDGES_K7,
    N_POS_ROOTS_SO7,
    N_VERTICES_K7,
    RANK_SO7,
)


# ---------------------------------------------------------------------------
# Spin(7) rep-class ladder N_EDGES_K7 − k = {18, 16, 14, 12}
# ---------------------------------------------------------------------------

class TestSpin7Ladder:
    """Lock in the {N_EDGES_K7 − k} ladder for k ∈ {3, 5, 7, 9}."""

    def test_18_eq_21_minus_3(self):
        assert N_EDGES_K7 - RANK_SO7 == 18
        assert ladder_k_for_count(18) == RANK_SO7

    def test_16_eq_21_minus_5(self):
        assert N_EDGES_K7 - H_V_SO7 == 16
        assert ladder_k_for_count(16) == H_V_SO7

    def test_14_eq_21_minus_7(self):
        assert N_EDGES_K7 - N_VERTICES_K7 == 14
        assert ladder_k_for_count(14) == N_VERTICES_K7

    def test_12_eq_21_minus_9(self):
        assert N_EDGES_K7 - N_POS_ROOTS_SO7 == 12
        assert ladder_k_for_count(12) == N_POS_ROOTS_SO7

    def test_full_ladder_sequence(self):
        """{18, 16, 14, 12} all on the ladder; the {3, 5, 7, 9} integers
        are the first four Spin(7) rep-class invariants."""
        ladder_ks = sorted([RANK_SO7, H_V_SO7, N_VERTICES_K7, N_POS_ROOTS_SO7])
        assert ladder_ks == [3, 5, 7, 9]
        counts = [N_EDGES_K7 - k for k in ladder_ks]
        assert counts == [18, 16, 14, 12]

    def test_non_ladder_count_returns_none(self):
        """20-electron count is not on the substrate ladder."""
        assert ladder_k_for_count(20) is None
        assert ladder_k_for_count(13) is None  # ladder probe at k=8 not enabled

    def test_32_not_on_d_block_ladder(self):
        """32e is f-block (K_7_TRIANGLES − RANK_SO7), not d-block ladder."""
        assert ladder_k_for_count(32) is None


# ---------------------------------------------------------------------------
# Canonical substrate forms
# ---------------------------------------------------------------------------

class TestCanonicalForms:
    @pytest.mark.parametrize("ne,expected_form", [
        (18, "N_EDGES_K7 − RANK_SO7 = 21 − 3"),
        (16, "N_EDGES_K7 − H_V_SO7 = 21 − 5"),
        (14, "N_EDGES_K7 − N_VERTICES_K7 = 21 − 7"),
        (12, "N_EDGES_K7 − N_POS_ROOTS_SO7 = 21 − 9"),
        (32, "K_7_TRIANGLES − RANK_SO7 = 35 − 3"),
    ])
    def test_substrate_canonical_form(self, ne, expected_form):
        assert substrate_canonical_form(ne) == expected_form

    def test_non_predicted_count_returns_none(self):
        assert substrate_canonical_form(20) is None
        assert substrate_canonical_form(17) is None
        assert substrate_canonical_form(50) is None

    def test_is_substrate_predicted_stable(self):
        for ne in [12, 14, 16, 18, 32]:
            assert is_substrate_predicted_stable(ne)
        for ne in [10, 13, 15, 17, 19, 20, 22, 24, 30]:
            assert not is_substrate_predicted_stable(ne)

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            substrate_canonical_form(-1)  # type: ignore[arg-type]  # via electron_count_class
            electron_count_class(-1)


# ---------------------------------------------------------------------------
# Electron-count class enum
# ---------------------------------------------------------------------------

class TestElectronCountClass:
    def test_18e(self):
        assert electron_count_class(18) == ElectronCountClass.EIGHTEEN_E

    def test_16e(self):
        assert electron_count_class(16) == ElectronCountClass.SIXTEEN_E

    def test_14e(self):
        assert electron_count_class(14) == ElectronCountClass.FOURTEEN_E

    def test_12e(self):
        assert electron_count_class(12) == ElectronCountClass.TWELVE_E

    def test_32e(self):
        assert electron_count_class(32) == ElectronCountClass.THIRTYTWO_E

    @pytest.mark.parametrize("ne", [10, 13, 15, 17, 19, 20, 22, 24, 30])
    def test_other(self, ne):
        assert electron_count_class(ne) == ElectronCountClass.OTHER

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            electron_count_class(-1)


# ---------------------------------------------------------------------------
# Canonical 18-electron set
# ---------------------------------------------------------------------------

class TestEighteenElectronSet:
    CANONICAL_18E = [
        "Cr(CO)6", "Fe(CO)5", "Ni(CO)4", "Mn2(CO)10", "Fe(Cp)2",
        "HMn(CO)5", "Co(NH3)6 3+", "V(CO)6 -",
    ]

    def test_all_have_18_electrons(self):
        for formula in self.CANONICAL_18E:
            entry = transition_metal_entry(formula)
            assert entry.electron_count == 18

    def test_all_are_18e_class(self):
        for formula in self.CANONICAL_18E:
            entry = transition_metal_entry(formula)
            assert entry.rule_class == ElectronCountClass.EIGHTEEN_E

    def test_all_have_substrate_form(self):
        for formula in self.CANONICAL_18E:
            entry = transition_metal_entry(formula)
            assert entry.substrate_form == "N_EDGES_K7 − RANK_SO7 = 21 − 3"
            assert entry.on_substrate_ladder


# ---------------------------------------------------------------------------
# Canonical 16-electron set
# ---------------------------------------------------------------------------

class TestSixteenElectronSet:
    CANONICAL_16E = [
        "Pt(NH3)2Cl2", "RhCl(PPh3)3", "IrCl(CO)(PPh3)2",
        "Ni(CN)4 2-", "Pd(PPh3)4",
    ]

    def test_all_have_16_electrons(self):
        for formula in self.CANONICAL_16E:
            entry = transition_metal_entry(formula)
            assert entry.electron_count == 16

    def test_all_are_16e_class(self):
        for formula in self.CANONICAL_16E:
            entry = transition_metal_entry(formula)
            assert entry.rule_class == ElectronCountClass.SIXTEEN_E
            assert entry.substrate_form == "N_EDGES_K7 − H_V_SO7 = 21 − 5"

    def test_wilkinson_square_planar(self):
        """Wilkinson's catalyst RhCl(PPh3)3."""
        e = transition_metal_entry("RhCl(PPh3)3")
        assert e.geometry == "square planar"
        assert e.d_or_f_config == "d8"

    def test_vaska_square_planar(self):
        """Vaska's complex IrCl(CO)(PPh3)2."""
        e = transition_metal_entry("IrCl(CO)(PPh3)2")
        assert e.geometry == "square planar"


# ---------------------------------------------------------------------------
# Canonical 14-electron set
# ---------------------------------------------------------------------------

class TestFourteenElectronSet:
    CANONICAL_14E = ["W(CO)3(PCy3)2", "Pt(PCy3)2", "TiCp2Cl2"]

    def test_all_have_14_electrons(self):
        for formula in self.CANONICAL_14E:
            assert transition_metal_entry(formula).electron_count == 14

    def test_all_are_14e_class(self):
        for formula in self.CANONICAL_14E:
            e = transition_metal_entry(formula)
            assert e.rule_class == ElectronCountClass.FOURTEEN_E
            assert e.substrate_form == "N_EDGES_K7 − N_VERTICES_K7 = 21 − 7"


# ---------------------------------------------------------------------------
# 12-electron probe set (independent test from pre-reg)
# ---------------------------------------------------------------------------

class TestTwelveElectronProbe:
    """Pre-reg locked in k=N_POS_ROOTS_SO7=9 as the ladder extension to 12e.
    Validate against documented bent lanthanide metallocenes."""

    def test_cp2sm_is_12e(self):
        e = transition_metal_entry("Cp2Sm")
        assert e.electron_count == 12
        assert e.rule_class == ElectronCountClass.TWELVE_E

    def test_cp2eu_is_12e(self):
        e = transition_metal_entry("Cp2Eu")
        assert e.electron_count == 12

    def test_12e_substrate_form_uses_n_pos_roots(self):
        e = transition_metal_entry("Cp2Sm")
        assert e.substrate_form == "N_EDGES_K7 − N_POS_ROOTS_SO7 = 21 − 9"


# ---------------------------------------------------------------------------
# 32-electron f-block set
# ---------------------------------------------------------------------------

class TestThirtyTwoElectronSet:
    def test_cerocene(self):
        e = transition_metal_entry("Ce(COT)2")
        assert e.electron_count == 32
        assert e.rule_class == ElectronCountClass.THIRTYTWO_E
        assert e.substrate_form == "K_7_TRIANGLES − RANK_SO7 = 35 − 3"

    def test_uranocene(self):
        e = transition_metal_entry("U(COT)2")
        assert e.electron_count == 32
        assert e.geometry == "sandwich"

    def test_32e_reuses_periodic_table_form(self):
        """32 = K_7_TRIANGLES − RANK_SO7 is the SAME identification as the
        periodic-table A.3 shell-32 result."""
        e = transition_metal_entry("Ce(COT)2")
        assert "K_7_TRIANGLES" in e.substrate_form
        assert "RANK_SO7" in e.substrate_form


# ---------------------------------------------------------------------------
# Cross-arc: Fe(CO)_5 ≡ B_5 H_5^{2-} trigonal bipyramid
# ---------------------------------------------------------------------------

class TestFeCOFiveB5CrossArc:
    """Fe(CO)_5 (18e organometallic, trigonal bipyramid) and B_5 H_5^{2-}
    (Wade's n=5 closo, trigonal bipyramid deltahedron) share IDENTICAL
    substrate V/E/F identifications: V=5=H_V_SO7, E=9=N_POS_ROOTS_SO7,
    F=6=H_COXETER_SO7."""

    def test_fe_co_5_is_trig_bipyramid_18e(self):
        e = transition_metal_entry("Fe(CO)5")
        assert e.geometry == "trig bipyramid"
        assert e.electron_count == 18

    def test_trig_bipyramid_v_e_f_substrate_canonical(self):
        """V=5, E=9, F=6 — all substrate-canonical."""
        V, E, F = 5, 9, 6
        assert V == H_V_SO7
        assert E == N_POS_ROOTS_SO7
        assert F == H_COXETER_SO7

    def test_b_5_deltahedron_matches_trig_bipyramid(self):
        """B_5 H_5^{2-} closo deltahedron is the trigonal bipyramid:
        V=5, E=9, F=6, identical to Fe(CO)_5 coordination polyhedron."""
        n = 5  # B_5 closo
        r = closo_borane_substrate_canonical(n)
        assert r.n_canonical_label == "H_V_SO7"
        # Edge count
        assert deltahedron_edge_count(n) == N_POS_ROOTS_SO7   # = 9
        # Face count (from Euler: F = 2n - 4)
        assert deltahedron_face_count(n) == H_COXETER_SO7      # = 6

    def test_cross_arc_substrate_identity(self):
        """The substrate algebra picks out the same V/E/F integers for
        Fe(CO)_5 (18e organometallic) and B_5 H_5^{2-} (Wade's n=5 closo).
        This is the chemistry-sector cross-arc unification."""
        fe = transition_metal_entry("Fe(CO)5")
        b5 = closo_borane_substrate_canonical(5)
        # Both share the trigonal-bipyramid geometry
        assert fe.geometry == "trig bipyramid"
        # B_5 sits at the Spin(7) rep-class ladder (h_v, h_Coxeter) position
        assert b5.spin7_ladder_position == "(h_v, h_Coxeter) = (5, 6)"


# ---------------------------------------------------------------------------
# Form D verdict — reference table integrity
# ---------------------------------------------------------------------------

class TestFormDVerdict:
    """Lock in the Form-D verdict counts in the reference table."""

    def test_reference_set_size(self):
        """Reference contains all canonical classes: 18e (8) + 16e (5) +
        14e (3) + 12e (2) + 32e (2) = 20 entries."""
        assert len(TRANSITION_METAL_REFERENCE) == 20

    def test_18e_count(self):
        n = sum(1 for e in TRANSITION_METAL_REFERENCE.values()
                if e.rule_class == ElectronCountClass.EIGHTEEN_E)
        assert n == 8

    def test_16e_count(self):
        n = sum(1 for e in TRANSITION_METAL_REFERENCE.values()
                if e.rule_class == ElectronCountClass.SIXTEEN_E)
        assert n == 5

    def test_14e_count(self):
        n = sum(1 for e in TRANSITION_METAL_REFERENCE.values()
                if e.rule_class == ElectronCountClass.FOURTEEN_E)
        assert n == 3

    def test_12e_count(self):
        n = sum(1 for e in TRANSITION_METAL_REFERENCE.values()
                if e.rule_class == ElectronCountClass.TWELVE_E)
        assert n == 2

    def test_32e_count(self):
        n = sum(1 for e in TRANSITION_METAL_REFERENCE.values()
                if e.rule_class == ElectronCountClass.THIRTYTWO_E)
        assert n == 2

    def test_all_entries_on_substrate_ladder(self):
        """All entries in the reference table have substrate-canonical form."""
        for e in TRANSITION_METAL_REFERENCE.values():
            assert e.substrate_form is not None
            assert e.on_substrate_ladder

    def test_unknown_complex_raises(self):
        with pytest.raises(KeyError):
            transition_metal_entry("XYZ(unknown)")
