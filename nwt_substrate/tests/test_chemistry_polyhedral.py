"""Tests for substrate polyhedral-cluster identifications (Wade's rules /
3D aromaticity).

Resolves chemistry Tier-B.5 ([[wade-rules-3d-aromaticity-resolution]]):

  Form A WEAK — closo "+1" SEP from generic χ topology (k = 3 − χ).
  Form B PARTIAL — B_5 → B_8 (n, n+1) traverses Spin(7) rep-class ladder
                    EXACTLY (h_v=5, h=6, dim_V=7, dim_S=8, N_POS_ROOTS=9).
  Form C FAILS — substrate silent on closo/nido/arachno beyond χ.
  Form D WINS — partial substrate identification + χ topology fallback.
  Form E REJECTED — p(random 4-consecutive double-hit) ≈ 1.5×10⁻⁴.

  Bonus: B_9 tricapped-trig-prism deltahedron has E = 21 = N_EDGES_K7.
"""
from __future__ import annotations

import pytest

from nwt_substrate.chemistry import (
    CLOSO_POLYHEDRA,
    ClosoCanonicalResult,
    WadeClass,
    closo_borane_sep_count,
    closo_borane_substrate_canonical,
    deltahedron_edge_count,
    deltahedron_euler_chi,
    deltahedron_face_count,
    wade_classification,
)
from nwt_substrate.isa.constants import (
    DIM_OCTONION,
    H_COXETER_SO7,
    H_V_SO7,
    K8_PARTITION,
    N_EDGES_K7,
    N_POS_ROOTS_SO7,
    N_VERTICES_K7,
)


# ---------------------------------------------------------------------------
# Wade-Mingos classification
# ---------------------------------------------------------------------------

class TestWadeClassification:
    def test_closo_b12(self):
        """B_12 H_12^{2-} icosahedron: 12 vertices, 13 SEPs → closo."""
        assert wade_classification(12, 13) == WadeClass.CLOSO

    def test_closo_b6(self):
        """B_6 H_6^{2-} octahedron: 6 vertices, 7 SEPs → closo."""
        assert wade_classification(6, 7) == WadeClass.CLOSO

    def test_nido_b5h9(self):
        """B_5 H_9 pentaborane: 5 vertices, 7 SEPs → nido."""
        assert wade_classification(5, 7) == WadeClass.NIDO

    def test_arachno_b4h10(self):
        """B_4 H_10 tetraborane: 4 vertices, 7 SEPs → arachno."""
        assert wade_classification(4, 7) == WadeClass.ARACHNO

    def test_hypho(self):
        """4 vertices, 8 SEPs → hypho (k=4)."""
        assert wade_classification(4, 8) == WadeClass.HYPHO

    def test_offset_below_closo_raises(self):
        with pytest.raises(ValueError):
            wade_classification(6, 6)   # k=0, below closo minimum

    def test_offset_above_hypho_raises(self):
        with pytest.raises(ValueError):
            wade_classification(5, 10)  # k=5, outside known classes

    def test_too_few_vertices_raises(self):
        with pytest.raises(ValueError):
            wade_classification(3, 4)


# ---------------------------------------------------------------------------
# Closo SEP and deltahedron counts (generic χ topology — Form A weak)
# ---------------------------------------------------------------------------

class TestClosoCounts:
    @pytest.mark.parametrize("n,expected_sep", [
        (5, 6), (6, 7), (7, 8), (8, 9),
        (9, 10), (10, 11), (11, 12), (12, 13),
    ])
    def test_sep_count_matches_canonical_set(self, n, expected_sep):
        assert closo_borane_sep_count(n) == expected_sep

    @pytest.mark.parametrize("n,expected_edges", [
        (5, 9), (6, 12), (7, 15), (8, 18),
        (9, 21), (10, 24), (11, 27), (12, 30),
    ])
    def test_deltahedron_edges(self, n, expected_edges):
        """E = 3n − 6 from Euler V − E + F = 2 and 3F = 2E."""
        assert deltahedron_edge_count(n) == expected_edges

    @pytest.mark.parametrize("n,expected_faces", [
        (5, 6), (6, 8), (7, 10), (8, 12),
        (9, 14), (10, 16), (11, 18), (12, 20),
    ])
    def test_deltahedron_faces(self, n, expected_faces):
        """F = 2n − 4."""
        assert deltahedron_face_count(n) == expected_faces

    @pytest.mark.parametrize("n", list(range(4, 21)))
    def test_euler_characteristic_always_2(self, n):
        """All closo deltahedra have χ = 2 (genus 0)."""
        assert deltahedron_euler_chi(n) == 2
        V = n
        E = deltahedron_edge_count(n)
        F = deltahedron_face_count(n)
        assert V - E + F == 2

    def test_too_few_vertices_raises(self):
        with pytest.raises(ValueError):
            closo_borane_sep_count(3)
        with pytest.raises(ValueError):
            deltahedron_edge_count(3)


# ---------------------------------------------------------------------------
# Substrate-canonical identifications — Form D content
# ---------------------------------------------------------------------------

class TestSpin7RepClassLadder:
    """Lock in the 4-of-4 B_5–B_8 ladder match documented in the
    resolution memo."""

    def test_b5_h_v_h_coxeter(self):
        """B_5 trigonal bipyramid: (5, 6) = (H_V_SO7, H_COXETER_SO7)."""
        r = closo_borane_substrate_canonical(5)
        assert r.n_canonical_label == "H_V_SO7"
        assert r.sep_canonical_label == "H_COXETER_SO7"
        assert r.spin7_ladder_position == "(h_v, h_Coxeter) = (5, 6)"
        assert r.double_canonical is True
        assert r.on_spin7_ladder is True

    def test_b6_h_coxeter_dim_v(self):
        """B_6 octahedron: (6, 7) = (H_COXETER_SO7, N_VERTICES_K7)."""
        r = closo_borane_substrate_canonical(6)
        assert r.n_canonical_label == "H_COXETER_SO7"
        assert r.sep_canonical_label == "N_VERTICES_K7"
        assert r.spin7_ladder_position == "(h_Coxeter, dim_V) = (6, 7)"
        assert r.double_canonical is True

    def test_b7_dim_v_dim_s(self):
        """B_7 pentagonal bipyramid: (7, 8) = (dim_V, dim_S)."""
        r = closo_borane_substrate_canonical(7)
        assert r.n_canonical_label == "N_VERTICES_K7"
        assert r.sep_canonical_label == "DIM_OCTONION"
        assert r.spin7_ladder_position == "(dim_V, dim_S) = (7, 8)"
        assert r.double_canonical is True

    def test_b8_dim_s_n_pos_roots(self):
        """B_8 snub disphenoid: (8, 9) = (dim_S, N_POS_ROOTS)."""
        r = closo_borane_substrate_canonical(8)
        assert r.n_canonical_label == "DIM_OCTONION"
        assert r.sep_canonical_label == "N_POS_ROOTS_SO7"
        assert r.spin7_ladder_position == "(dim_S, N_POS_ROOTS) = (8, 9)"
        assert r.double_canonical is True

    def test_ladder_actual_integer_values(self):
        """Crucial: the ladder integers MUST equal the substrate constants
        (locked at canonical-constant import time)."""
        assert H_V_SO7 == 5
        assert H_COXETER_SO7 == 6
        assert N_VERTICES_K7 == 7
        assert DIM_OCTONION == 8
        assert N_POS_ROOTS_SO7 == 9


class TestClosoB9KSubseven:
    """Bonus finding: B_9 tricapped trig prism has 21 edges = N_EDGES_K7."""

    def test_b9_n_pos_roots(self):
        """B_9: n=9=N_POS_ROOTS_SO7 (vertex side hit; SEP side empty)."""
        r = closo_borane_substrate_canonical(9)
        assert r.n_canonical_label == "N_POS_ROOTS_SO7"
        assert r.sep_canonical_label is None
        assert r.double_canonical is False
        assert r.on_spin7_ladder is False  # ladder runs only B_5–B_8

    def test_b9_edge_count_is_k7_edges(self):
        """B_9 deltahedron has E = 21 = N_EDGES_K7 = dim(so(7)). EXACT."""
        r = closo_borane_substrate_canonical(9)
        assert deltahedron_edge_count(9) == N_EDGES_K7
        assert r.edge_canonical_label == "N_EDGES_K7"

    def test_b9_notes_flag_k7_followup(self):
        """The K_7-edge-count finding is flagged in notes for follow-up."""
        r = closo_borane_substrate_canonical(9)
        assert any("K_7" in note for note in r.notes)


class TestClosoB12Trefoil:
    """B_12 icosahedron: (12, 13) = (K_8 partition[2], trefoil p²+q²)."""

    def test_b12_double_canonical_via_separate_identifications(self):
        r = closo_borane_substrate_canonical(12)
        assert r.n_canonical_label == "K8_PARTITION[2]"
        assert r.sep_canonical_label == "trefoil(p²+q²) = 2²+3² = 13"
        assert r.double_canonical is True
        assert r.on_spin7_ladder is False  # not on Spin(7) ladder

    def test_b12_uses_k8_partition_entry(self):
        """12 = K8_PARTITION[2] — the 12-edge sector of K_8."""
        assert K8_PARTITION[2] == 12

    def test_b12_trefoil_equation(self):
        """13 = 2² + 3² = trefoil quadratic-form integer."""
        assert 2 * 2 + 3 * 3 == 13

    def test_b12_notes_mention_separate_identifications(self):
        r = closo_borane_substrate_canonical(12)
        assert any("K_8" in note and "trefoil" in note for note in r.notes)


class TestClosoEmptyMiddle:
    """B_10 and B_11 are partial-or-zero canonical hits (Form D edge)."""

    def test_b10_zero_canonical(self):
        """B_10 bicapped sq antiprism: (10, 11) — neither integer canonical."""
        r = closo_borane_substrate_canonical(10)
        assert r.n_canonical_label is None
        assert r.sep_canonical_label is None
        assert r.double_canonical is False

    def test_b11_sep_canonical_only(self):
        """B_11 octadecahedron: (11, 12) — sep=12 canonical via K_8."""
        r = closo_borane_substrate_canonical(11)
        assert r.n_canonical_label is None
        assert r.sep_canonical_label == "K8_PARTITION[2]"
        assert r.double_canonical is False


class TestCanonicalPolyhedraTable:
    """The CLOSO_POLYHEDRA registry has the canonical 8-cluster set."""

    def test_full_set_present(self):
        assert set(CLOSO_POLYHEDRA.keys()) == {5, 6, 7, 8, 9, 10, 11, 12}

    def test_b6_is_octahedron(self):
        assert CLOSO_POLYHEDRA[6] == ("octahedron", "O_h")

    def test_b12_is_icosahedron(self):
        assert CLOSO_POLYHEDRA[12] == ("icosahedron", "I_h")


# ---------------------------------------------------------------------------
# Cross-consistency: Form-D verdict
# ---------------------------------------------------------------------------

class TestFormDLadderCount:
    """The Form-D verdict requires exactly 4 of 8 primitive double-hits at
    B_5–B_8 and exactly 5 of 8 with derived (adding B_12 via trefoil 13).
    Lock these in to catch any future regression."""

    def test_primitive_double_hit_count_is_four(self):
        # B_5, B_6, B_7, B_8 on Spin(7) rep-class ladder
        double = sum(
            1 for n in CLOSO_POLYHEDRA
            if (r := closo_borane_substrate_canonical(n)).double_canonical
            and r.spin7_ladder_position is not None
        )
        assert double == 4

    def test_total_double_hit_count_is_five(self):
        # B_5–B_8 (Spin(7) ladder) + B_12 (K_8 + trefoil)
        double = sum(
            1 for n in CLOSO_POLYHEDRA
            if closo_borane_substrate_canonical(n).double_canonical
        )
        assert double == 5

    def test_zero_canonical_count_is_one(self):
        """Only B_10 has zero substrate-canonical identifications among canonical set."""
        zero = sum(
            1 for n in CLOSO_POLYHEDRA
            if (r := closo_borane_substrate_canonical(n)).n_canonical_label is None
            and r.sep_canonical_label is None
        )
        assert zero == 1
