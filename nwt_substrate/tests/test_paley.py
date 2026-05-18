"""Tests for Paley (q, k, λ) bifold and AGL(1, q) → residual breaking."""
from __future__ import annotations

import pytest

from nwt_substrate.algebra import (
    AGLBreakingAnalysis,
    PaleyDesign,
    agl_breaking_analysis,
    is_paley_admissible,
    is_prime,
    k7_paley_bifold,
    non_residues,
    paley_design,
    quadratic_residues,
)


# ---------------------------------------------------------------------------
# Primality and admissibility
# ---------------------------------------------------------------------------

class TestPrimality:
    def test_small_primes(self):
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
            assert is_prime(p)

    def test_small_composites(self):
        for n in [0, 1, 4, 6, 8, 9, 10, 12, 14, 15, 21, 25, 27, 49]:
            assert not is_prime(n)

    def test_paley_admissible_q7(self):
        assert is_paley_admissible(7)

    def test_paley_admissible_family(self):
        for q in [7, 11, 19, 23, 31, 43, 47]:
            assert is_paley_admissible(q), f"{q} should be Paley-admissible"

    def test_paley_inadmissible(self):
        # q ≡ 1 (mod 4) primes (Paley-Hadamard family needs q ≡ 3 (mod 4))
        for q in [2, 5, 13, 17, 29, 37, 41]:
            assert not is_paley_admissible(q), f"{q} should NOT be Paley-admissible"


# ---------------------------------------------------------------------------
# Quadratic residues
# ---------------------------------------------------------------------------

class TestQuadraticResidues:
    def test_qr_p7(self):
        assert quadratic_residues(7) == {1, 2, 4}

    def test_nr_p7(self):
        assert non_residues(7) == {3, 5, 6}

    def test_qr_plus_nr_equals_F_star(self):
        for p in [7, 11, 19, 23]:
            qr = quadratic_residues(p)
            nr = non_residues(p)
            assert qr | nr == set(range(1, p))
            assert qr & nr == set()
            assert len(qr) == (p - 1) // 2
            assert len(nr) == (p - 1) // 2

    def test_qr_p11(self):
        # squares mod 11: 1,4,9,5,3 → {1,3,4,5,9}
        assert quadratic_residues(11) == {1, 3, 4, 5, 9}

    def test_qr_p_must_be_prime(self):
        with pytest.raises(ValueError):
            quadratic_residues(9)


# ---------------------------------------------------------------------------
# Paley (7, 3, 1) design
# ---------------------------------------------------------------------------

class TestPaley7Bifold:
    def test_k7_bifold_returns_design(self):
        d = k7_paley_bifold()
        assert isinstance(d, PaleyDesign)
        assert d.q == 7

    def test_k7_parameters(self):
        d = k7_paley_bifold()
        assert d.k == 3
        assert d.lambda_ == 1

    def test_k7_qr_nr(self):
        d = k7_paley_bifold()
        assert d.qr == frozenset({1, 2, 4})
        assert d.nr == frozenset({3, 5, 6})

    def test_k7_blocks_count(self):
        d = k7_paley_bifold()
        # q QR-blocks and q NR-blocks (translates of QR and NR around F_q)
        assert d.num_blocks() == 7
        assert len(d.nr_blocks) == 7

    def test_k7_blocks_size(self):
        d = k7_paley_bifold()
        for block in d.qr_blocks:
            assert len(block) == 3
        for block in d.nr_blocks:
            assert len(block) == 3

    def test_k7_first_block_is_qr(self):
        d = k7_paley_bifold()
        # QR + 0 = QR itself
        assert d.qr_blocks[0] == frozenset({1, 2, 4})

    def test_k7_translates_are_shifts(self):
        d = k7_paley_bifold()
        # QR + 1 = {2, 3, 5}
        assert d.qr_blocks[1] == frozenset({2, 3, 5})
        # QR + 3 = {4, 5, 0} = {0, 4, 5}
        assert d.qr_blocks[3] == frozenset({0, 4, 5})

    def test_k7_2_subsets_in_exactly_lambda_blocks(self):
        # Steiner triple system property: every pair lies in exactly λ = 1 block
        d = k7_paley_bifold()
        # collect the unique blocks (translates of QR — there are 7 of them,
        # and for a Steiner system on 7 points with k=3, b=7)
        unique_blocks = set(d.qr_blocks)
        assert len(unique_blocks) == 7
        # every 2-element subset of {0,..,6} appears in exactly one block
        from itertools import combinations
        pair_counts: dict[frozenset[int], int] = {}
        for block in unique_blocks:
            for pair in combinations(sorted(block), 2):
                key = frozenset(pair)
                pair_counts[key] = pair_counts.get(key, 0) + 1
        for pair, count in pair_counts.items():
            assert count == 1, f"pair {sorted(pair)} appears {count} times, expected 1"
        # there are C(7,2) = 21 pairs
        assert len(pair_counts) == 21


# ---------------------------------------------------------------------------
# AGL(1, q) → residual symmetry-breaking
# ---------------------------------------------------------------------------

class TestAGLBreaking:
    def test_q7_three_generations(self):
        """The headline NWT result: AGL(1, 7) → Z_3 under spin-axis breaking."""
        a = agl_breaking_analysis(7)
        assert a.agl_order == 42
        assert a.vertex_0_choices == 7
        assert a.polarity_choices == 2
        assert a.broken_subgroup_order == 14
        assert a.residual_order == 3
        assert a.residual_label == "Z_3"

    def test_q11_five_generations(self):
        a = agl_breaking_analysis(11)
        assert a.agl_order == 110
        assert a.broken_subgroup_order == 22
        assert a.residual_order == 5

    def test_q19_nine_generations(self):
        a = agl_breaking_analysis(19)
        assert a.agl_order == 342
        assert a.broken_subgroup_order == 38
        assert a.residual_order == 9

    def test_q23_eleven_generations(self):
        a = agl_breaking_analysis(23)
        assert a.agl_order == 23 * 22
        assert a.residual_order == 22 // 2
        assert a.residual_order == 11

    def test_residual_is_half_qminus1(self):
        """General formula: residual = (q - 1) / 2 for any Paley-admissible q."""
        for q in [7, 11, 19, 23, 31, 43]:
            a = agl_breaking_analysis(q)
            assert a.residual_order == (q - 1) // 2

    def test_inadmissible_raises(self):
        with pytest.raises(ValueError):
            agl_breaking_analysis(5)  # 5 mod 4 == 1, not admissible


# ---------------------------------------------------------------------------
# Cross-checks against NWT substrate constants
# ---------------------------------------------------------------------------

class TestNWTSubstrateConsistency:
    def test_k7_paley_matches_substrate_K7(self):
        """The Paley (7, 3, 1) bifold has 21 distinct edges in its
        incidence structure (the 7 lines × 3 points / 1 line-sharing
        = 21 distinct pairs), matching N_EDGES_K7 = 21."""
        from nwt_substrate.isa.constants import N_EDGES_K7
        # every pair of K_7 vertices is on a Paley line
        d = k7_paley_bifold()
        unique_blocks = set(d.qr_blocks)
        from itertools import combinations
        pair_counts: dict[frozenset[int], int] = {}
        for block in unique_blocks:
            for pair in combinations(sorted(block), 2):
                pair_counts[frozenset(pair)] = pair_counts.get(frozenset(pair), 0) + 1
        assert len(pair_counts) == N_EDGES_K7

    def test_q7_residual_matches_three_generations(self):
        """The residual Z_3 has the same order as the three-fermion-generation
        count and the K_7 σ-orbit size used in Paper 20."""
        a = agl_breaking_analysis(7)
        # rank(so(7)) = 3 also; same canonical integer recurs
        from nwt_substrate.isa.constants import RANK_SO7
        assert a.residual_order == RANK_SO7

    def test_agl_full_order_matches_42(self):
        """|AGL(1,7)| = 42 — same as the Z_7 × F_7^* symmetry order
        of the Paley (7, 3, 1) design used in NWT's stabilizer-code work."""
        a = agl_breaking_analysis(7)
        assert a.agl_order == 42
