"""Tests for nwt_substrate.qcd.confinement -- substrate confinement theorem.

P5b closure (Galasyn 2026-05-23): color confinement on the substrate is
a STRUCTURAL CONSEQUENCE of the K_7 graph + SU(3) ⊂ G_2 ⊂ Spin(7) chain
+ closed-walk constraint, not a dynamical mechanism.

The 6 non-polar K_7 vertices split into a color triplet {e_1, e_2, e_3}
and antitriplet {e_5, e_6, e_7} under the v̂ = e_4 SU(3)-singlet
direction.  Every closed walk on K_7 has a well-defined net color
N_c(W); the SU(3)-singlet subset (N_c ≡ 0 mod 3) is the physical
particle spectrum.
"""

import pytest

from nwt_substrate.qcd.confinement import (
    POLAR_VERTEX,
    TRIPLET_VERTICES,
    ANTITRIPLET_VERTICES,
    VERTEX_INDICES,
    color_charge,
    Walk,
    net_color,
    is_singlet_walk,
    enumerate_closed_walks,
    verify_confinement,
    verify_confinement_range,
    GLUEBALL_CANDIDATE_WALK,
    MESON_CANDIDATE_WALK,
    BARYON_CANDIDATE_WALK,
)


# ============================================================
# K_7 vertex partition
# ============================================================

def test_polar_vertex_is_e4():
    """v̂ = e_4 fixes SU(3) inside G_2; corresponds to K_7 polar vertex."""
    assert POLAR_VERTEX == 4


def test_triplet_is_e1_e2_e3():
    """Color triplet (3) = {e_1, e_2, e_3} (+i eigenspace of J = R_{e_4})."""
    assert TRIPLET_VERTICES == (1, 2, 3)


def test_antitriplet_is_e5_e6_e7():
    """Anti-color triplet (3̄) = {e_5, e_6, e_7} (-i eigenspace of J)."""
    assert ANTITRIPLET_VERTICES == (5, 6, 7)


def test_partition_covers_all_seven():
    """The 7 K_7 vertices partition into 1 + 3 + 3."""
    assert set([POLAR_VERTEX]).isdisjoint(TRIPLET_VERTICES)
    assert set([POLAR_VERTEX]).isdisjoint(ANTITRIPLET_VERTICES)
    assert set(TRIPLET_VERTICES).isdisjoint(ANTITRIPLET_VERTICES)
    all_v = {POLAR_VERTEX} | set(TRIPLET_VERTICES) | set(ANTITRIPLET_VERTICES)
    assert all_v == set(VERTEX_INDICES)


# ============================================================
# Color charge of K_7 vertices
# ============================================================

def test_triplet_color_charge_plus_one():
    for v in TRIPLET_VERTICES:
        assert color_charge(v) == +1


def test_antitriplet_color_charge_minus_one():
    for v in ANTITRIPLET_VERTICES:
        assert color_charge(v) == -1


def test_polar_color_charge_zero():
    assert color_charge(POLAR_VERTEX) == 0


def test_color_charge_rejects_invalid_vertex():
    with pytest.raises(ValueError):
        color_charge(0)
    with pytest.raises(ValueError):
        color_charge(8)


# ============================================================
# Walk basics
# ============================================================

def test_walk_no_self_loops():
    """K_7 has no self-loops; walk (1,1) should raise."""
    with pytest.raises(ValueError):
        Walk((1, 1))


def test_walk_rejects_invalid_vertex():
    with pytest.raises(ValueError):
        Walk((1, 8, 1))


def test_walk_length():
    """A walk on n vertices has length n - 1 edges."""
    assert Walk((1, 2, 3, 1)).length == 3
    assert Walk((1, 5)).length == 1


def test_walk_is_closed():
    """A walk is closed iff first == last vertex (and length ≥ 1)."""
    assert Walk((1, 2, 3, 1)).is_closed is True
    assert Walk((1, 2, 3)).is_closed is False


# ============================================================
# Net color
# ============================================================

def test_net_color_of_triplet_3_walk():
    """Walk (1, 2, 3, 1) visits {1, 2, 3} = three triplet vertices.

    Net color = +3 (closing vertex not double-counted).
    """
    assert net_color(Walk((1, 2, 3, 1))) == 3


def test_net_color_of_meson_walk():
    """Walk (1, 5, 1): +1 (triplet) + -1 (antitriplet) = 0."""
    assert net_color(Walk((1, 5, 1))) == 0


def test_net_color_polar_only():
    """Walk through polar vertex contributes 0 to net color."""
    assert net_color(Walk((4, 1, 4))) == 1   # polar(0) + e_1(+1) closure-drop
    assert net_color(Walk((4, 5, 4))) == -1


# ============================================================
# Singlet (substrate physical-state) condition
# ============================================================

def test_glueball_candidate_is_singlet():
    """Triangle on color triplet {1,2,3}: N_c = 3 ≡ 0 (mod 3) ✓."""
    assert is_singlet_walk(GLUEBALL_CANDIDATE_WALK)


def test_meson_candidate_is_singlet():
    """qq̄ analog: +1 + -1 = 0 ✓."""
    assert is_singlet_walk(MESON_CANDIDATE_WALK)


def test_baryon_candidate_is_singlet():
    """3-quark singlet: +1 + +1 + +1 = 3 ≡ 0 ✓."""
    assert is_singlet_walk(BARYON_CANDIDATE_WALK)


def test_non_singlet_walk_detected():
    """Walk (1, 2, 1) has N_c = +2 ≢ 0 (mod 3) — would be color-charged."""
    w = Walk((1, 2, 1))
    assert net_color(w) == 2
    assert not is_singlet_walk(w)


# ============================================================
# Closed-walk enumeration on K_7
# ============================================================

def test_closed_walks_length_2_from_v1():
    """L=2 from v=1: 6 walks (1→x→1, x ∈ {2,3,4,5,6,7})."""
    walks = enumerate_closed_walks(length=2, start=1)
    assert len(walks) == 6
    # All start AND end at 1
    for w in walks:
        assert w.vertices[0] == 1
        assert w.vertices[-1] == 1
        assert w.is_closed


def test_closed_walks_length_3_from_v1():
    """L=3 from v=1: 6·5 = 30 walks (no self-loops, K_7 complete)."""
    walks = enumerate_closed_walks(length=3, start=1)
    assert len(walks) == 30


def test_closed_walks_length_4_from_v1():
    """L=4 from v=1: at each interior step 6 neighbours, last must
    return to 1; counted by transfer-matrix powers."""
    walks = enumerate_closed_walks(length=4, start=1)
    # Total closed walks of length 4 from v=1 on K_7 (complete graph):
    # 6^3 walks of length 3 from v, of which 1/(7-1)=1/6 return... no,
    # walks of length L on K_7 returning to start: (transfer-matrix
    # closed-walk count) — we just check the enumerator's correctness
    # against the simulation.
    expected = 186  # verified by direct count: 6·31 = 186
    assert len(walks) == expected


# ============================================================
# Confinement theorem verification
# ============================================================

def test_confinement_check_length_2():
    """L=2 from v=1: 6 walks total, 3 are singlets (x ∈ {5,6,7})."""
    chk = verify_confinement(length=2, start=1)
    assert chk.total_walks == 6
    assert chk.singlet_walks == 3
    assert chk.fraction_singlet == 0.5


def test_confinement_check_length_3_singlet_count():
    """L=3 from v=1: 30 walks, exactly 8 singlets."""
    chk = verify_confinement(length=3, start=1)
    assert chk.total_walks == 30
    assert chk.singlet_walks == 8


def test_confinement_check_long_walks_approach_one_third():
    """For L → ∞, the singlet fraction approaches 1/3 (uniform
    distribution over Z_3 mod-classes).  This is the substrate-
    monism reading of color confinement: only the 1/3 SU(3)-singlet
    sector corresponds to physical particles."""
    chks = verify_confinement_range(max_length=7, start=1)
    fractions = [chk.fraction_singlet for chk in chks]
    # By L=6 the fraction should be within 2 % of 1/3
    assert abs(fractions[-2] - 1.0 / 3.0) < 0.02
    # By L=7 within 0.5 %
    assert abs(fractions[-1] - 1.0 / 3.0) < 0.005


def test_confinement_range_covers_requested_lengths():
    chks = verify_confinement_range(max_length=4, start=1)
    assert [chk.length for chk in chks] == [2, 3, 4]


def test_every_singlet_walk_has_net_color_divisible_by_three():
    """Cross-check: enumerate ALL singlet walks at L=5 and verify
    each has net color exactly divisible by 3."""
    walks = enumerate_closed_walks(length=5, start=1)
    singlets = [w for w in walks if is_singlet_walk(w)]
    for w in singlets:
        assert net_color(w) % 3 == 0


# ============================================================
# Substrate-monism reading
# ============================================================

def test_no_open_walks_no_free_quarks_principle():
    """The substrate confinement principle: walks must close for the
    state to be a substrate particle.  An open walk like (1, 5) is
    NOT closed and so cannot be a particle — irrespective of its
    color content.  This encodes 'no free quarks' substrate-mechanically.
    """
    open_walk = Walk((1, 5))
    assert open_walk.is_closed is False
    # Color content of an OPEN walk would be ±1 (color or anti-color),
    # which the substrate physical-state constraint forbids by demanding
    # closure.  The walk itself isn't a particle.


def test_polar_only_walks_are_singlets_trivially():
    """Walks confined to the polar singlet (just {e_4}) are length-0
    self-points — not valid walks; closest is (4, x, 4) with x non-polar.
    These are singlet iff the visited non-polar vertex is in a Z_3-
    balanced pattern, which a single visit is not."""
    # (4, 1, 4): visits e_1 once → N_c = +1, NOT singlet.
    assert net_color(Walk((4, 1, 4))) == 1
    assert not is_singlet_walk(Walk((4, 1, 4)))


# ============================================================
# Open-walk semantics
# ============================================================

def test_net_color_on_open_walk_counts_every_vertex():
    """For an OPEN walk, net_color counts every vertex (no closure-drop)."""
    assert net_color(Walk((1, 2, 3))) == 3
    assert net_color(Walk((1, 5))) == 0
    assert net_color(Walk((4, 1))) == 1
