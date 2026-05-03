"""
Binary icosahedral group 2I = SL(2, F_5), order 120, and its representation theory.

2I has 9 irreducible representations with dimensions
{1, 2, 2-bar, 3, 3-bar, 4, 4-bar, 5, 6}, summing to
1 + 4 + 4 + 9 + 9 + 16 + 16 + 25 + 36 = 120.

Used by:
  - particles.charge: n_q is dim of a 2I irrep selected by carrier knot
  - compositions.hopf_linked: Hopf links carry 2I rep theory
  - gravity.K7_wilson: K_7 graph state and 2I action
"""

from __future__ import annotations

from itertools import product
from collections import defaultdict


def sl2_f5_elements():
    """All 120 elements of SL(2, F_5) as 4-tuples (a, b, c, d) with ad - bc = 1 mod 5."""
    elts = []
    for a, b, c, d in product(range(5), repeat=4):
        if (a * d - b * c) % 5 == 1:
            elts.append((a, b, c, d))
    return elts


def mul_sl2(x, y):
    """Multiply two SL(2, F_5) elements (mod 5 matrix product)."""
    a, b, c, d = x
    e, f, g, h = y
    return (
        (a * e + b * g) % 5,
        (a * f + b * h) % 5,
        (c * e + d * g) % 5,
        (c * f + d * h) % 5,
    )


def inverse_sl2(x):
    """Inverse in SL(2, F_5): (a,b,c,d) -> (d, -b, -c, a) mod 5."""
    a, b, c, d = x
    return (d % 5, (-b) % 5, (-c) % 5, a % 5)


def order_of_element(x, max_order: int = 200) -> int:
    """Order of an SL(2, F_5) element."""
    identity = (1, 0, 0, 1)
    g = x
    n = 1
    while g != identity:
        g = mul_sl2(g, x)
        n += 1
        if n > max_order:
            raise RuntimeError(f"Order overflow on {x}")
    return n


def conjugacy_classes(elements=None):
    """Partition 2I elements into conjugacy classes."""
    if elements is None:
        elements = sl2_f5_elements()
    classes = []
    seen = set()
    for g in elements:
        if g in seen:
            continue
        cls = set()
        for h in elements:
            cls.add(mul_sl2(mul_sl2(h, g), inverse_sl2(h)))
        classes.append(frozenset(cls))
        seen |= cls
    return classes


def irrep_dimensions_2I() -> list[int]:
    """
    Irreducible representation dimensions of 2I.

    2I has 9 irreps; sum-of-squares = group order 120.
    Distinct dimension values: {1, 2, 3, 4, 5, 6}.
    """
    return [1, 2, 2, 3, 3, 4, 4, 5, 6]


def distinct_irrep_dimensions_2I() -> list[int]:
    """Distinct 2I irrep dimensions {1, 2, 3, 4, 5, 6}."""
    return sorted(set(irrep_dimensions_2I()))


def element_order_distribution() -> dict[int, int]:
    """How many 2I elements have each order."""
    elements = sl2_f5_elements()
    counts = defaultdict(int)
    for g in elements:
        counts[order_of_element(g)] += 1
    return dict(counts)
