"""Tests for nwt_substrate.algebra.icosahedral."""

from nwt_substrate.algebra.icosahedral import (
    sl2_f5_elements,
    mul_sl2,
    inverse_sl2,
    order_of_element,
    conjugacy_classes,
    irrep_dimensions_2I,
    distinct_irrep_dimensions_2I,
)


def test_2I_has_120_elements():
    """|SL(2, F_5)| = 120."""
    assert len(sl2_f5_elements()) == 120


def test_identity_is_identity():
    """(1,0,0,1) acts as identity."""
    e = (1, 0, 0, 1)
    for x in sl2_f5_elements()[:10]:
        assert mul_sl2(e, x) == x
        assert mul_sl2(x, e) == x


def test_inverse_is_inverse():
    """g * g^{-1} = e for all g."""
    e = (1, 0, 0, 1)
    for g in sl2_f5_elements():
        ginv = inverse_sl2(g)
        assert mul_sl2(g, ginv) == e
        assert mul_sl2(ginv, g) == e


def test_2I_has_9_conjugacy_classes():
    """2I has 9 conjugacy classes (and hence 9 irreps)."""
    classes = conjugacy_classes()
    assert len(classes) == 9


def test_irrep_dim_squared_sum_is_120():
    """Sum of squares of irrep dimensions equals group order."""
    dims = irrep_dimensions_2I()
    assert sum(d * d for d in dims) == 120


def test_irrep_dimensions_match_2I():
    """2I has irrep dims {1, 2, 2-bar, 3, 3-bar, 4, 4-bar, 5, 6}."""
    dims = irrep_dimensions_2I()
    assert sorted(dims) == [1, 2, 2, 3, 3, 4, 4, 5, 6]


def test_distinct_irrep_dimensions():
    """Distinct irrep dims are {1, 2, 3, 4, 5, 6}."""
    assert distinct_irrep_dimensions_2I() == [1, 2, 3, 4, 5, 6]


def test_element_order_2():
    """The unique element of order 2 in 2I is the central -I = (4, 0, 0, 4)."""
    minus_I = (4, 0, 0, 4)
    assert order_of_element(minus_I) == 2
