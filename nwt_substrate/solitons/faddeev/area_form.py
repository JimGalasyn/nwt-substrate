"""Faithful (geometric) area form for a lattice S^2 field.

The naive same-index product F_ij = n . (d_i n x d_j n) built from finite
differences is only O(dx)-accurate as a topological density (the two derivatives
live on different half-bonds), carries no lattice barrier, and lets a discrete
Faddeev hopfion unwind to vacuum at every resolution.  The faithful density is
the SIGNED SPHERICAL AREA the field sweeps over each lattice plaquette.  Since
|n| = 1, both d_i n and d_j n are tangent to S^2 so d_i n x d_j n is parallel to
n and F_ij equals the pullback area 2-form; its lattice version is the solid
angle of the spherical quadrilateral spanned by the four corner n-vectors, split
into two triangles (Van Oosterom-Strackee):

    Omega(a, b, c) = 2 atan2( a.(b x c),  1 + a.b + b.c + c.a )
    F_ij(x) dx^2   = Omega(A, B, C) + Omega(A, C, D),
        A = n(x), B = n(x + e_i), C = n(x + e_i + e_j), D = n(x + e_j).

Summed over a closed surface this gives 2*pi * (integer winding): an EXACT
topological density (Berg-Luscher), which is what supplies the lattice
unwinding barrier the naive product lacks.
"""
from __future__ import annotations

import numpy as np


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return np.stack([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])


def solid_angle_triangle(a, b, c):
    """Signed solid angle of the spherical triangle (a, b, c).

    a, b, c are unit 3-vectors stacked along axis 0, shape (3, ...).  Returns
    the per-point solid angle Omega = 2 atan2(a.(b x c), 1 + a.b + b.c + c.a)
    (Van Oosterom-Strackee), in (-2*pi, 2*pi]."""
    num = _dot(a, _cross(b, c))
    den = 1.0 + _dot(a, b) + _dot(b, c) + _dot(c, a)
    return 2.0 * np.arctan2(num, den)


def plaquette_area_form(n, i, j):
    """Faithful area form F_ij * dx^2 on the (i, j) lattice plaquette.

    n is the unit field stacked as (3, N, N, N) (axis 0 is the S^2 component).
    Returns the per-cell solid angle (= F_ij * dx^2) of the n-quadrilateral,
    summed over its two triangles, with periodic (rolled) neighbours."""
    A = n
    B = np.roll(n, -1, axis=1 + i)          # n(x + e_i)
    D = np.roll(n, -1, axis=1 + j)          # n(x + e_j)
    C = np.roll(B, -1, axis=1 + j)          # n(x + e_i + e_j)
    return solid_angle_triangle(A, B, C) + solid_angle_triangle(A, C, D)
