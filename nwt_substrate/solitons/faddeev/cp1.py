"""CP^1 representation of the Faddeev-Skyrme S^2 field.

The Faddeev(-Niemi) model's unit 3-vector field n on S^2 is written in terms of
a 2-component complex spinor Z = (Z1, Z2) as the manifestly gauge-invariant

    n^a = Z* sigma^a Z / (Z* Z) ,

i.e. n = (2 Re(Z1* Z2), 2 Im(Z1* Z2), |Z1|^2 - |Z2|^2) / (|Z1|^2 + |Z2|^2).
This is the L_3 (Skyrme-Faddeev) sector of the NWT Paper 16 master Lagrangian;
its solitons are the knotted/linked Hopf tubes (hopfions) carrying the Hopf
charge Q_H (see `charge.whitehead_hopf_charge`).

The reverse map (a lift n -> Z) is deliberately NOT provided here: the obvious
single-patch lift Z1 = (n1 - i n2)/sqrt(2(1 - n3)) is 0/0 on the hopfion core
ring (n3 -> +1) and is exactly what `seeds.rational_hopfion` is built to avoid.
"""
from __future__ import annotations

import numpy as np


def n_from_Z(Z1, Z2):
    """Unit S^2 field n = Z(dagger) sigma Z / Z(dagger)Z from the CP^1 spinor.

    Z1, Z2 are complex arrays of the same shape; returns (n1, n2, n3) real
    arrays with n1^2 + n2^2 + n3^2 = 1 (Z need not be normalised)."""
    Z1 = np.asarray(Z1)
    Z2 = np.asarray(Z2)
    norm = np.abs(Z1) ** 2 + np.abs(Z2) ** 2
    n1 = 2.0 * np.real(np.conj(Z1) * Z2) / norm
    n2 = 2.0 * np.imag(np.conj(Z1) * Z2) / norm
    n3 = (np.abs(Z1) ** 2 - np.abs(Z2) ** 2) / norm
    return n1, n2, n3
