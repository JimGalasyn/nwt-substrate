"""
Polarization vectors for gauge bosons and other vector states.

For a massless vector with 4-momentum k:
  - k^2 = 0
  - 2 physical polarizations satisfying k.eps = 0  (transverse)

For a massive vector with 4-momentum k and mass m:
  - k^2 = m^2
  - 3 physical polarizations: 2 transverse + 1 longitudinal,
    all satisfying k.eps = 0.
"""

from __future__ import annotations

import numpy as np


def transverse_polarizations(k_mu: np.ndarray) -> tuple:
    """
    Two real linear transverse polarization 4-vectors for massless k_mu.

    Builds an orthonormal frame (k_hat, x_hat, y_hat) in 3D, then sets
    eps_1^mu = (0, x_hat) and eps_2^mu = (0, y_hat).

    Both satisfy k.eps = 0 (mostly-minus metric):
        k.eps = k^0 eps^0 - k.spatial . eps.spatial
              = 0 - k_hat . (transverse) * |k|
              = 0
    """
    k_spatial = np.real(k_mu[1:])
    k_norm = np.linalg.norm(k_spatial)
    if k_norm < 1e-15:
        raise ValueError("Cannot build transverse polarizations for k=0.")
    k_hat = k_spatial / k_norm

    # Pick a reference direction not parallel to k_hat to take cross-product
    if abs(k_hat[2]) < 0.9:
        ref = np.array([0.0, 0.0, 1.0])
    else:
        ref = np.array([1.0, 0.0, 0.0])
    x_hat = np.cross(k_hat, ref)
    x_hat /= np.linalg.norm(x_hat)
    y_hat = np.cross(k_hat, x_hat)

    eps1 = np.zeros(4, dtype=complex)
    eps2 = np.zeros(4, dtype=complex)
    eps1[1:] = x_hat
    eps2[1:] = y_hat
    return eps1, eps2


def massive_polarizations(k_mu: np.ndarray, m: float) -> tuple:
    """
    Three polarization 4-vectors for a massive vector boson.

    Two transverse (as for massless) plus one longitudinal:
        eps_L^mu = (|k| / m, k_hat * E / m)
    where E = k^0.  Satisfies k.eps_L = (k^0 |k|/m - |k| E/m) = 0  and
    eps_L^2 = -1 in mostly-minus metric.
    """
    k_spatial = np.real(k_mu[1:])
    k_norm = np.linalg.norm(k_spatial)
    E = k_mu[0]
    if k_norm < 1e-15:
        # Rest-frame:  three transverse polarizations along x, y, z
        eps_x = np.array([0.0, 1.0, 0.0, 0.0], dtype=complex)
        eps_y = np.array([0.0, 0.0, 1.0, 0.0], dtype=complex)
        eps_z = np.array([0.0, 0.0, 0.0, 1.0], dtype=complex)
        return eps_x, eps_y, eps_z

    eps1, eps2 = transverse_polarizations(k_mu)
    k_hat = k_spatial / k_norm
    eps_L = np.zeros(4, dtype=complex)
    eps_L[0] = k_norm / m
    eps_L[1:] = (E / m) * k_hat
    return eps1, eps2, eps_L
