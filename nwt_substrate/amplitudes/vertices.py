"""
Vertex factors for substrate-algebraic Feynman diagrams.

Standard model vertex content reproduced from Cl(1,3) Dirac structure:

  QED (photon-fermion):       -i e gamma^mu        ->  qed_vertex(mu, gammas, e)
  V-A weak (W-fermion):       -i (g/sqrt(2)) gamma^mu (1 - gamma^5) / 2
                              ->  va_vertex(mu, gammas, g)
  Yukawa (Higgs-fermion):     -i (m_f / v) I_8     ->  yukawa_vertex(m_f, v)

These factors slot between spinors u_bar and u in matrix-element
calculations, and into propagator chains for internal lines.
"""

from __future__ import annotations

import numpy as np

from ..algebra.dirac import make_gamma5


# Default electric charge (sign convention: e- has charge -e, so vertex is -i e gamma^mu)
ELECTRIC_CHARGE = 0.30282212  # e = sqrt(4 pi alpha) at low energy, dimensionless
# Weak coupling (rough, approximate)
G_WEAK = 0.6536          # g_W ~ e / sin(theta_W) for sin^2(theta_W) ~ 0.231


def qed_vertex(mu: int, gammas, e: float = ELECTRIC_CHARGE) -> np.ndarray:
    """
    QED vertex factor:  -i e gamma^mu  (8x8 complex matrix).

    Sign:  charged-fermion vertex is -i e gamma^mu where e > 0 is the unit
    charge.  An e- (charge -e) gives this; an e+ flips sign.
    """
    return -1j * e * gammas[mu]


def va_vertex(mu: int, gammas, g: float = G_WEAK) -> np.ndarray:
    """
    V-A weak vertex (charged current):  -i (g / sqrt(2)) gamma^mu P_L
    where P_L = (1 - gamma^5) / 2.

    Used for W^+- couplings to left-handed fermion currents.
    """
    g5 = make_gamma5(gammas)
    I = np.eye(8, dtype=complex)
    P_L = (I - g5) / 2.0
    return -1j * (g / np.sqrt(2.0)) * gammas[mu] @ P_L


def yukawa_vertex(m_f: float, v: float = 246.0) -> np.ndarray:
    """
    Yukawa vertex:  -i (m_f / v) I_8.

    Higgs to fermion coupling.  v is the Higgs VEV in GeV (default 246).
    """
    return -1j * (m_f / v) * np.eye(8, dtype=complex)


def gluon_vertex(mu: int, gammas, gs: float = 1.22) -> np.ndarray:
    """
    Color-stripped gluon-quark vertex:  -i g_s gamma^mu  (without color matrix).

    For an actual amplitude, multiply by the SU(3) generator T^a (lambda^a / 2).
    Color factors are convention-dependent; we leave them to the caller.
    """
    return -1j * gs * gammas[mu]


def three_gluon_vertex_lorentz(k1: np.ndarray, k2: np.ndarray, k3: np.ndarray,
                                eta: np.ndarray | None = None) -> np.ndarray:
    """
    Lorentz tensor structure of the 3-gluon vertex (color factor stripped).

    Standard Feynman rule (Peskin eq. 16.5, all momenta incoming, k1+k2+k3=0):

        V^{mu nu rho}(k1, k2, k3) = (k1 - k2)^rho * eta^{mu nu}
                                   + (k2 - k3)^mu  * eta^{nu rho}
                                   + (k3 - k1)^nu  * eta^{rho mu}

    The full 3-gluon vertex is  V^{mu nu rho}_{abc}  =  g_s * f^{abc} * V_Lorentz.

    Returns a (4, 4, 4) tensor V[mu, nu, rho] in mostly-minus metric.
    Bose-symmetric under (a, mu, k1) <-> (b, nu, k2) <-> (c, rho, k3) when
    combined with f^{abc} (which is totally antisymmetric).
    """
    if eta is None:
        eta = np.diag([1.0, -1.0, -1.0, -1.0])

    V = np.zeros((4, 4, 4), dtype=complex)
    diff_12 = k1 - k2
    diff_23 = k2 - k3
    diff_31 = k3 - k1

    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                V[mu, nu, rho] = (
                    diff_12[rho] * eta[mu, nu]
                    + diff_23[mu]  * eta[nu, rho]
                    + diff_31[nu]  * eta[rho, mu]
                )
    return V


def three_gluon_vertex_full(k1: np.ndarray, k2: np.ndarray, k3: np.ndarray,
                             a: int, b: int, c: int,
                             gs: float, f_abc: np.ndarray | None = None) -> np.ndarray:
    """
    Full 3-gluon vertex including the SU(3) structure-constant color factor.

        V^{mu nu rho}_{abc}  =  g_s * f^{abc} * V_Lorentz(k1, k2, k3)

    Returns a (4, 4, 4) tensor for color indices (a, b, c).
    """
    if f_abc is None:
        from ..algebra.su3 import structure_constants
        f_abc = structure_constants()

    V_lorentz = three_gluon_vertex_lorentz(k1, k2, k3)
    return gs * float(f_abc[a, b, c]) * V_lorentz


def four_gluon_vertex_lorentz(eta: np.ndarray | None = None) -> np.ndarray:
    """
    Color-stripped 4-gluon vertex Lorentz tensors (Peskin eq. 16.6).

    The full vertex has THREE distinct color structures, each with its own
    Lorentz tensor.  We return the three Lorentz tensors as a (3, 4, 4, 4, 4)
    array; the caller multiplies by the appropriate f-product color factor:

        V_{abcd}^{mu nu rho sigma}  =  -i g_s^2 [
            f^{abe} f^{cde} (eta^{mu rho} eta^{nu sigma} - eta^{mu sigma} eta^{nu rho})
          + f^{ace} f^{bde} (eta^{mu nu} eta^{rho sigma} - eta^{mu sigma} eta^{nu rho})
          + f^{ade} f^{bce} (eta^{mu nu} eta^{rho sigma} - eta^{mu rho} eta^{nu sigma})
        ]

    The 4-gluon vertex is needed for gg -> gg at order alpha_s^2.

    Returns the 3 Lorentz tensors L[k][mu, nu, rho, sigma], k in {0, 1, 2}.
    """
    if eta is None:
        eta = np.diag([1.0, -1.0, -1.0, -1.0])

    L = np.zeros((3, 4, 4, 4, 4))
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for sig in range(4):
                    L[0, mu, nu, rho, sig] = (
                        eta[mu, rho] * eta[nu, sig] - eta[mu, sig] * eta[nu, rho]
                    )
                    L[1, mu, nu, rho, sig] = (
                        eta[mu, nu] * eta[rho, sig] - eta[mu, sig] * eta[nu, rho]
                    )
                    L[2, mu, nu, rho, sig] = (
                        eta[mu, nu] * eta[rho, sig] - eta[mu, rho] * eta[nu, sig]
                    )
    return L
