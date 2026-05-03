#!/usr/bin/env python3
"""
Paper 18 -- G1: matter on background metric, linearization, T_mu_nu.

This script sets up the integrand for the Paper 18 Sakharov-induced
gravity calculation.  We:

  (a) write the L1+L2+L3 matter action covariantly with g_mu_nu;
  (b) linearize g_mu_nu = eta_mu_nu + h_mu_nu around flat space;
  (c) identify the homogeneous BPS vacuum and verify T_mu_nu = 0
      there (Minkowski vacuum is stress-energy-free);
  (d) expand T_mu_nu to SECOND order in matter fluctuations around
      the vacuum;
  (e) print the linear graviton-matter coupling -(1/2) h^mu_nu T_mu_nu
      for each sector (psi, A, n);
  (f) verify T^mu_mu (the trace) matches the conformal-anomaly /
      classical-trace expectation.

The output of (d)-(e) is the integrand that G2 will loop over.

==========================================================================
 THE COVARIANT MATTER ACTION
==========================================================================

  S_matter[psi, A, n; g]
    = integral d^4 x sqrt(-g) [
          g^{mu nu} (D_mu psi)^* (D_nu psi)                            [psi kinetic]
        - (1/4) g^{mu alpha} g^{nu beta} F_{mu nu} F_{alpha beta}      [gauge kinetic]
        - (lambda/4) (|psi|^2 - v^2)^2                                 [Higgs potential]
        + (f_pi^2 / 2) g^{mu nu} (d_mu n^a)(d_nu n^a)                  [Skyrme kinetic]
        - (1/(4 e_Sk^2)) (g^{mu alpha} g^{nu beta} - g^{mu beta} g^{nu alpha})
                           (d_mu n^a)(d_nu n^b)(d_alpha n^a)(d_beta n^b)  [Skyrme quartic]
        + theta Q_H[n]                                                  [topological]
      ]

with constraint  n^a n^a = 1,  D_mu psi = (d_mu - i e A_mu) psi,
F_{mu nu} = d_mu A_nu - d_nu A_mu, lambda = e^2/2 (BPS critical).

==========================================================================
 LINEARIZATION g = eta + h
==========================================================================

  g^{mu nu}  =  eta^{mu nu} - h^{mu nu} + O(h^2)
  sqrt(-g)   =  1 + (1/2) h + O(h^2)         where h = eta^{mu nu} h_{mu nu}

So to first order in h:

  S_matter  =  S_matter^{(0)}[psi, A, n]                                  [flat-space action]
            +  integral d^4 x [(1/2) h L_flat - h^{mu nu} T^{(M)}_{mu nu}^{flat}]

where T^{(M)}_{mu nu} is the matter stress-energy in flat space, and the
(1/2) h L_flat part contributes only to the trace-mode coupling.  We
recover the standard form:

  S_int^{(1)}  =  -(1/2) integral d^4 x  h^{mu nu} T^{(M)}_{mu nu}

after using  h^{mu nu} eta_{mu nu} = h  and combining.

==========================================================================
 BPS HOMOGENEOUS VACUUM
==========================================================================

  psi_0   =  v  (real positive)
  A_mu    =  0
  n^a_0   =  delta^{a,3}  (north pole of S^2)

All gradients vanish, |psi_0|^2 = v^2 so the Higgs potential vanishes,
and there is no quartic Skyrme contribution since (d n)^4 = 0.  Therefore

   L_flat^{vacuum} = 0
   T^{(M)}_{mu nu}^{vacuum} = 0

confirming that the matter sector contributes no cosmological constant
classically at the vacuum.  Quantum loops will generate one (Sakharov
predicts Lambda ~ Lambda_UV^4 / (16 pi G)) -- that is the cosmological
constant problem and is tracked separately.

==========================================================================
 SECOND-ORDER FLUCTUATIONS
==========================================================================

Parametrize fluctuations:

  psi  =  v + phi_1 + i phi_2                  (phi_1, phi_2 real scalars)
  A_mu =  A_mu                                 (already a fluctuation)
  n^a  =  (n_1, n_2, sqrt(1 - n_1^2 - n_2^2))  (n_1, n_2 are S^2 Goldstones)
       ~=  (n_1, n_2, 1 - (n_1^2 + n_2^2)/2)   (to second order)

Then to SECOND order in fluctuations, the flat-space matter Lagrangian
splits into kinetic+mass quadratic forms for each fluctuation type
(the contents of S^{(0)}_{matter} for the matter loop in G2).

The script computes the second-order matter stress-energy

  T^{(M),(2)}_{mu nu}  =  (psi sector) + (A sector) + (n sector) + (cross)

which is what couples to h^{mu nu} at the matter-graviton vertex.
"""

from __future__ import annotations

import sympy as sp
from sympy import (
    Symbol, symbols, Function, Rational, sqrt, simplify, expand,
    Matrix, eye, Indexed, IndexedBase, Sum, Idx,
)


# ==========================================================================
# 0. Symbol setup.
# ==========================================================================

# Spacetime indices use 4 components (mu, nu, alpha, beta = 0..3).
DIM = 4

# Physical parameters.
e, v, lam, f_pi, e_Sk, theta = symbols(
    'e v lambda f_pi e_Sk theta', positive=True
)

# At the BPS critical coupling, lambda = e^2 / 2.
# We keep lambda symbolic so the script works at any coupling.

# Coordinates / momenta.
mu, nu, alpha, beta = symbols('mu nu alpha beta', integer=True)

# Field-dependent symbols (one per spacetime point; we'll work with
# their gradients as separate symbols since we're doing the algebra
# at a single point).
phi1, phi2 = symbols('phi_1 phi_2')          # real fluctuations of psi
n1, n2 = symbols('n_1 n_2')                  # S^2 Goldstones
A0, A1, A2, A3 = symbols('A_0 A_1 A_2 A_3')  # gauge field components

# Gradients: we use a covariant representation as 4-component column
# vectors (d_mu f), since at second order in fluctuations everything is
# bilinear in gradients.  This avoids index-tensor pain.
def grad(name: str) -> Matrix:
    """Return a 4x1 column matrix of gradient components d_mu name."""
    return Matrix(symbols(f'{name}_t {name}_x {name}_y {name}_z'))


# Flat metric (mostly-plus signature, eta = diag(-1, 1, 1, 1)).
ETA = sp.diag(-1, 1, 1, 1)


def lower(vec_up: Matrix) -> Matrix:
    """Lower an index using eta (mostly-plus convention)."""
    return ETA * vec_up


def dot(a: Matrix, b: Matrix) -> sp.Expr:
    """eta^{mu nu} a_mu b_nu  =  -a_t b_t + a_x b_x + a_y b_y + a_z b_z."""
    return (lower(a).T * b)[0, 0]


# ==========================================================================
# 1. Flat-space Lagrangian at SECOND ORDER in fluctuations.
# ==========================================================================

def L_psi_flat_quadratic() -> sp.Expr:
    """L_psi at second order in (phi_1, phi_2, A_mu).

    Starting point: |D_mu psi|^2 - (lambda/4)(|psi|^2 - v^2)^2.

    With  psi = v + phi_1 + i phi_2:

      |D_mu psi|^2  =  (d_mu phi_1)(d^mu phi_1)
                    + (d_mu phi_2 - e v A_mu)(d^mu phi_2 - e v A^mu)
                    + O(phi^3, phi^2 A, phi A^2, A^3)

      (|psi|^2 - v^2)^2  =  4 v^2 phi_1^2 + O(phi^3)

    so

      L_psi^{(2)}  =  (d phi_1)^2 + (d phi_2 - e v A)^2 - (lambda v^2) phi_1^2

    The phi_2 - eA combination Higgses A_mu (the Goldstone phi_2 is
    eaten); the radial mode phi_1 has mass m_psi^2 = lambda v^2.
    """
    dphi1 = grad('phi1')
    dphi2 = grad('phi2')
    A = Matrix([A0, A1, A2, A3])

    kin_phi1 = dot(dphi1, dphi1)
    # phi_2 covariant derivative: d_mu phi_2 - e v A_mu
    Dphi2 = dphi2 - e * v * A
    kin_phi2 = dot(Dphi2, Dphi2)

    higgs_mass = -lam * v**2 * phi1**2

    return simplify(kin_phi1 + kin_phi2 + higgs_mass)


def L_A_flat_quadratic() -> sp.Expr:
    """L_A at second order: -(1/4) F_mu_nu F^mu_nu with F = dA - dA.

    For the gauge field at second order (no fluctuation around A=0
    means each A is itself the fluctuation), this is just the Maxwell
    Lagrangian.  We don't expand it explicitly here -- we'll defer the
    F-tensor manipulation to G2.  For G1 we record symbolically.
    """
    # Symbolic placeholder; the full F^2 expansion is standard.
    return Symbol('L_A_quadratic_FF')   # placeholder; computed in G2


def L_n_flat_quadratic() -> sp.Expr:
    """L_n at second order in (n_1, n_2).

    With n^a = (n_1, n_2, 1 - (n_1^2+n_2^2)/2):
      d_mu n^a d^mu n^a  =  d_mu n_1 d^mu n_1 + d_mu n_2 d^mu n_2 + O(n^4)

    So the Skyrme kinetic term contributes (f_pi^2/2)(d n_1)^2 + (f_pi^2/2)(d n_2)^2
    at second order.  The quartic term is O(n^4), zero at this order.
    The Hopf term is topological, no metric dependence.
    """
    dn1 = grad('n1')
    dn2 = grad('n2')
    return simplify((f_pi**2 / 2) * (dot(dn1, dn1) + dot(dn2, dn2)))


# ==========================================================================
# 2. Stress-energy tensor at second order in fluctuations.
# ==========================================================================

def T_psi_quadratic() -> Matrix:
    """T^psi_{mu nu} at second order in (phi_1, phi_2, A).

    For L_psi = (d phi_1)^2 + (d phi_2 - e v A)^2 - lambda v^2 phi_1^2,
    by direct construction or via the metric-variation formula:

      T^psi_{mu nu}  =  2 d_mu phi_1 d_nu phi_1
                      + 2 (d_mu phi_2 - e v A_mu)(d_nu phi_2 - e v A_nu)
                      - eta_{mu nu} L_psi^{(2)}
    """
    dphi1 = grad('phi1')
    dphi2 = grad('phi2')
    A = Matrix([A0, A1, A2, A3])
    Dphi2 = dphi2 - e * v * A

    L = L_psi_flat_quadratic()
    T = sp.zeros(DIM, DIM)
    for m in range(DIM):
        for n in range(DIM):
            # 2 d_mu phi_1 d_nu phi_1 + 2 D_mu phi_2 D_nu phi_2 - eta_mu_nu L
            T[m, n] = (
                2 * dphi1[m] * dphi1[n]
                + 2 * Dphi2[m] * Dphi2[n]
                - ETA[m, n] * L
            )
    return T


def T_n_quadratic() -> Matrix:
    """T^n_{mu nu} at second order in (n_1, n_2)."""
    dn1 = grad('n1')
    dn2 = grad('n2')
    L = L_n_flat_quadratic()
    T = sp.zeros(DIM, DIM)
    for m in range(DIM):
        for n in range(DIM):
            T[m, n] = (
                f_pi**2 * dn1[m] * dn1[n]
                + f_pi**2 * dn2[m] * dn2[n]
                - ETA[m, n] * L
            )
    return T


# ==========================================================================
# 3. Verification: vacuum T_mu_nu = 0
# ==========================================================================

def verify_vacuum_T_zero() -> dict:
    """At the homogeneous vacuum all field gradients vanish and phi=0,
    so T_psi and T_n trivially vanish.  We verify by substitution."""
    T_psi = T_psi_quadratic()
    T_n = T_n_quadratic()

    # Vacuum substitutions: all gradients -> 0, all fluctuations -> 0,
    # all gauge-field components -> 0.
    vacuum_subs = {phi1: 0, phi2: 0, n1: 0, n2: 0,
                    A0: 0, A1: 0, A2: 0, A3: 0}
    for var in ['phi1', 'phi2', 'n1', 'n2']:
        for c in 'txyz':
            vacuum_subs[Symbol(f'{var}_{c}')] = 0

    T_psi_vac = T_psi.subs(vacuum_subs)
    T_n_vac = T_n.subs(vacuum_subs)

    return {
        'T_psi_vacuum_is_zero': T_psi_vac == sp.zeros(DIM, DIM),
        'T_n_vacuum_is_zero': T_n_vac == sp.zeros(DIM, DIM),
    }


# ==========================================================================
# 4. Trace and TT-decomposition at second order.
# ==========================================================================

def trace_T(T: Matrix) -> sp.Expr:
    """eta^{mu nu} T_{mu nu}."""
    return simplify(sum(ETA[m, m] * T[m, m] for m in range(DIM)))


# ==========================================================================
# 5. Graviton-matter coupling.
# ==========================================================================

def graviton_coupling_psi() -> sp.Expr:
    """The linear matter-graviton coupling -(1/2) h^{mu nu} T^psi_{mu nu}.

    We return the symbolic expression; the actual h_mu_nu indices
    (TT, scalar, vector pieces) are extracted in G2 via momentum-space
    analysis of the matter loop.
    """
    h = sp.MatrixSymbol('h', DIM, DIM)
    T = T_psi_quadratic()
    # -1/2 h^{mu nu} T_{mu nu} (h is contravariant by convention here)
    coupling = -Rational(1, 2) * sum(
        h[m, n] * T[m, n] for m in range(DIM) for n in range(DIM)
    )
    return coupling


def graviton_coupling_n() -> sp.Expr:
    h = sp.MatrixSymbol('h', DIM, DIM)
    T = T_n_quadratic()
    return -Rational(1, 2) * sum(
        h[m, n] * T[m, n] for m in range(DIM) for n in range(DIM)
    )


# ==========================================================================
# 6. Output and sanity checks.
# ==========================================================================

def section(title: str) -> None:
    print()
    print('=' * 76)
    print(f' {title}')
    print('=' * 76)


def main() -> None:
    section('Paper 18 -- G1: setup of matter action and T_mu_nu')

    # ---- 1. Quadratic flat-space Lagrangians ----
    section('Quadratic flat-space matter Lagrangians')
    L_psi = L_psi_flat_quadratic()
    L_n = L_n_flat_quadratic()
    print('\n  L_psi^{(2)} (second order in phi_1, phi_2, A):')
    print(f'    = {L_psi}')
    print('\n  L_n^{(2)}   (second order in n_1, n_2):')
    print(f'    = {L_n}')
    print('\n  L_A^{(2)}   = -(1/4) F_{mu nu} F^{mu nu}')
    print('              (Maxwell -- standard, deferred to G2 for explicit')
    print('               momentum-space form)')

    # ---- 2. Stress-energy at second order ----
    section('Stress-energy T_mu_nu at second order')
    T_psi = T_psi_quadratic()
    T_n = T_n_quadratic()

    print('\n  T^psi_{00}  (energy density of psi-sector fluctuations):')
    print(f'    = {simplify(T_psi[0, 0])}')

    print('\n  T^psi_{xx}  (xx pressure of psi-sector fluctuations):')
    print(f'    = {simplify(T_psi[1, 1])}')

    print('\n  T^psi_{0x}  (energy flux):')
    print(f'    = {simplify(T_psi[0, 1])}')

    print('\n  T^n_{00}   (Skyrme energy density):')
    print(f'    = {simplify(T_n[0, 0])}')

    # ---- 3. Trace check ----
    section('Trace eta^mu_nu T_mu_nu  (conformal-anomaly proxy)')
    tr_psi = trace_T(T_psi)
    tr_n = trace_T(T_n)
    print(f'\n  tr T^psi  =  {tr_psi}')
    print(f'\n  tr T^n    =  {tr_n}')
    print('\n  (For a free massless scalar tr T = 0; mass terms break')
    print('   conformal invariance and contribute -m^2 phi^2.  Here the')
    print('   tr T^psi contains the lambda v^2 phi_1^2 mass term.)')

    # ---- 4. Vacuum check ----
    section('Vacuum check: T_mu_nu(vacuum) = 0')
    res = verify_vacuum_T_zero()
    for k, v_ok in res.items():
        status = 'PASS' if v_ok else 'FAIL'
        print(f'  {k}: {status}')

    # ---- 5. Graviton-matter coupling ----
    section('Linear graviton-matter coupling -(1/2) h^mu_nu T_mu_nu')
    cpsi = graviton_coupling_psi()
    cn = graviton_coupling_n()
    print('\n  S_int^{(1)} from psi sector (sketch):')
    print(f'    = {cpsi}')
    print('\n  S_int^{(1)} from n sector (sketch):')
    print(f'    = {cn}')
    print('\n  (The matter loop in G2 will integrate these against the')
    print('   graviton propagator to produce the EH kinetic term.)')

    # ---- 6. Summary ----
    section('Summary of G1')
    print("""
  G1 deliverables:
    1. Matter Lagrangian written covariantly with g_mu_nu.            [done]
    2. Linearization g = eta + h, stress-energy formula identified.   [done]
    3. Homogeneous BPS vacuum: psi_0 = v, A_mu = 0, n^a = delta^{a,3} [done]
    4. T_mu_nu(vacuum) = 0 verified (no classical Lambda).            [verified]
    5. Quadratic-order T_mu_nu computed for psi and n sectors.        [done]
    6. Graviton-matter coupling -(1/2) h^mu_nu T_mu_nu identified.    [done]

  Outputs ready for G2:
    - L_psi^{(2)}  with massive radial mode (m_psi^2 = lambda v^2)
                    and Higgsed phase mode (eaten by A_mu).
    - L_n^{(2)}    with two massless Goldstones (n_1, n_2).
    - L_A^{(2)}    Maxwell, deferred to G2 momentum-space evaluation.
    - T^psi_{mu nu}^{(2)}, T^n_{mu nu}^{(2)} computed symbolically.

  The matter content for the Sakharov loop:
    psi sector:  1 massive scalar (radial Higgs) + 1 massive vector
                 (Higgsed gauge), m^2 = e^2 v^2 for both at BPS.
    n sector:    2 massless Goldstones with derivative interactions.

  Open issues for G2:
    - Skyrme quartic affects loop at higher order; OK to ignore at G2.
    - Explicit Maxwell loop on h-background: standard, transcribe.
    - Heat-kernel a_2 coefficient on flat space: well-known, give
      coefficient of 1/(16 pi G) in Sakharov form.
    - Match to Paper 17's (8/7)^2 alpha^21 (ℏc/m_e^2)^{-1}:  expected
      to FAIL at first order without K_7 Wilson-amplitude.  See
      paper18-einstein-derivation-plan.md risk register.
""")


if __name__ == '__main__':
    main()
