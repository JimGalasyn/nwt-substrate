#!/usr/bin/env python3
"""
Paper 18 -- G2: matter loop contribution to 1/(16 pi G) via heat-kernel a_2.

G1 set up the matter integrand:  matter content = 1 massive real scalar
(radial Higgs, m^2 = lambda v^2), 1 massive Higgsed vector (3
polarisations, m^2 = e^2 v^2), 2 massless Goldstones.

G2 computes the long-wavelength one-loop matter contribution to the
graviton kinetic term, equivalent to the Seeley-DeWitt a_2 coefficient
on a flat background with an explicit UV cutoff Lambda_UV.

==========================================================================
 STANDARD SAKHAROV RESULT
==========================================================================

For a single matter field with kinetic operator  -box + m^2 + xi R  on
a curved background, the heat-kernel expansion gives the one-loop
effective action

  Gamma_1-loop[g]  =  -(1/2) integral_0^inf  (dt/t)  Tr[e^{-t K}]
                   =  -(1/2) integral d^d x sqrt(g) [a_0 + a_1 t + a_2 t^2 + ...]
                       integrated against (4 pi t)^{-d/2}

In d=4, after regulating with cutoff Lambda_UV (proper-time t_min ~
1/Lambda_UV^2):

  Gamma_1-loop  =  - (1/(32 pi^2)) integral d^4 x sqrt(g) [
                       (Lambda_UV^4/2) a_0
                    +  (Lambda_UV^2)   a_1
                    +  ln(Lambda_UV^2/m^2) a_2  +  finite ]

For minimally coupled massive scalar (xi = 0):
   a_0 = 1
   a_1 = R/6 - m^2
   a_2 = (1/180) (R^2 - R_mu_nu^2 + R_mu_nu_rho_sig^2) + (1/2)(m^2 - R/6)^2
        + (1/30) box R + ...

The R-coefficient in the long-wavelength expansion is the EH coupling:

  Gamma  ~  - integral d^4 x sqrt(g) (1/(16 pi G)) R
   ==>  1/(16 pi G)  =  - (1/(32 pi^2)) [- (Lambda_UV^2)/6 + (m^2/6) ln(Lambda_UV^2/m^2) + ...]
                      =  Lambda_UV^2/(192 pi^2)  -  (m^2)/(192 pi^2) ln(...)  + ...

So per scalar species, with cutoff Lambda_UV:

   Delta(1/(16 pi G))_{scalar}  =  Lambda_UV^2 / (192 pi^2)
                                  - m^2 ln(Lambda_UV^2/m^2) / (192 pi^2)

==========================================================================
 EFFECTIVE NUMBER OF DOF FOR NWT
==========================================================================

In 4D flat-space Sakharov, the leading Lambda^2 contribution scales
with the effective number of bosonic DOF (each contributes
Lambda^2 / 192 pi^2, with possible different coefficients for spin):

  Massive real scalar          : 1 DOF -> coefficient 1
  Massive Higgsed vector       : 3 DOF -> coefficient 3 * (vector norm.)
  Massless real (Goldstone)    : 1 DOF -> coefficient 1
  FP ghost (negative metric)   : -2 per gauge field

Standard formula (Vassilevich review / Birrell-Davies):
   coefficient for spin-J in d=4: c_J = (-1)^{2J} (2J+1) - "ghost subtraction"
   spin-0 : +1
   spin-1 (massive) : +3
   spin-2 graviton on flat space: +5  (only relevant if we self-include grav)

For NWT matter content (L1+L2+L3 around BPS vacuum):
   1 massive scalar       contributes  +1
   1 massive Higgsed vector  +3
   2 massless Goldstones  +2
   total                  +6

(The Skyrme quartic gives derivative interactions that contribute at
higher loop order, ignored at this stage.  A_mu's Faddeev-Popov ghosts
have been absorbed into the Higgsed-vector counting.)

==========================================================================
 TARGET: PAPER 17's STRUCTURAL G
==========================================================================

Paper 17's structural derivation:

   m_e/m_Pl  =  (8/7) (1 + alpha/7 + 3 alpha^2) alpha^{21/2}

(with the alpha^2 coefficient verified on Heron to 3.07% via the
K_9/K_7 ZNE ratio test).  In natural units hbar = c = 1:

   G m_e^2  =  (m_e/m_Pl)^2  =  (8/7)^2 (1 + alpha/7)^2 alpha^{21}
   1/(16 pi G)  =  m_Pl^2 / (16 pi)
                =  m_e^2 / [16 pi (8/7)^2 (1 + alpha/7)^2 alpha^{21}]

Setting Sakharov's expression for 1/(16 pi G) equal to this and
solving for Lambda_UV:

   N_DOF Lambda_UV^2 / 192 pi^2  =  m_e^2 / [16 pi (8/7)^2 (1+alpha/7)^2 alpha^{21}]
   Lambda_UV^2  =  (192 pi^2 / 16 pi N_DOF) * m_e^2 / [(8/7)^2 (1+alpha/7)^2 alpha^{21}]
   Lambda_UV    =  m_e * sqrt[(12 pi / N_DOF) / ((8/7)^2 (1+alpha/7)^2 alpha^{21})]

which we evaluate numerically.

==========================================================================
 INTERPRETATION
==========================================================================

If Lambda_UV comes out at a known NWT scale (e.g., E_GUT = 7.4e15 GeV
from L1b's UV completion, or the Compton-scale of the heaviest
fundamental excitation), Sakharov's naive coefficient is consistent
with Paper 17 up to a multiplicative factor that is ABSORBED into
the K_7 amplitude structure.

If Lambda_UV is unrelated to any known NWT scale, then the K_7
amplitude is doing nontrivial work that Sakharov does not capture --
exactly the gap Phase 5 flagged.
"""

from __future__ import annotations

import math


# ==========================================================================
# 1. NWT physical inputs.
# ==========================================================================

ALPHA = 7.2973525693e-3                # CODATA fine-structure constant
M_E_GEV = 0.51099895069e-3            # m_e in GeV
E_GUT_GEV = 7.4e15                    # NWT GUT scale (Paper 14)
M_PL_GEV = 1.220890e19                # reduced Planck mass; we use full M_Pl
HBAR_C_GEV_FM = 0.1973                 # ℏc in GeV*fm
G_NEWTON = 6.67430e-11                 # m^3 kg^-1 s^-2  (SI, for cross-checks)


def paper17_G_natural() -> dict:
    """Paper 17 structural prediction G m_e^2 = (m_e/m_Pl)^2.

    Returns dimensionless G in natural units where hbar = c = m_e = 1.
    """
    nlo = 1.0 + ALPHA / 7.0  # NLO bracket; ignore alpha^2 NNLO at LO
    prefac = (8.0 / 7.0) ** 2 * nlo ** 2
    alpha21 = ALPHA ** 21
    return dict(
        prefactor=prefac,
        alpha21=alpha21,
        G_natural=prefac * alpha21,
        m_e_over_m_Pl_squared=prefac * alpha21,
        inv_16piG_natural=1.0 / (16.0 * math.pi * prefac * alpha21),
    )


# ==========================================================================
# 2. Sakharov heat-kernel coefficient: per-DOF Lambda^2 contribution.
# ==========================================================================

# Per-bosonic-DOF coefficient for 1/(16 pi G) is Lambda^2 / (192 pi^2).
# This is the standard Sakharov result; see Visser 2002, "Sakharov's
# induced gravity: a modern perspective".

SAKHAROV_PER_DOF_COEFF = 1.0 / (192.0 * math.pi ** 2)   # times Lambda^2


def matter_dof_count() -> dict:
    """Effective DOF count for NWT matter (psi, A, n) at BPS vacuum.

    Returns the spin-resolved count and total DOF.
    """
    # 1 massive real radial Higgs (phi_1)
    # 1 massive Higgsed vector (3 polarisations after eating Goldstone phi_2)
    # 2 massless Goldstones (n_1, n_2 from S^2 Skyrme field)
    return dict(
        n_massive_scalar=1,            # phi_1 with m^2 = lambda v^2
        n_massive_vector=1,            # A_mu (3 DOF after Higgs)
        n_massless_scalar=2,           # n_1, n_2
        # Total bosonic DOF that contribute to Lambda^2 in Sakharov:
        N_DOF=1 + 3 + 2,               # = 6
    )


# ==========================================================================
# 3. Solve for the Lambda_UV that reproduces Paper 17's G.
# ==========================================================================

def solve_lambda_uv_natural(N_dof: int, target_inv_16piG: float
                              ) -> float:
    """Solve  N_dof * Lambda^2 / (192 pi^2)  =  target_inv_16piG
    for Lambda in natural units where m_e = 1.
    """
    return math.sqrt(target_inv_16piG * 192.0 * math.pi ** 2 / N_dof)


def lambda_uv_in_GeV(Lambda_natural: float) -> float:
    """Convert Lambda from m_e units to GeV."""
    return Lambda_natural * M_E_GEV


# ==========================================================================
# 4. Comparison to known NWT scales.
# ==========================================================================

KNOWN_SCALES_GEV = {
    'm_electron':   M_E_GEV,
    'GUT (NWT)':    E_GUT_GEV,
    'Planck':       M_PL_GEV,
    'm_e / alpha':  M_E_GEV / ALPHA,
    'm_e / alpha^2': M_E_GEV / ALPHA ** 2,
    'm_e * alpha^(-21/2)': M_E_GEV / ALPHA ** 10.5,
    'sqrt(8 pi/3) m_Pl (cosmological)': math.sqrt(8 * math.pi / 3) * M_PL_GEV,
}


def closest_scale(Lambda_GeV: float) -> tuple[str, float]:
    """Find the closest known scale (in log space)."""
    best = (None, float('inf'))
    for name, scale in KNOWN_SCALES_GEV.items():
        ratio = abs(math.log10(Lambda_GeV / scale))
        if ratio < best[1]:
            best = (name, ratio)
    return best


# ==========================================================================
# 5. Output.
# ==========================================================================

def section(title: str) -> None:
    print()
    print('=' * 76)
    print(f' {title}')
    print('=' * 76)


def main() -> None:
    section('Paper 18 -- G2: matter-loop 1/(16 pi G) via Sakharov')

    # ---- Paper 17 target ----
    section('Step 1: Paper 17 target for 1/(16 pi G)')
    p17 = paper17_G_natural()
    print(f"\n  Paper 17 structural prediction (LO + NLO bracket):")
    print(f"    (8/7)^2 (1 + alpha/7)^2  =  {p17['prefactor']:.10f}")
    print(f"    alpha^21                  =  {p17['alpha21']:.6e}")
    print(f"    G in natural units        =  {p17['G_natural']:.6e}")
    print(f"    1/(16 pi G) in nat units  =  {p17['inv_16piG_natural']:.6e}")
    print(f"    (m_e/m_Pl)^2              =  {p17['m_e_over_m_Pl_squared']:.6e}")
    print(f"    sqrt(m_e/m_Pl)^2 = m_e/m_Pl = {math.sqrt(p17['G_natural']):.6e}")
    print(f"    For comparison, CODATA m_e/m_Pl = 4.18549e-23")

    # ---- DOF counting ----
    section('Step 2: NWT matter content at BPS vacuum')
    dof = matter_dof_count()
    print(f"\n  After Higgs mechanism on the BPS vacuum:")
    print(f"    massive real scalar  (phi_1, Higgs radial)  : "
          f"{dof['n_massive_scalar']} DOF")
    print(f"    massive Higgsed vector (A_mu w/ eaten phi_2): "
          f"{dof['n_massive_vector']} DOFs (3 each)")
    print(f"    massless Goldstones (n_1, n_2 of S^2)       : "
          f"{dof['n_massless_scalar']} DOF")
    print(f"    --------------------")
    print(f"    Total bosonic DOF for Sakharov Lambda^2 coeff: "
          f"N_DOF = {dof['N_DOF']}")

    # ---- Sakharov result ----
    section('Step 3: Sakharov per-DOF coefficient')
    print(f"\n  Standard formula:")
    print(f"    Delta(1/(16 pi G))_per_DOF  =  Lambda_UV^2 / (192 pi^2)")
    print(f"                                =  {SAKHAROV_PER_DOF_COEFF:.6e} * "
          f"Lambda_UV^2")
    print(f"\n  Total NWT matter contribution (sum over species):")
    print(f"    Delta(1/(16 pi G))_total = {dof['N_DOF']} * "
          f"Lambda_UV^2 / (192 pi^2)")
    print(f"                              = "
          f"{dof['N_DOF'] * SAKHAROV_PER_DOF_COEFF:.6e} * Lambda_UV^2")

    # ---- Solve for required cutoff ----
    section('Step 4: solve for the Lambda_UV that reproduces Paper 17')
    Lambda_natural = solve_lambda_uv_natural(
        dof['N_DOF'], p17['inv_16piG_natural']
    )
    Lambda_GeV = lambda_uv_in_GeV(Lambda_natural)
    print(f"\n  Required Lambda_UV (in m_e units):  {Lambda_natural:.6e}")
    print(f"  Required Lambda_UV (in GeV):        {Lambda_GeV:.6e}")

    # ---- Compare to known scales ----
    section('Step 5: compare to known NWT scales')
    print(f"\n  Lambda_UV vs known scales:")
    print(f"  {'scale':>30} {'value (GeV)':>16} {'ratio':>10}")
    print('  ' + '-' * 60)
    for name, scale in KNOWN_SCALES_GEV.items():
        ratio = Lambda_GeV / scale
        print(f"  {name:>30} {scale:>16.4e} {ratio:>10.4e}")

    closest_name, log_dist = closest_scale(Lambda_GeV)
    print(f"\n  Closest known scale (in log):  {closest_name}  "
          f"(log10 distance {log_dist:.2f})")

    # ---- Interpretation ----
    section('Step 6: interpretation')
    print(f"""
  The Sakharov-induced 1/(16 pi G) from the NWT matter content
  matches Paper 17's structurally-derived value if and only if the UV
  cutoff is

      Lambda_UV  =  {Lambda_GeV:.3e} GeV
                 =  {Lambda_GeV / E_GUT_GEV:.3e} * E_GUT
                 =  {Lambda_GeV / M_PL_GEV:.3e} * M_Pl
                 =  {Lambda_natural:.3e} * m_e

  This number is the SCALE AT WHICH new physics has to enter to
  Wilson-renormalize the matter-loop contribution down from its naive
  Sakharov estimate.

  Naive Sakharov gives 1/G ~ N_DOF * Lambda^2 / (192 pi^2).
  Paper 17 gives          1/G  =  (8/7)^{{-2}} (1+alpha/7)^{{-2}} alpha^{{-21}} m_e^2 / (16 pi).
  Equating:               Lambda^2  =  (12 pi / N_DOF) / [(8/7)^2 (1+alpha/7)^2 alpha^21] * m_e^2.

  The ratio of Paper 17's G to naive Sakharov G with cutoff at any
  given scale is exactly the suppression that Paper 17's K_7
  Wilson-amplitude structure must supply.

  G3 will fold in K_7 amplitude factors and verify the suppression
  reproduces Paper 17 self-consistently.

  Per Phase 5's note (nwt_zeta_phase5_1overG.py): "naive Casimir
  framework computes the NAIVE coefficient without [Wilson-amplitude]
  suppression."  G2 quantifies that gap; G3 closes it.
""")


if __name__ == '__main__':
    main()
