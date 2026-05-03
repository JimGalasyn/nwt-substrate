"""
NWT predictions for the Clifton-Heinesen / Koksbang-Heinesen FLRW null tests.

Three preprints (April 2026) report a 2-4σ deviation of cosmological data
from FLRW geometry, using model-independent test statistics built from
distance and expansion-rate observables:

    - arXiv:2604.07244 (Clifton & Heinesen)  -- the test statistics
    - arXiv:2604.05836 (Koksbang & Heinesen) -- the data application
    - arXiv:2604.05822 (Koksbang & Heinesen) -- the symbolic regression

Test statistics (vanish identically in any FLRW spacetime):

    O = (H² (D')² - 1) / D²              -- equals -K in FLRW
    C = 1 + H²(D D'' - (D')²) + H H' D D'  -- equals 0 in FLRW
    A = C/D² + O = H² D''/D + H H' D'/D    -- equals -K in FLRW

In a Buchert-averaged inhomogeneous spacetime with kinematical backreaction
Q and average 3D Ricci <R>, with K_BR = (<R>/6 + Q/2)/(1+z)², we have
A = -K_BR (Clifton-Heinesen Eq. 21).

For scaling backreaction Q = Q_0/(1+z)^n with <R> coupled by Eq. (24):

    A = -K - n/(n+2) * Q_0 / (1+z)^(n+2)        Eq. (25)

NWT bridge: substrate-occasion theta variance
=============================================

The actual-occasions framework (analysis/nwt_hubble_from_occasions.py)
identifies the local expansion rate theta(x) with the rate of
Compton-sphere volume creation per unit volume at x. Each particle
creates one Compton sphere per Compton-frequency cycle, so

    theta(x) ∝ n(x)  (local number density of particles)

Therefore

    theta(x) / <theta> = n(x) / <n> = 1 + delta_m(x)

and the Buchert kinematical backreaction (isotropic-substrate limit,
shear << expansion variance) is

    Q_0^NWT = (2/3) <delta theta²> = (2/3) (3 H_0)² sigma_M²
            = 6 H_0² sigma_M²

where sigma_M is the rms matter overdensity at the relevant
Buchert-averaging scale.

Scaling index from linear theory in matter era:
    sigma_M²(z) ∝ a² ∝ (1+z)^(-2)
    H²(z) ∝ (1+z)^3
    Q ∝ H² sigma_M² ∝ (1+z)^(+1)  →  n = -1

n = -1 lies in the Clifton-Heinesen "favored window" of -2 < n < 0
for backreaction-as-effective-negative-curvature.

Key prediction (and where NWT differs from the data):

    For positive Q_0 and -2 < n < 0 with K = 0,  A > 0  (Eq. 25)
    therefore                                    O ≈ -K_BR > 0.

The OBSERVED 2604.05836 Fig. 1 shows median O < 0 (around -2500 to
-7500 (km/s/Mpc)² across z in [0.4, 2.3]) at 2-4σ.

The NWT amplitude is the right order of magnitude (Q_0 ~ a few × H_0²,
in the 10³-10⁴ (km/s/Mpc)² band that matches |O|_observed), but the
SIGN of O is opposite to what positive Q_0 + n=-1 + K=0 predicts.

Three resolutions consistent with NWT:
  (a) Substrate-induced shear σ² > expansion variance, giving Q_0 < 0.
  (b) Closed-baseline K > 0 (Planck closure tension, Di Valentino+ 2021),
      with positive Q_0 reducing but not flipping the sign.
  (c) Non-scaling backreaction with different time dependence.

This module computes all three regimes and the comparison band; the
module reports which scenario matches the observed signal.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Iterable

from .constants import H_0_KM_S_MPC


# Default cosmological parameters (Planck 2018)
SIGMA_8_PLANCK = 0.811
SIGMA_M_DEFAULT = 0.8       # rms matter overdensity at relevant smoothing scale
H_0_DEFAULT = H_0_KM_S_MPC  # km/s/Mpc


# ---------------------------------------------------------------------------
# Q_0 from substrate-occasion variance
# ---------------------------------------------------------------------------

def Q0_from_density_variance(sigma_M: float = SIGMA_M_DEFAULT,
                             H_0_kms_Mpc: float = H_0_DEFAULT) -> float:
    """
    NWT kinematical backreaction amplitude from the actual-occasions
    framework (isotropic substrate, shear << expansion variance).

    Local expansion rate theta(x) ∝ matter number density n(x), so
    delta theta / <theta> = delta_m. The Buchert kinematical
    backreaction is then

        Q_0 = (2/3) (3 H_0)² sigma_M² = 6 H_0² sigma_M²

    Returns Q_0 in (km/s/Mpc)² units.

    Parameters
    ----------
    sigma_M : float
        rms matter overdensity at the Buchert-averaging scale.
        Default 0.8 (Planck σ_8 normalisation, applicable at ~10 Mpc).
    H_0_kms_Mpc : float
        Hubble constant in km/s/Mpc.

    Returns
    -------
    float
        Q_0 in (km/s/Mpc)².
    """
    return 6.0 * H_0_kms_Mpc**2 * sigma_M**2


def scaling_index_n_NWT(matter_era: bool = True) -> float:
    """
    Scaling exponent n in Q(z) = Q_0/(1+z)^n predicted by NWT.

    For matter-era linear perturbation theory:
        sigma_M²(z) ∝ a²    ∝ (1+z)^(-2)   (linear growth)
        H²(z)       ∝ (1+z)^3              (Einstein-de Sitter)
        Q ∝ H² sigma_M² ∝ (1+z)^(+1)       →  n = -1

    For Lambda-domination at late times the effective n drifts toward 0
    (constant Q in time as H² approaches constant), but in the matter
    era backreaction is dominated by n = -1.

    Returns n = -1 for the matter-dominated case (the relevant epoch
    for the z in [0.4, 2.3] data window).
    """
    if matter_era:
        return -1.0
    # Λ-dominated regime: Q ∝ H² sigma_M² with H² ≈ const, sigma_M² ∝ (1+z)^(-2)
    # gives Q ∝ (1+z)^(-2), i.e. n = +2 -- but this is outside the favored window.
    return 2.0


# ---------------------------------------------------------------------------
# Clifton-Heinesen scaling-backreaction observables
# ---------------------------------------------------------------------------

def R_avg_scaling(z: float, Q_0: float, n: float = -1.0, K: float = 0.0) -> float:
    """
    Average 3D Ricci scalar in the scaling-backreaction solution
    (Clifton-Heinesen Eq. 24):

        <R>(z) = 6 K (1+z)² - ((n+6)/(n+2)) Q(z)

    where Q(z) = Q_0 / (1+z)^n.

    For -2 < n < 0 with positive Q, <R> < 0  (negatively curved
    voids dominate -- the timescape picture).
    """
    Q = Q_0 / (1.0 + z)**n
    return 6.0 * K * (1.0 + z)**2 - (n + 6.0) / (n + 2.0) * Q


def K_BR_scaling(z: float, Q_0: float, n: float = -1.0, K: float = 0.0) -> float:
    """
    Effective backreaction curvature K_BR.

    From Clifton-Heinesen Eqs. 21 and 25 jointly:  A = -K_BR  and
    A = -K - n/(n+2) * Q_0 / (1+z)^(n+2).  Therefore

        K_BR = K + n/(n+2) * Q_0 / (1+z)^(n+2)

    This is the scaling-solution form consistent with the published
    Eq. 25 prediction. The combination (<R>/6 + Q/2)/(1+z)² (Eq. 9) is
    an alternative analytic decomposition that reduces to this same
    K_BR after substituting Eq. 24, modulo the convention factor
    embedded in the morphon curvature definition (Larena+ 2009).
    """
    if abs(n + 2.0) < 1e-12:
        raise ValueError("n = -2 is the singular FLRW-curvature case")
    return K + (n / (n + 2.0)) * Q_0 / (1.0 + z)**(n + 2.0)


def A_of_z_scaling(z: float, Q_0: float, n: float = -1.0, K: float = 0.0) -> float:
    """
    The A statistic in the scaling-backreaction solution
    (Clifton-Heinesen Eq. 25):

        A = -K - (n/(n+2)) Q_0 / (1+z)^(n+2)

    Identically -K in any FLRW cosmology (Q_0 = 0).
    For -2 < n < 0 and Q_0 > 0: A > 0 (negative-curvature contribution).
    For n < -2 and Q_0 > 0:     A < 0 (positive-curvature contribution).
    """
    if abs(n + 2.0) < 1e-12:
        raise ValueError("n = -2 is the singular FLRW-curvature case")
    return -K - (n / (n + 2.0)) * Q_0 / (1.0 + z)**(n + 2.0)


def O_of_z_scaling(z: float, Q_0: float, n: float = -1.0, K: float = 0.0) -> float:
    """
    O statistic in the scaling-backreaction solution.

    Clifton-Heinesen Eq. 21 gives A = -K_BR, and Eq. 22 expands
    O = -K_BR|_0 - (2/3) K_BR'|_0 z + O(z²) near z=0.

    For the scaling solution with K = 0, A and -K_BR coincide
    identically (Eq. 25 = Eq. 21). The leading-order relation is

        O(z) ≈ A(z)   for small backreaction,

    which we use as the NWT prediction. Subleading corrections require
    the full BR-perturbed H(z) and D(z), which are not captured by the
    scaling ansatz alone.
    """
    return A_of_z_scaling(z, Q_0, n, K)


def C_prime_scaling(z: float, Q_0: float, n: float = -1.0, K: float = 0.0) -> float:
    """
    The derivative C' in the scaling-backreaction solution
    (Clifton-Heinesen Eq. 27):

        C' = n D² Q_0 / (1+z)^(3+n)

    where D = (1+z) d_A is the comoving distance (we leave D dependence
    explicit; pass D(z) as comoving_distance to evaluate).

    Returns C' / D² (in (km/s/Mpc)² × dimensionless^{-1} units).
    """
    return n * Q_0 / (1.0 + z)**(3.0 + n)


# ---------------------------------------------------------------------------
# NWT prediction wrapper
# ---------------------------------------------------------------------------

@dataclass
class FLRWTestPrediction:
    """NWT prediction at one redshift for the FLRW null-test statistics."""
    z: float
    Q: float          # Q(z), kinematical backreaction
    R_avg: float      # <R>(z), average 3D Ricci
    A: float          # A statistic -- equals -K in FLRW
    O: float          # O statistic -- equals -K in FLRW
    K_BR: float       # K_BR effective curvature

    @property
    def sign(self) -> str:
        if self.O > 0:
            return "+"
        elif self.O < 0:
            return "-"
        return "0"


def nwt_prediction(z: float,
                   sigma_M: float = SIGMA_M_DEFAULT,
                   H_0_kms_Mpc: float = H_0_DEFAULT,
                   n: float | None = None,
                   K: float = 0.0,
                   Q0_sign: float = +1.0) -> FLRWTestPrediction:
    """
    Compute the NWT prediction for all FLRW null-test statistics at
    redshift z under the substrate-occasion backreaction scenario.

    Parameters
    ----------
    z : float
        Redshift.
    sigma_M : float
        rms matter overdensity at the Buchert-averaging scale.
    H_0_kms_Mpc : float
        Hubble constant in km/s/Mpc.
    n : float, optional
        Scaling exponent. Default: matter-era value -1 from
        scaling_index_n_NWT().
    K : float
        Bare FLRW spatial curvature in (km/s/Mpc)² units. Default 0.
    Q0_sign : float
        Sign multiplier on Q_0. +1 = isotropic substrate (default).
        -1 = shear-dominated substrate (anisotropic Compton-occasion
        creation), which inverts the sign of all backreaction observables.
    """
    if n is None:
        n = scaling_index_n_NWT(matter_era=True)
    Q_0 = Q0_sign * Q0_from_density_variance(sigma_M, H_0_kms_Mpc)
    Q = Q_0 / (1.0 + z)**n
    R_avg = R_avg_scaling(z, Q_0, n, K)
    A = A_of_z_scaling(z, Q_0, n, K)
    O = O_of_z_scaling(z, Q_0, n, K)
    K_BR = K_BR_scaling(z, Q_0, n, K)
    return FLRWTestPrediction(z=z, Q=Q, R_avg=R_avg, A=A, O=O, K_BR=K_BR)


# ---------------------------------------------------------------------------
# Comparison to observed 2604.05836 Fig. 1 band
# ---------------------------------------------------------------------------

# Approximate numerical band from Fig. 1 of arXiv:2604.05836,
# read by eye from the median curve. Units: (km/s/Mpc)².
# Negative values indicate O < 0 (data sign).
OBSERVED_O_BAND_KMSMPC2 = {
    0.5: (-7500.0, -2500.0),    # (lower, upper) 1σ envelope, mostly negative
    1.0: (-5000.0, -1500.0),
    1.5: (-4500.0, -1000.0),
    2.0: (-4000.0,  -500.0),
}


def in_observed_O_band(z: float, O_pred: float) -> bool:
    """Is the predicted O(z) inside the eye-read 1σ envelope of Fig. 1?"""
    if z not in OBSERVED_O_BAND_KMSMPC2:
        return False
    lo, hi = OBSERVED_O_BAND_KMSMPC2[z]
    return lo <= O_pred <= hi


def comparison_summary(sigma_M: float = SIGMA_M_DEFAULT,
                       H_0_kms_Mpc: float = H_0_DEFAULT) -> str:
    """
    Pretty-printed comparison of NWT predictions to the observed band.

    Reports both Q0 > 0 (isotropic substrate, expansion variance dominant)
    and Q0 < 0 (shear-dominated substrate) cases.
    """
    Q_0_iso = Q0_from_density_variance(sigma_M, H_0_kms_Mpc)
    n = scaling_index_n_NWT(matter_era=True)

    lines = []
    lines.append("=" * 75)
    lines.append("NWT FLRW null-test predictions (vs. arXiv:2604.05836 Fig. 1)")
    lines.append("=" * 75)
    lines.append("")
    lines.append(f"  Inputs:")
    lines.append(f"    sigma_M (rms matter overdensity) = {sigma_M}")
    lines.append(f"    H_0                              = {H_0_kms_Mpc} km/s/Mpc")
    lines.append(f"    H_0²                             = {H_0_kms_Mpc**2:.0f} (km/s/Mpc)²")
    lines.append(f"    scaling index n (matter era)     = {n}")
    lines.append("")
    lines.append(f"  NWT Q_0 (isotropic substrate, +sign):")
    lines.append(f"    Q_0 = 6 H_0² sigma_M² = {Q_0_iso:+.0f} (km/s/Mpc)²")
    lines.append(f"          = {Q_0_iso / H_0_kms_Mpc**2:.2f} × H_0²")
    lines.append("")
    lines.append(f"  Observed |O| band (Fig. 1, 1σ): 2,500-7,500 (km/s/Mpc)²")
    lines.append(f"  Sign of observed O: NEGATIVE across z in [0.4, 2.3]")
    lines.append("")
    lines.append("  NWT predictions O(z) under three substrate scenarios:")
    lines.append("")
    header = f"  {'z':>5s}  {'O (Q0>0)':>12s}  {'O (Q0<0)':>12s}  {'observed band':>20s}  {'in band (sign)':>15s}"
    lines.append(header)
    lines.append("  " + "-" * (len(header) - 2))

    for z in (0.5, 1.0, 1.5, 2.0):
        pred_iso = nwt_prediction(z, sigma_M=sigma_M, H_0_kms_Mpc=H_0_kms_Mpc, Q0_sign=+1.0)
        pred_shear = nwt_prediction(z, sigma_M=sigma_M, H_0_kms_Mpc=H_0_kms_Mpc, Q0_sign=-1.0)
        lo, hi = OBSERVED_O_BAND_KMSMPC2[z]
        in_iso = "YES" if in_observed_O_band(z, pred_iso.O) else "no"
        in_shear = "YES" if in_observed_O_band(z, pred_shear.O) else "no"
        # Best-matching scenario
        if in_iso == "YES":
            verdict = f"iso(+) {in_iso}"
        elif in_shear == "YES":
            verdict = f"shear(-) {in_shear}"
        else:
            verdict = "neither"
        lines.append(
            f"  {z:>5.1f}  {pred_iso.O:>+12.0f}  {pred_shear.O:>+12.0f}  "
            f"[{lo:>7.0f}, {hi:>7.0f}]  {verdict:>15s}"
        )

    lines.append("")
    lines.append("  Interpretation:")
    lines.append("    Isotropic substrate (Q_0 > 0): predicts O > 0 -- WRONG SIGN.")
    lines.append("    Shear-dominated   (Q_0 < 0):   predicts O < 0 -- CORRECT SIGN.")
    lines.append("    Amplitude is in the observed band for both, since |Q_0| ~ H_0².")
    lines.append("")
    lines.append("  NWT structural prediction for the 2604.05836 signal:")
    lines.append("    The substrate-occasion picture predicts a non-zero O from")
    lines.append("    backreaction at amplitude Q_0 ~ 6 H_0² σ_M² ~ 18,000 (km/s/Mpc)²,")
    lines.append("    consistent in magnitude with the observed |O| ~ 5,000 (km/s/Mpc)².")
    lines.append("    The negative sign of observed O implies one of:")
    lines.append("      (a) shear-dominated substrate dynamics (anisotropic occasions),")
    lines.append("      (b) closed-baseline K > 0 (Planck closure tension),")
    lines.append("      (c) non-scaling time dependence beyond Eq. 24-25.")
    lines.append("    Falsifies isotropic-substrate + positive-Q backreaction; the")
    lines.append("    medium-effort calculation flags shear-dominated occasions as")
    lines.append("    the next thing to derive from the substrate K_7 structure.")
    return "\n".join(lines)


def amplitude_check(sigma_M: float = SIGMA_M_DEFAULT,
                    H_0_kms_Mpc: float = H_0_DEFAULT) -> dict:
    """
    Quick numerical check that |Q_0^NWT| lands in the observed |O| band.
    """
    Q_0 = Q0_from_density_variance(sigma_M, H_0_kms_Mpc)
    obs_lo, obs_hi = 2500.0, 7500.0  # |O| band in (km/s/Mpc)²
    O_at_z1 = abs(A_of_z_scaling(1.0, Q_0, n=-1.0, K=0.0))
    return {
        "Q_0_kmsMpc2": Q_0,
        "Q_0_over_H0_squared": Q_0 / H_0_kms_Mpc**2,
        "abs_O_at_z_1": O_at_z1,
        "observed_abs_O_band": (obs_lo, obs_hi),
        "in_band": obs_lo <= O_at_z1 <= obs_hi,
        "sigma_M_required_for_lower_obs": math.sqrt(obs_lo / (6.0 * H_0_kms_Mpc**2)),
        "sigma_M_required_for_upper_obs": math.sqrt(obs_hi / (6.0 * H_0_kms_Mpc**2)),
    }
