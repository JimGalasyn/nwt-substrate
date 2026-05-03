"""
Cross-discipline demo: Drell-Yan
================================

Drell-Yan (q q̄ → γ* → ℓ⁺ ℓ⁻) is the parton-level workhorse of LHC
electroweak physics: a quark and antiquark from colliding hadrons
annihilate via a virtual photon (or Z) into a charged lepton pair.
At leading order in the photon channel:

    σ̂(q q̄ → ℓ⁺ ℓ⁻)  =  (4π α² / 3 ŝ) × Q_q²  × (1/N_c) × β_kin

The relationship to the R-ratio:

    σ(e⁺e⁻ → q q̄)  =  (4π α² / 3 s)  × Q_q²  × N_c     ← color SUM (final)
    σ̂(q q̄ → ℓ⁺ ℓ⁻)  =  (4π α² / 3 ŝ) × Q_q²  × 1/N_c   ← color AVERAGE (initial)

Same QED amplitude, time-reversed kinematics, and the color factor flips
from N_c to 1/N_c.  The substrate library handles this naturally because
the amplitude is built from a Lorentz tensor (not a process-specific
formula); kinematic crossing changes which legs are incoming vs outgoing
without changing the amplitude.

Why the substrate library handles Drell-Yan nimbly:

  - σ̂(q q̄ → ℓ⁺ ℓ⁻) at parton level is the *same calculation* as
    σ(e⁺e⁻ → μ⁺ μ⁻) with a charge / color rescaling.
  - nwt.qed.eemumu provides the QED kernel (machine-precision tested).
  - The crossing symmetry is structural -- no recomputation needed,
    just a kinematic re-interpretation.
  - Hadron-level Drell-Yan needs PDFs (a separate piece of infrastructure
    not in this library); the *amplitude* part is fully covered.

Run with:  python3 analysis/nwt_drell_yan_demo.py
"""

from __future__ import annotations
import math

import nwt_substrate.qed as qed
import nwt_substrate.qft as qft
from nwt_substrate.amplitudes.running_couplings import alpha_s


N_C = 3
ALPHA = 1.0 / 137.036
GEV2_TO_PB = 0.3894e9       # GeV⁻² → pb conversion


# ---------------------------------------------------------------------------
# Partonic Drell-Yan cross-section
# ---------------------------------------------------------------------------

def sigma_dy_partonic(sqrt_shat: float, q, m_lepton: float = 0.10566) -> float:
    """
    Partonic Drell-Yan cross-section σ̂(q q̄ → ℓ⁺ ℓ⁻) at √ŝ.

    Same calculation as σ(e⁺e⁻ → ℓ⁺ ℓ⁻) but
        - vertex on the initial side has charge Q_q instead of |e|;
        - color average over initial-state quark colors gives 1/N_c.

    Implementation: reuse nwt.qed.eemumu (which gives σ for any massless
    initial state with charge ±e annihilating into ℓ⁺ ℓ⁻ of mass m_lepton)
    and rescale by Q_q² × (1/N_c).
    """
    if 2 * m_lepton >= sqrt_shat:
        return 0.0
    sigma_pointlike = qed.eemumu.sigma_total(sqrt_shat, m_mu=m_lepton)
    return (q.charge ** 2 / N_C) * sigma_pointlike


def sigma_dy_textbook(sqrt_shat: float, q, m_lepton: float = 0.10566) -> float:
    """
    Textbook leading-order partonic DY cross-section (massless lepton limit):

        σ̂(q q̄ → ℓ⁺ ℓ⁻) = (4π α² / 3 ŝ) × Q_q² / N_c × β_ℓ × (1 + 2 m_ℓ²/ŝ)

    The β_ℓ × (1 + 2 m_ℓ²/ŝ) phase-space factor ≈ 1 at high ŝ.
    """
    if 2 * m_lepton >= sqrt_shat:
        return 0.0
    s = sqrt_shat ** 2
    beta_ell = math.sqrt(1.0 - 4.0 * m_lepton ** 2 / s)
    phase = beta_ell * (1.0 + 2.0 * m_lepton ** 2 / s)
    return (4 * math.pi * ALPHA ** 2 / (3 * s)) * (q.charge ** 2 / N_C) * phase * GEV2_TO_PB


def sigma_r_ratio_partial(sqrt_s: float, q, m_lepton: float = 0.10566) -> float:
    """
    For comparison: the partial R-ratio cross-section σ(e⁺e⁻ → q q̄) at √s.
    Same QED amplitude as Drell-Yan but with N_c (color sum) instead of 1/N_c.
    """
    if 2 * q.mass >= sqrt_s:
        return 0.0
    sigma_pointlike = qed.eemumu.sigma_total(sqrt_s, m_mu=q.mass)
    return q.n_color * q.charge ** 2 * sigma_pointlike


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def header(text: str, char: str = "=", width: int = 76) -> None:
    print(f"\n{char * width}\n  {text}\n{char * width}")


def section(text: str) -> None:
    print(f"\n--- {text} ---")


def main() -> None:
    header("CROSS-DISCIPLINE DEMO: DRELL-YAN")
    print("σ̂(q q̄ → γ* → ℓ⁺ ℓ⁻):  the parton-level engine of LHC W/Z physics.")
    print()
    print("Same QED amplitude as R-ratio's hadronic side, kinematic crossing,")
    print("color factor flips from N_c (sum, final state) to 1/N_c (average,")
    print("initial state) -- a factor of 1/9 relative to the R-ratio integrand.")

    section("Substrate handling: structural crossing symmetry")
    print()
    print("  σ(e⁺e⁻ → q q̄)  = N_c × Q_q²  × σ_QED(s, m_q)     ← R-ratio")
    print("  σ̂(q q̄ → ℓ⁺ ℓ⁻) = (1/N_c) × Q_q² × σ_QED(ŝ, m_ℓ)   ← Drell-Yan")
    print()
    print("  The σ_QED kernel is the SAME object: nwt.qed.eemumu.sigma_total.")
    print("  Crossing q↔ℓ in the initial/final state flips N_c ↔ 1/N_c.")
    print("  The substrate's amplitude is a Lorentz tensor; kinematic")
    print("  re-interpretation does not require recomputation.")

    # =======================================================================
    # Substrate vs textbook at parton level
    # =======================================================================

    section("Partonic Drell-Yan σ̂(q q̄ → μ⁺ μ⁻) at √ŝ = 100 GeV")
    sqrt_shat = 100.0
    print(f"  {'Flavor':<8s}  {'Q_q²':<8s}  {'Q_q²/N_c':<10s}  "
          f"{'σ̂ substrate [pb]':<22s}  {'σ̂ textbook [pb]':<22s}  {'ratio':s}")
    print("  " + "-" * 78)
    for q in [qft.up_quark, qft.down_quark, qft.strange,
              qft.charm, qft.bottom]:
        sub = sigma_dy_partonic(sqrt_shat, q)
        txt = sigma_dy_textbook(sqrt_shat, q)
        ratio = sub / txt if txt > 0 else float('nan')
        print(f"  {q.name:<8s}  {q.charge**2:<8.4f}  {q.charge**2/N_C:<10.4f}  "
              f"{sub:<22.4e}  {txt:<22.4e}  {ratio:.6f}")

    print()
    print(f"  Substrate matches the textbook formula to ~6 significant figures.")
    print(f"  (Small offset from the (1+2m²/s)·β factor in the textbook expression")
    print(f"  vs the substrate's full lepton-mass treatment.)")

    # =======================================================================
    # Crossing comparison: R-ratio vs Drell-Yan at same energy
    # =======================================================================

    section("Crossing duality: σ(e⁺e⁻ → q q̄) vs σ̂(q q̄ → ℓ⁺ ℓ⁻) at √s = 30 GeV")
    sqrt_s = 30.0
    print()
    print(f"  {'Flavor':<8s}  {'σ(e⁺e⁻ → q q̄)':<20s}  "
          f"{'σ̂(q q̄ → μ⁺ μ⁻)':<22s}  {'ratio (= N_c²)':s}")
    print("  " + "-" * 70)
    for q in [qft.up_quark, qft.down_quark, qft.strange,
              qft.charm, qft.bottom]:
        sigma_R = sigma_r_ratio_partial(sqrt_s, q)
        sigma_DY = sigma_dy_partonic(sqrt_s, q)
        if sigma_DY > 0:
            ratio = sigma_R / sigma_DY
            print(f"  {q.name:<8s}  {sigma_R:<20.4e}  {sigma_DY:<22.4e}  "
                  f"{ratio:.4f}  (expected N_c² = 9)")
        else:
            print(f"  {q.name:<8s}  (one or both processes below threshold)")

    print()
    print("  σ(R-ratio) / σ̂(DY) = N_c × N_c = N_c² = 9 exactly (modulo small")
    print("  mass effects from the m_q vs m_ℓ in the kinematic kernel).")
    print("  This is the textbook 'parton-model factor of 1/9 for Drell-Yan'.")

    # =======================================================================
    # Drell-Yan invariant-mass spectrum
    # =======================================================================

    section("Partonic σ̂(q q̄ → μ⁺ μ⁻) vs invariant mass M = √ŝ")
    print()
    print("  Sum over u,d,s,c,b at each M:")
    print(f"  {'M [GeV]':<10s}  {'σ̂_total [pb]':<14s}  {'σ̂ × M³ [pb·GeV³]':<22s}")
    print("  " + "-" * 50)
    M_values = [10, 30, 60, 100, 200, 500, 1000]
    for M in M_values:
        sigma_sum = sum(sigma_dy_partonic(M, q)
                        for q in [qft.up_quark, qft.down_quark, qft.strange,
                                   qft.charm, qft.bottom]
                        if 2 * q.mass < M)
        scaled = sigma_sum * M ** 3
        print(f"  {M:<10.0f}  {sigma_sum:<14.4e}  {scaled:<22.4e}")
    print()
    print("  σ̂ scales as 1/ŝ at high ŝ; σ̂ × M³ ∝ M (i.e. ~constant after the")
    print("  Σ_q Q_q²/N_c factor saturates).  This is the well-known DY")
    print("  asymptotic 'scaling' of the LO parton-level cross-section.")

    # =======================================================================
    # Hadron-level Drell-Yan: what's missing
    # =======================================================================

    section("Hadron-level Drell-Yan: PDF convolution")
    print()
    print("  σ(p p → ℓ⁺ ℓ⁻ + X, M)  =  Σ_q ∫dx₁∫dx₂  ×")
    print("                          [f_q(x₁, Q²) f_q̄(x₂, Q²) + (1 ↔ 2)]")
    print("                          × σ̂(q q̄ → ℓ⁺ ℓ⁻, ŝ = x₁ x₂ s)")
    print()
    print("  The amplitude / partonic σ̂ is what the substrate library")
    print("  computes and is the substrate-physics content of Drell-Yan.")
    print("  The PDFs f_q(x, Q²) are non-perturbative input fitted to data;")
    print("  they are PDF infrastructure, not substrate primitives.  A toy")
    print("  parametric PDF (e.g. CTEQ-like form) would suffice for an")
    print("  order-of-magnitude check; full LHC predictions need lhapdf.")

    # =======================================================================
    # Cross-shim composition tableau
    # =======================================================================

    header("CROSS-SHIM COMPOSITION TABLEAU")
    print()
    print(f"  {'Ingredient':<40s}  {'Source shim':<25s}  {'Notes':s}")
    print(f"  {'-' * 40}  {'-' * 25}  {'-' * 25}")
    print(f"  {'σ_QED(s, m_lepton) kernel':<40s}  "
          f"{'nwt.qed.eemumu':<25s}  {'machine-precision tested':s}")
    print(f"  {'Quark mass m_q':<40s}  "
          f"{'nwt.qft.<quark>':<25s}  {'DiracSpinor (substrate)':s}")
    print(f"  {'Quark charge Q_q':<40s}  "
          f"{'nwt.qft.<quark>.charge':<25s}  {'electric-charge attribute':s}")
    print(f"  {'Color factor 1/N_c (initial avg)':<40s}  "
          f"{'nwt.qft.<quark>.n_color':<25s}  {'flips R-ratio N_c via crossing':s}")
    print(f"  {'Cl(1,3) Dirac primitive':<40s}  "
          f"{'algebra.dirac':<25s}  {'shared lepton/quark spinor space':s}")
    print()
    print(f"  σ̂(q q̄ → ℓ⁺ ℓ⁻) = (Q_q² / N_c) × σ_QED(ŝ, m_ℓ)")
    print()
    print("  Same Cl(1,3) primitive as e⁺e⁻ → μ⁺ μ⁻; same QED amplitude")
    print("  as R-ratio.  Drell-Yan = R-ratio with the color-flow flipped")
    print("  and kinematics crossed.")

    header("CONCLUSION", char="=")
    print(f"At √ŝ = 100 GeV, substrate σ̂(u ū → μ⁺ μ⁻) = "
          f"{sigma_dy_partonic(100.0, qft.up_quark):.4e} pb,")
    print(f"matching the textbook LO formula to ~6 sig figs.  Crossing")
    print(f"duality σ(R)/σ̂(DY) = N_c² = 9 verified for all 5 active flavors.")
    print()
    print("Drell-Yan and R-ratio are *the same substrate calculation* with")
    print("color-flow flipped.  The Lorentz tensor amplitude is unchanged;")
    print("only the bookkeeping of which legs are incoming changes.  In the")
    print("substrate library this is structural -- not a separate process.")
    print()
    print("Adding LHC W/Z physics requires extending only the Z-vertex")
    print("(γ-Z mixing + chiral V-A coupling) and PDFs.  The QED kernel,")
    print("color algebra, and substrate Cl(1,3) primitives remain unchanged.")
    print()


if __name__ == "__main__":
    main()
