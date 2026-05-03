"""
Cross-discipline demo: the Z resonance
======================================

Completes the cross-community DY trilogy (Phases Q.8 R-ratio, Q.9 partonic
DY, Q.10 hadronic DY) by adding the Z boson via the new nwt.electroweak
shim.  Produces:

  1. The Z partial-width table (predicts Γ_Z = 2.42 GeV vs PDG 2.495 GeV)
  2. The Z resonance lineshape σ(e⁺e⁻ → μ⁺μ⁻) at peak ≈ 2 nb (PDG match)
  3. The famous γ-Z interference dip below the pole
  4. R-ratio updated with Z contribution at the Z peak (~10⁴-fold enhancement)

Substrate connection: the Cl(0,7) internal SU(2) (bivector triplet) is
vector under Lorentz; the SM's chiral SU(2)_L is implemented at the
vertex-factor level via V-A couplings (g_V = T_3 - 2Q sin²θ_W, g_A = T_3).
The library's QED amplitudes (Compton, eemumu, etc.) keep using α(0); the
electroweak shim uses the G_F-derived α(M_Z) ≈ 1/127.9 which absorbs the
QED running 0 → M_Z and gives PDG-matching Z observables.

Run with:  python3 analysis/nwt_z_resonance_demo.py
"""

from __future__ import annotations
import math

import nwt_substrate.electroweak as ew


def header(text: str, char: str = "=", width: int = 76) -> None:
    print(f"\n{char * width}\n  {text}\n{char * width}")


def section(text: str) -> None:
    print(f"\n--- {text} ---")


def main() -> None:
    header("CROSS-DISCIPLINE DEMO: THE Z RESONANCE")
    print("Adds the Z boson (γ-Z interference, chiral V-A couplings) to the")
    print("R-ratio / Drell-Yan picture from Phases Q.8-Q.10.")

    section("Substrate constants and EW couplings")
    print(f"  M_Z          = {ew.M_Z} GeV          (PDG)")
    print(f"  Γ_Z          = {ew.GAMMA_Z} GeV          (PDG)")
    print(f"  M_W          = {ew.M_W} GeV          (PDG)")
    print(f"  sin²θ_W      = {ew.SIN2_THETA_W:.5f}        (effective)")
    print(f"  G_F          = 1.1664 × 10⁻⁵ GeV⁻²  (Fermi constant, PDG)")
    print(f"  g_Z²         = 4√2 G_F M_Z² = {ew.G_Z**2:.5f}")
    print(f"  α(M_Z)       = {ew.E_CHARGE_Z**2/(4*math.pi):.6f} = 1/{4*math.pi/ew.E_CHARGE_Z**2:.2f}")
    print(f"  α(0) [QED]   = {ew.ALPHA_QED:.6f} = 1/{1/ew.ALPHA_QED:.3f}")
    print()
    print("  → α(M_Z) ≈ 1/128 includes the QED running 0 → M_Z (~6% enhancement);")
    print("    the substrate's nwt.qed amplitudes use α(0)=1/137 for low-energy")
    print("    Compton / e+e-→μμ; the electroweak shim uses the G_F-derived")
    print("    coupling for Z-pole physics.")

    section("SM fermion Z-couplings (g_V = T_3 - 2 Q sin²θ_W, g_A = T_3)")
    print(ew.coupling_summary())

    section("Z partial widths (LO)")
    print(ew.width_summary())

    section("Z resonance lineshape: σ(e⁺e⁻ → μ⁺μ⁻) vs √s")
    print()
    print(f"  {'√s [GeV]':<10s}  {'σ(γ+Z) [pb]':<14s}  "
          f"{'σ(γ-only) [pb]':<16s}  {'ratio':<10s}  {'comment':s}")
    print("  " + "-" * 76)
    for sqrt_s, comment in [
        (40,  "low-energy QED dominant"),
        (60,  "approaching Z pole"),
        (80,  "destructive γ-Z interference (the 'dip')"),
        (88,  "rising to Z pole"),
        (90,  "near peak"),
        (91.2,"Z PEAK -- PDG ≈ 2 nb"),
        (92,  "just past peak"),
        (94,  "falling from peak"),
        (100, "above pole; constructive interference"),
        (200, "high-energy continuum"),
        (500, "asymptotic regime"),
    ]:
        sigma_full = ew.sigma_total(sqrt_s, "mu", m_f=ew.FERMION_MASS_GEV["mu"])
        sigma_qed = ew.sigma_qed_only(sqrt_s, "mu", m_f=ew.FERMION_MASS_GEV["mu"])
        ratio = sigma_full / sigma_qed if sigma_qed > 0 else float('inf')
        print(f"  {sqrt_s:<10.1f}  {sigma_full:<14.3e}  {sigma_qed:<16.3e}  "
              f"{ratio:<10.3f}  {comment:s}")

    section("Z peak vs PDG textbook value")
    sigma_peak_substrate = ew.sigma_total(ew.M_Z, "mu",
                                           m_f=ew.FERMION_MASS_GEV["mu"])
    # PDG narrow-resonance approximation: σ_peak = 12π × Γ_ee × Γ_μμ / (M_Z² Γ_Z²)
    Gamma_ee = ew.partial_width_Z("e")
    Gamma_mumu = ew.partial_width_Z("mu")
    sigma_peak_BW = (12.0 * math.pi * Gamma_ee * Gamma_mumu
                     / (ew.M_Z ** 2 * ew.GAMMA_Z ** 2)) * 0.3894e9
    print(f"  Substrate σ_peak (full lineshape):  {sigma_peak_substrate:.0f} pb")
    print(f"  Narrow-resonance (12π Γee Γμμ / M²Γ²): {sigma_peak_BW:.0f} pb")
    print(f"  PDG observed σ_peak:                 ≈ 2000 pb (2 nb)")
    print(f"  agreement (substrate / PDG):         {sigma_peak_substrate/2000.0*100:.1f}%")

    section("R-ratio at the Z peak (with Z contribution)")
    print()
    print("  At √s ≪ M_Z and √s ≫ M_Z, the photon channel dominates and")
    print("  R ≈ N_c Σ Q_q² ≈ 11/3 (Phase Q.8).  At the Z pole, Z exchange")
    print("  enhances R by ~10⁴ because both σ_had and σ_μμ get the same")
    print("  Z-resonance enhancement (so R itself stays roughly constant)")
    print("  but the absolute cross-sections explode.")
    print()
    print(f"  {'√s [GeV]':<10s}  {'R_full(Z)':<12s}  {'R_QED-only':<12s}  "
          f"{'σ_had/σ_μμ':s}")
    print("  " + "-" * 70)
    for sqrt_s in [30, 60, 91.2, 200]:
        sigma_mumu_full = ew.sigma_total(sqrt_s, "mu",
                                          m_f=ew.FERMION_MASS_GEV["mu"])
        sigma_mumu_qed = ew.sigma_qed_only(sqrt_s, "mu",
                                            m_f=ew.FERMION_MASS_GEV["mu"])
        sigma_had_full = sum(
            ew.sigma_total(sqrt_s, q, m_f=ew.FERMION_MASS_GEV[q])
            for q in ["u", "d", "s", "c", "b"]
            if 2 * ew.FERMION_MASS_GEV[q] < sqrt_s
        )
        sigma_had_qed = sum(
            ew.sigma_qed_only(sqrt_s, q, m_f=ew.FERMION_MASS_GEV[q])
            for q in ["u", "d", "s", "c", "b"]
            if 2 * ew.FERMION_MASS_GEV[q] < sqrt_s
        )
        R_full = sigma_had_full / sigma_mumu_full
        R_qed = sigma_had_qed / sigma_mumu_qed
        print(f"  {sqrt_s:<10.1f}  {R_full:<12.4f}  {R_qed:<12.4f}  "
              f"{sigma_had_full/sigma_mumu_full:.4f} (Z), "
              f"{sigma_had_qed/sigma_mumu_qed:.4f} (γ-only)")

    print()
    print("  Note: at the Z pole, R ≈ 20 (not 11/3) because Γ(Z → had)/Γ(Z → ℓℓ)")
    print("  involves the (g_V² + g_A²) coupling structure rather than just")
    print("  electric charges Q_f².  The Z prefers to decay hadronically.")

    # =======================================================================
    # Cross-shim composition
    # =======================================================================

    header("CROSS-SHIM COMPOSITION TABLEAU")
    print()
    print(f"  {'Ingredient':<40s}  {'Source':<28s}  {'Notes':s}")
    print(f"  {'-'*40}  {'-'*28}  {'-'*25}")
    print(f"  {'Z mass M_Z, width Γ_Z':<40s}  "
          f"{'nwt.electroweak.constants':<28s}  {'PDG':s}")
    print(f"  {'g_V, g_A (chiral V-A couplings)':<40s}  "
          f"{'nwt.electroweak.couplings':<28s}  {'T_3 - 2 Q sin²θ_W':s}")
    print(f"  {'Z propagator (Breit-Wigner)':<40s}  "
          f"{'nwt.electroweak.process.chi':<28s}  {'1/(s - M_Z² + i M_Z Γ_Z)':s}")
    print(f"  {'γ + Z + interference σ_total':<40s}  "
          f"{'nwt.electroweak.sigma_total':<28s}  {'Halzen-Martin Eq 13.55-59':s}")
    print(f"  {'Cl(1,3) Dirac primitive':<40s}  "
          f"{'algebra.dirac':<28s}  {'shared with QED, QCD, QFT':s}")
    print(f"  {'Z partial widths Γ(Z → ff̄)':<40s}  "
          f"{'nwt.electroweak.decays':<28s}  {'(g_V²+g_A²) N_c β':s}")
    print()
    print("  The chiral V-A coupling structure is implemented at the Lagrangian")
    print("  vertex level, separate from the substrate's vector internal SU(2).")
    print("  Substrate primitives (Cl(1,3) spinors, Z propagator, Lorentz tensor")
    print("  amplitude) are shared with the QED, QCD, QFT, and string shims.")

    header("CONCLUSION", char="=")
    print(f"Z resonance peak σ(e⁺e⁻ → μ⁺μ⁻) at √s = M_Z:")
    print(f"  Substrate prediction  = {sigma_peak_substrate:.0f} pb")
    print(f"  PDG observed          = ~2000 pb (~2 nb)")
    print(f"  Agreement: {sigma_peak_substrate/2000.0*100:.1f}%")
    print()
    print(f"Total Z width: substrate {ew.total_width_Z():.4f} GeV vs PDG 2.4952 GeV")
    print(f"  Agreement: {ew.total_width_Z()/2.4952*100:.2f}%  (LO; residual from")
    print(f"  radiative corrections + QCD corrections to hadronic widths)")
    print()
    print("With the electroweak shim added, the cross-discipline demo arc")
    print("(Phases Q.8 R-ratio, Q.9 partonic DY, Q.10 hadronic DY, Q.11 Z")
    print("resonance) covers the full LHC W/Z-physics workhorse at LO,")
    print("using primitives shared across 6 view-shims.")
    print()


if __name__ == "__main__":
    main()
