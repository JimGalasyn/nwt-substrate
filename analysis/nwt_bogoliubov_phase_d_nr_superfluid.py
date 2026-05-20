"""Bogoliubov spectrum Phase D — NR superfluid + Lorentz-invariance
constraint reproduces ξ_substrate = λ̄_C.

Phase C's negative result (no linear σ-orbit projection of the
relativistic Bogoliubov spectrum gives ξ = λ̄_C) is RESOLVED by Paper 5
§II: ξ = λ̄_C is enforced by the Lorentz-invariance condition c_s = c
in the NON-RELATIVISTIC LIMIT of the abelian-Higgs condensate.

Two healing-length scales coexist:
  ξ_substrate = λ̄_C        NR superfluid (this module)
  ξ_σ ≈ 4.67 λ̄_C           Higgs Compton (Phase B relativistic gap)
  ξ_A ≈ 3.30 λ̄_C           Gauge Compton (Phase B Higgs-mechanism)

These are different observables.  The framework's ξ_substrate is the
NR healing length, enforced by relativistic consistency (c_s = c).

This script:
  1. Verifies that with m* = m_e/√2 and g n_0 = m_e c²/√2:
       ξ = λ̄_C to machine precision  AND  c_s = c to machine precision
  2. Compares to the Phase B/C relativistic gap-mode scales
  3. Lays groundwork for Phase D-2: σ-orbit Wilson-loop dynamics in
     the NR framework (vortex-ring energies for torus-knot bound
     states, Paper 5 / Paper 11-13).

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_d_nr_superfluid.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.abelian_higgs import (
    AbelianHiggsParams, M_E_EV, LAMBDA_C_M, LAMBDA_C_FM,
)
from nwt_substrate.condensate.nr_superfluid import (
    substrate_natural_NR, verify_substrate_natural, NRSuperfluidParams,
)
from nwt_substrate.gravity.constants import (
    C_LIGHT_M_S as C_LIGHT_SI, M_ELECTRON_KG as M_E_KG, HBAR_J_S as HBAR_SI,
)
from nwt_substrate.isa.constants import ALPHA_NWT


OUT_DIR = Path(__file__).parent / "phase_d_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE D — NR superfluid limit + Lorentz constraint")
    print("=" * 78)
    print()

    # ---- 1. Verify substrate-natural NR parameters ----------------------
    chk = verify_substrate_natural()
    p = chk["params"]
    print("Paper 5 §II substrate-natural NR superfluid (verification):")
    print(f"  m*  (effective boson mass) = m_e/√2 = "
          f"{p.m_star_kg*1e30:.4e} × 10⁻³⁰ kg "
          f"= {M_E_KG/(math.sqrt(2)*1.602176634e-19/9e16):.4f} ... hmm")
    print(f"  m*  = {p.m_star_kg} kg = "
          f"{p.m_star_kg/M_E_KG:.4f} m_e")
    print(f"  ρ_0 = {p.rho_0_kg_m3:.3e} kg/m³  (Paper 5 calibration)")
    print(f"  n_0 = {p.n_0_per_m3:.3e} 1/m³")
    print(f"  g   = {p.g_SI:.3e} J·m³")
    print(f"  μ_chem = g n_0 = {chk['mu_chem_eV']/1e3:.3f} keV  "
          f"(= m_e c²/√2 = {M_E_EV/math.sqrt(2)/1e3:.3f} keV)")
    print()
    print(f"Lorentz-invariance consistency checks:")
    print(f"  ξ_NR (healing length)   = {chk['xi_m']*1e15:.4f} fm")
    print(f"  λ̄_C (electron Compton)  = {LAMBDA_C_M*1e15:.4f} fm")
    print(f"  ξ / λ̄_C = {chk['xi_over_lambda_C']:.10f}  "
          f"{'★ EXACT' if chk['xi_match'] else '✗ MISMATCH'}")
    print()
    print(f"  c_s (sound speed) = {chk['c_s_m_s']:.3e} m/s")
    print(f"  c (speed of light) = {C_LIGHT_SI:.3e} m/s")
    print(f"  c_s / c = {chk['c_s_over_c']:.10f}  "
          f"{'★ EXACT' if chk['c_s_match'] else '✗ MISMATCH'}")
    print()
    print(f"  ★ Paper 5 §II: ξ = λ̄_C ⇔ c_s = c.  Both consequences of the")
    print(f"    same Lorentz-invariance constraint g n_0 = m* c².")
    print()

    # ---- 2. Compare to Phase B relativistic gap scales ------------------
    ah = AbelianHiggsParams.substrate_natural()
    print("=" * 78)
    print("RECONCILIATION — three length scales in the condensate")
    print("=" * 78)
    print()
    print(f"  {'scale':<32} {'value (fm)':>12} {'/ λ̄_C':>9} {'source'}")
    print("  " + "-" * 76)
    print(f"  {'ξ_substrate (NR healing)':<32} "
          f"{chk['xi_m']*1e15:>12.4f} {chk['xi_over_lambda_C']:>9.4f}  "
          f"Paper 5 §II")
    print(f"  {'λ̄_C (e⁻ Compton)':<32} "
          f"{LAMBDA_C_FM:>12.4f} {1.0:>9.4f}  "
          f"definition")
    print(f"  {'ξ_σ (Higgs Compton)':<32} "
          f"{ah.xi_sigma_m*1e15:>12.4f} "
          f"{ah.xi_sigma_m/LAMBDA_C_M:>9.4f}  Phase B rel.")
    print(f"  {'ξ_A (gauge Compton)':<32} "
          f"{ah.xi_A_m*1e15:>12.4f} "
          f"{ah.xi_A_m/LAMBDA_C_M:>9.4f}  Phase B rel.")
    print()
    print(f"  ★ ξ_substrate ≠ ξ_σ ≠ ξ_A — these are DISTINCT observables:")
    print(f"      ξ_substrate is the NR healing length of the bulk condensate")
    print(f"      ξ_σ, ξ_A are Compton wavelengths of relativistic gap modes")
    print(f"    Phase C's negative result (no projection of m_σ/m_A → λ̄_C)")
    print(f"    was correct — they are different quantities by construction.")
    print()

    # ---- 3. Energy-scale ladder for reference ---------------------------
    print("=" * 78)
    print("MASS / ENERGY LADDER (substrate-natural)")
    print("=" * 78)
    print()
    print(f"  {'quantity':<35} {'value':>16}  {'note'}")
    print("  " + "-" * 80)
    print(f"  {'α_NWT (Paper 17)':<35} "
          f"{ALPHA_NWT:>16.6e}  '1/(25π√3 + 1)'")
    print(f"  {'m_e (electron rest energy)':<35} "
          f"{M_E_EV/1e3:>13.3f} keV  Paper 16 substrate VEV")
    print(f"  {'m*  (effective NR boson)':<35} "
          f"{M_E_EV/math.sqrt(2)/1e3:>13.3f} keV  m_e/√2")
    print(f"  {'μ_chem (chemical potential)':<35} "
          f"{chk['mu_chem_eV']/1e3:>13.3f} keV  g n_0 = m_e c²/√2")
    print(f"  {'m_σ (relativistic Higgs gap)':<35} "
          f"{ah.m_sigma_eV/1e3:>13.3f} keV  √(2πα) m_e (Phase B)")
    print(f"  {'m_A (relativistic gauge gap)':<35} "
          f"{ah.m_A_eV/1e3:>13.3f} keV  √(4πα) m_e (Phase B)")
    print(f"  {'μ_BPS (line tension)':<35} "
          f"{ah.mu_BPS_eV_per_fm/1e3:>13.3f} keV/fm  "
          f"2π m_e²/(ℏc), Phase B")
    print()

    # ---- 4. Plot ---------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) Length scales bar chart
    ax = axes[0]
    scales_fm = [
        ('ξ_substrate\n(NR healing)', chk['xi_m']*1e15),
        ('λ̄_C\n(e⁻ Compton)', LAMBDA_C_FM),
        ('ξ_σ\n(Higgs Compton)', ah.xi_sigma_m*1e15),
        ('ξ_A\n(gauge Compton)', ah.xi_A_m*1e15),
    ]
    names = [s[0] for s in scales_fm]
    vals = [s[1] for s in scales_fm]
    colors = ['C2', 'C3', 'C0', 'C1']
    bars = ax.bar(names, vals, color=colors)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, v * 1.05,
                f'{v:.0f} fm', ha='center', fontsize=9)
    ax.axhline(LAMBDA_C_FM, color='C3', ls='--', alpha=0.5,
               label=f'λ̄_C = {LAMBDA_C_FM:.0f} fm')
    ax.set_ylabel('length (fm)')
    ax.set_title('Three length scales in the substrate condensate\n'
                 '(Phase C resolved: distinct observables, not all = λ̄_C)')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    # (b) Mass-scale ladder
    ax = axes[1]
    masses_keV = [
        ('m_σ\n(Higgs gap)', ah.m_sigma_eV/1e3),
        ('m_A\n(gauge gap)', ah.m_A_eV/1e3),
        ('m*\n(NR effective)', M_E_EV/math.sqrt(2)/1e3),
        ('μ_chem\n(g n_0)', chk['mu_chem_eV']/1e3),
        ('m_e\n(electron)', M_E_EV/1e3),
    ]
    names = [s[0] for s in masses_keV]
    vals = [s[1] for s in masses_keV]
    bars = ax.bar(names, vals,
                   color=['C0', 'C1', 'C4', 'C5', 'C3'])
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, v * 1.02,
                f'{v:.0f} keV', ha='center', fontsize=9)
    ax.axhline(M_E_EV/1e3, color='C3', ls='--', alpha=0.5,
               label=f'm_e = {M_E_EV/1e3:.0f} keV')
    ax.set_ylabel('energy (keV)')
    ax.set_title('Mass ladder in substrate condensate')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')

    fig.suptitle(
        f"Phase D — NR superfluid limit reproduces ξ_substrate = λ̄_C\n"
        f"(Paper 5 §II: ξ = λ̄_C ⇔ c_s = c, Lorentz invariance)",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_d_nr_superfluid.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_d_nr_superfluid.npz",
             alpha_nwt=ALPHA_NWT,
             xi_substrate_m=chk['xi_m'],
             lambda_C_m=LAMBDA_C_M,
             xi_substrate_over_lambda_C=chk['xi_over_lambda_C'],
             c_s_over_c=chk['c_s_over_c'],
             m_sigma_eV=ah.m_sigma_eV,
             m_A_eV=ah.m_A_eV,
             m_star_kg=p.m_star_kg,
             mu_chem_J=p.mu_chem_J,
             rho_0_kg_m3=p.rho_0_kg_m3,
             n_0_per_m3=p.n_0_per_m3,
             g_SI=p.g_SI)
    print(f"  data saved {OUT_DIR / 'phase_d_nr_superfluid.npz'}")

    # ---- 5. Headline ----------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase D NR-superfluid reconciliation")
    print("=" * 78)
    print()
    print(f"  ★ ξ_substrate = λ̄_C reproduced EXACTLY in NR superfluid limit")
    print(f"    with substrate-natural choice m* = m_e/√2, g n_0 = m_e c²/√2")
    print(f"    (Paper 5 §II Lorentz-invariance constraint).")
    print(f"  ★ Equivalently, c_s = c (sound speed matches speed of light)")
    print()
    print(f"  ★ Phase C's negative result EXPLAINED: the relativistic Phase B")
    print(f"    gap masses m_σ, m_A produce Compton wavelengths ξ_σ, ξ_A ≠ λ̄_C")
    print(f"    by construction — they're different observables from the bulk")
    print(f"    healing length.  No projection between them is required.")
    print()
    print(f"  ★ NEXT (Phase D-2): σ-orbit Wilson-loop dynamics for bound-state")
    print(f"    masses (electron from (2,1) torus knot per Paper 5 §III,")
    print(f"    proton from (1,4) knot per Paper 13).  σ-orbit topology")
    print(f"    selects which knot types are stable in each substrate channel.")


if __name__ == "__main__":
    main()
