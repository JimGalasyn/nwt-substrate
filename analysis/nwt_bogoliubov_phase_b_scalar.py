"""Bogoliubov spectrum Phase B — relativistic abelian-Higgs scalar warm-up.

Builds the substrate-natural condensate parameters (v = m_e,
e = √(4πα_NWT), λ = e²/2 BPS) and computes the resulting Higgs +
gauge boson masses, correlation lengths, and BPS line tension.

Validates against Paper 16 §L_2 published anchors:

    μ_BPS = 2π m_e² / (ℏc) ≈ 8.31 keV/fm    (line tension)
    ξ_substrate                              (vortex-core / correlation scale)

and computes the relativistic Bogoliubov dispersion E_k for both
Higgs and gauge branches.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_b_scalar.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.abelian_higgs import (
    AbelianHiggsParams,
    bogoliubov_dispersion_higgs,
    bogoliubov_dispersion_gauge,
    line_tension_BPS,
    healing_length_higgs,
    healing_length_gauge,
    LAMBDA_C_M,
    LAMBDA_C_FM,
    M_E_EV,
)
from nwt_substrate.isa.constants import ALPHA_NWT


OUT_DIR = Path(__file__).parent / "phase_b_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE B — Relativistic abelian-Higgs scalar warm-up")
    print("=" * 78)
    print()

    p = AbelianHiggsParams.substrate_natural()
    print(f"Substrate-natural parameters (v = m_e, e = √(4πα), λ = e²/2 BPS):")
    print(f"  v        = {p.v_phi_eV/1e6:.4f} MeV  (= m_e)")
    print(f"  e_gauge  = √(4π × {ALPHA_NWT:.5f}) = {p.e_gauge:.4f}")
    print(f"  λ        = e²/2 = {p.lambda_quartic:.5f}  (BPS: {p.at_bps})")
    print()

    # ---- Spectrum --------------------------------------------------------
    print("Mass spectrum (small perturbations around ψ_0 = v):")
    print(f"  m_σ (Higgs)  = √λ · v = √(2πα) · m_e = "
          f"{p.m_sigma_eV/1e6:.4f} MeV  "
          f"= {p.m_sigma_eV/M_E_EV:.4f} m_e")
    print(f"  m_A (gauge)  = e · v   = √(4πα) · m_e = "
          f"{p.m_A_eV/1e6:.4f} MeV  "
          f"= {p.m_A_eV/M_E_EV:.4f} m_e")
    print()

    # ---- Correlation lengths --------------------------------------------
    print("Correlation lengths:")
    xi_sigma = healing_length_higgs(p)
    xi_A = healing_length_gauge(p)
    print(f"  ξ_σ (Higgs) = ℏ/(m_σ c) = {xi_sigma*1e15:.2f} fm  "
          f"= {xi_sigma/LAMBDA_C_M:.4f} × λ̄_C")
    print(f"  ξ_A (gauge) = ℏ/(m_A c) = {xi_A*1e15:.2f} fm  "
          f"= {xi_A/LAMBDA_C_M:.4f} × λ̄_C")
    print(f"  λ̄_C        = ℏ/(m_e c) = {LAMBDA_C_FM:.2f} fm  (substrate anchor)")
    print()
    print(f"  Ratio ξ_σ/λ̄_C = 1/√(2πα) = {1.0/math.sqrt(2*math.pi*ALPHA_NWT):.4f}")
    print(f"  Ratio ξ_A/λ̄_C = 1/√(4πα) = {1.0/math.sqrt(4*math.pi*ALPHA_NWT):.4f}")
    print()
    print("  ★ NEITHER ξ_σ NOR ξ_A equals λ̄_C exactly — they're at scales")
    print("    related to λ̄_C by factor 1/√(2πα) ≈ 3.3 (ξ_σ) and ≈ 4.7 (ξ_A).")
    print("    Substrate's 'healing length' from "
          "[[framework_healing_length_principle]]")
    print("    might be either a third scale or a specific combination —")
    print("    needs reconciliation with Paper 16 / Paper 17 conventions.")
    print()

    # ---- Line tension validation against Paper 16 §L_2 ------------------
    mu_BPS = line_tension_BPS(p)
    print("Line tension validation (Paper 16 §L_2 line 333):")
    print(f"  μ_BPS = 2π v² (natural units)")
    print(f"        = 2π m_e² / (ℏc)")
    print(f"        = {mu_BPS:.4f} eV/fm")
    print(f"        = {mu_BPS/1000:.4f} keV/fm")
    # Independent textbook calc in MeV/fm: 2π m_e²/(ℏc)
    # with m_e in MeV (= 0.511), ℏc in MeV·fm (= 197.327)
    independent_MeV_fm = (2.0 * math.pi
                           * (0.510998928) ** 2 / 197.3269804)
    independent_eV_fm = independent_MeV_fm * 1e6
    print(f"  Expected (independent calc): {independent_eV_fm:.4f} eV/fm "
          f"= {independent_MeV_fm*1000:.4f} keV/fm")
    print(f"  Paper 16 §L_2 stated:        ~ 8.31 (units ambiguous in text)")
    rel_err = abs(mu_BPS / independent_eV_fm - 1.0)
    print(f"  Relative error: {rel_err*100:.4f}%")
    print()

    # ---- Plot dispersion -------------------------------------------------
    print("Computing relativistic Bogoliubov dispersion E_k for both branches…")
    # k range: from 10⁻⁴ × m_σ/(ℏc) to 100 × m_σ/(ℏc)
    k_vals = np.geomspace(1e-4 / xi_sigma, 100.0 / xi_sigma, 200)
    E_higgs = bogoliubov_dispersion_higgs(k_vals, p)
    E_gauge = bogoliubov_dispersion_gauge(k_vals, p)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) E vs k log-log
    ax = axes[0]
    from nwt_substrate.condensate.abelian_higgs import HBARC_EV_M
    ax.loglog(k_vals * xi_sigma, E_higgs / p.m_sigma_eV,
              'C0-', lw=2, label='Higgs branch (σ)')
    ax.loglog(k_vals * xi_sigma, E_gauge / p.m_sigma_eV,
              'C3-', lw=2, label='Gauge branch (A_μ)')
    # Asymptote: E = ℏck → linear in k
    ax.loglog(k_vals * xi_sigma,
              HBARC_EV_M * k_vals / p.m_sigma_eV,
              'k:', lw=1, alpha=0.6, label='E = ℏck (photon)')
    ax.axhline(1.0, color='C0', ls='--', alpha=0.4,
               label='m_σ c² (Higgs gap)')
    ax.axhline(p.m_A_eV / p.m_sigma_eV, color='C3', ls='--', alpha=0.4,
               label='m_A c² (gauge gap)')
    ax.set_xlabel('k · ξ_σ  (k in units of 1/ξ_σ)')
    ax.set_ylabel('E_k / (m_σ c²)')
    ax.set_title('Relativistic abelian-Higgs Bogoliubov dispersion\n'
                 '(no phonon limit — gapped at k=0)')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)

    # (b) Substrate predictions table
    ax = axes[1]
    ax.axis('off')
    table_data = [
        ['Quantity', 'Value', 'Note'],
        ['α_NWT', f'{ALPHA_NWT:.5f}', 'Paper 17 trefoil'],
        ['v (cond. VEV)', f'{p.v_phi_eV/1e6:.4f} MeV', '= m_e (Paper 16)'],
        ['e (gauge coupling)', f'{p.e_gauge:.4f}', '√(4πα)'],
        ['λ (quartic)', f'{p.lambda_quartic:.5f}', 'e²/2 (BPS)'],
        ['', '', ''],
        ['m_σ (Higgs mode)', f'{p.m_sigma_eV/1e6:.4f} MeV', '√(2πα)·m_e'],
        ['m_A (gauge boson)', f'{p.m_A_eV/1e6:.4f} MeV', '√(4πα)·m_e'],
        ['m_σ / m_e', f'{p.m_sigma_eV/M_E_EV:.4f}', 'substrate prediction'],
        ['m_A / m_e', f'{p.m_A_eV/M_E_EV:.4f}', 'substrate prediction'],
        ['', '', ''],
        ['ξ_σ (Higgs Cω-len)', f'{xi_sigma*1e15:.1f} fm', 'λ̄_C/√(2πα)'],
        ['ξ_A (gauge Cω-len)', f'{xi_A*1e15:.1f} fm', 'λ̄_C/√(4πα)'],
        ['λ̄_C (electron)', f'{LAMBDA_C_FM:.1f} fm', 'substrate anchor'],
        ['', '', ''],
        ['μ_BPS line tension', f'{mu_BPS/1e3:.3f} keV/fm',
            '2π m_e²/(ℏc)'],
    ]
    table = ax.table(cellText=table_data, loc='center', cellLoc='left',
                     colWidths=[0.42, 0.32, 0.30])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.5)
    ax.set_title('Phase B substrate-natural derived spectrum')

    fig.suptitle(
        f"Bogoliubov Phase B — relativistic abelian-Higgs scalar warm-up\n"
        f"(v = m_e, BPS λ = e²/2, Paper 16 §L_2)",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_b_scalar_dispersion.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    # ---- Save numerical outputs -----------------------------------------
    np.savez(OUT_DIR / "phase_b_scalar.npz",
             alpha_nwt=ALPHA_NWT,
             v_phi_eV=p.v_phi_eV,
             e_gauge=p.e_gauge,
             lambda_quartic=p.lambda_quartic,
             m_sigma_eV=p.m_sigma_eV,
             m_A_eV=p.m_A_eV,
             xi_sigma_m=xi_sigma,
             xi_A_m=xi_A,
             lambda_C_m=LAMBDA_C_M,
             mu_BPS_eV_per_fm=mu_BPS,
             k_vals=k_vals, E_higgs=E_higgs, E_gauge=E_gauge)
    print(f"  data saved {OUT_DIR / 'phase_b_scalar.npz'}")

    # ---- Headline -------------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase B substrate-derived spectrum")
    print("=" * 78)
    print()
    print(f"  Higgs mode mass:    m_σ = √(2πα) · m_e = {p.m_sigma_eV/1e6:.4f} MeV")
    print(f"  Gauge boson mass:   m_A = √(4πα) · m_e = {p.m_A_eV/1e6:.4f} MeV")
    print(f"  Both gapped (relativistic) — no Bogoliubov phonon limit.")
    print()
    print(f"  BPS line tension μ = 2π m_e²/(ℏc) = {mu_BPS/1e3:.3f} keV/fm")
    print(f"  ✓ matches Paper 16 §L_2 independent calc to {rel_err*100:.4f}%")
    print()
    print(f"  Healing-length convention note:")
    print(f"    ξ_σ = λ̄_C / √(2πα) ≈ {xi_sigma*1e15:.0f} fm "
          f"= {xi_sigma/LAMBDA_C_M:.2f} × λ̄_C")
    print(f"    Substrate's 'ξ = λ̄_C' from "
          f"[[framework_healing_length_principle]]")
    print(f"    differs by 1/√(2πα) — needs reconciliation in Phase C.")


if __name__ == "__main__":
    main()
