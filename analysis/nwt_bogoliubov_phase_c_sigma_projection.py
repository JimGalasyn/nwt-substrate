"""Bogoliubov spectrum Phase C — σ-orbit projection of the abelian-Higgs
Bogoliubov spectrum onto K_7 Heffter σ-orbits.

Phase B established m_σ = √(2πα)·m_e ≈ 109 keV (Higgs) and m_A = √(4πα)·m_e
≈ 155 keV (gauge), with correlation lengths ξ_σ ≈ 4.67 λ̄_C and
ξ_A ≈ 3.30 λ̄_C — both LONGER than the framework's claimed
ξ_substrate = λ̄_C.

Phase C tests whether the K_7 σ-orbit projection produces a substrate
healing length λ̄_C as a projection invariant.  Strategy:

  1. Enumerate the 7 σ-orbits of K_7 Heffter (3 edges each).
  2. Define candidate projection schemes (polar-weighted, BPS-inverse,
     etc.).
  3. For each scheme + each σ-orbit, compute effective mass and
     correlation length.
  4. Report which scheme + orbit reproduces m_e (giving ξ = λ̄_C).

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_c_sigma_projection.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.abelian_higgs import (
    AbelianHiggsParams, M_E_EV, LAMBDA_C_M, LAMBDA_C_FM,
)
from nwt_substrate.condensate.sigma_orbits import (
    SIGMA_ORBITS, verify_partition, orbit_invariants,
    effective_mass, effective_correlation_length,
    wilson_product_factor,
)
from nwt_substrate.isa.constants import ALPHA_NWT

OUT_DIR = Path(__file__).parent / "phase_c_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE C — K_7 σ-orbit projection of Bogoliubov spectrum")
    print("=" * 78)
    print()

    # ---- 1. Verify σ-orbit partition ------------------------------------
    chk = verify_partition()
    print("σ-orbit partition of K_7 (21 edges → 7 orbits × 3 edges):")
    print(f"  n_orbits = {chk['n_orbits']}, n_edges = {chk['n_edges']}, "
          f"covers K_7: {chk['covers_K7']}")
    assert chk['covers_K7'], (
        f"σ-orbit partition fails! missing={chk['missing']}, extra={chk['extra']}")
    print()

    # ---- 2. Enumerate orbits with invariants ----------------------------
    p = AbelianHiggsParams.substrate_natural()
    print(f"Substrate-natural abelian-Higgs spectrum (Phase B):")
    print(f"  m_σ = √(2πα) · m_e = {p.m_sigma_eV/1e3:.3f} keV  = "
          f"{p.m_sigma_eV/M_E_EV:.4f} m_e")
    print(f"  m_A = √(4πα) · m_e = {p.m_A_eV/1e3:.3f} keV  = "
          f"{p.m_A_eV/M_E_EV:.4f} m_e")
    print(f"  ξ_σ = {p.xi_sigma_m/LAMBDA_C_M:.4f} λ̄_C")
    print(f"  ξ_A = {p.xi_A_m/LAMBDA_C_M:.4f} λ̄_C")
    print(f"  TARGET: ξ_substrate = λ̄_C (framework claim)")
    print()

    print("K_7 Heffter σ-orbits:")
    print(f"  {'id':<3} {'edges':<37} {'polar':>5} {'cross':>5} {'role'}")
    for oid in range(7):
        o = SIGMA_ORBITS[oid]
        inv = orbit_invariants(oid)
        edges_str = ', '.join(f"{a}-{b}" for (a, b) in o["edges"])
        print(f"  {oid:<3} {edges_str:<37} {inv.polar_edges:>5} "
              f"{inv.cross_edges:>5} {o['physical_role']}")
    print()
    print(f"  Wilson product per σ-orbit (all 7): (1 - √α)³ = "
          f"{wilson_product_factor(0):.5f} = f_J  (Layer 3 Phase F)")
    print()

    # ---- 3. Test candidate projection schemes ---------------------------
    schemes = ["higgs", "gauge", "polar_weighted", "topological",
               "BPS_inverse", "wilson_normalized"]

    print("=" * 78)
    print("CANDIDATE PROJECTION SCHEMES — effective mass per σ-orbit")
    print("=" * 78)
    print()
    print(f"{'orbit':<6} {'role':<32} " + "  ".join(f"{s:>11}" for s in schemes))
    print(f"{'    ':<6} {'    ':<32} " + "  ".join(f"{'m/m_e':>11}" for _ in schemes))
    print("-" * 130)
    for oid in range(7):
        o = SIGMA_ORBITS[oid]
        role = o['physical_role'][:30]
        masses = [effective_mass(oid, p, s) / M_E_EV for s in schemes]
        cells = "  ".join(f"{m:>11.4f}" for m in masses)
        print(f"{oid:<6} {role:<32} {cells}")
    print()

    # ---- 4. Test for ξ = λ̄_C match ---------------------------------------
    print("=" * 78)
    print("CORRELATION LENGTHS — does any scheme + orbit give ξ = λ̄_C?")
    print("=" * 78)
    print()
    print(f"{'orbit':<6} {'role':<32} " + "  ".join(f"{s:>11}" for s in schemes))
    print(f"{'    ':<6} {'    ':<32} " + "  ".join(f"{'ξ/λ̄_C':>11}" for _ in schemes))
    print("-" * 130)
    best_match = None
    best_diff = float('inf')
    for oid in range(7):
        o = SIGMA_ORBITS[oid]
        role = o['physical_role'][:30]
        xis = []
        for s in schemes:
            xi = effective_correlation_length(oid, p, s)
            r = xi / LAMBDA_C_M
            xis.append(r)
            d = abs(r - 1.0)
            if d < best_diff:
                best_diff = d
                best_match = (oid, s, r, xi)
        cells = "  ".join(f"{x:>11.4f}" for x in xis)
        print(f"{oid:<6} {role:<32} {cells}")
    print()
    print(f"  BEST match to ξ = λ̄_C:  orbit {best_match[0]} ({SIGMA_ORBITS[best_match[0]]['name']})")
    print(f"    scheme = {best_match[1]}")
    print(f"    ξ      = {best_match[2]:.4f} λ̄_C")
    print(f"    deviation: {(best_match[2] - 1.0)*100:+.2f}%")
    print()

    # ---- 5. Substrate-natural mass scale per orbit (with Wilson factor) ----
    print("=" * 78)
    print("PHYSICALLY MOTIVATED SCHEMES (with Wilson factor f_J)")
    print("=" * 78)
    print()
    # Idea: the σ-orbit excitation has effective mass set by the Wilson product
    # of edges times some base mass.  Test: m_eff = m_σ × Wilson, m_eff = m_e × Wilson, etc.
    f_J = wilson_product_factor(0)
    print(f"  f_J = (1-√α)³ = {f_J:.6f}")
    print()

    candidates_named = {
        "m_σ × f_J": p.m_sigma_eV * f_J,
        "m_A × f_J": p.m_A_eV * f_J,
        "m_σ / f_J": p.m_sigma_eV / f_J,
        "m_A / f_J": p.m_A_eV / f_J,
        "(m_σ + m_A) × f_J": (p.m_sigma_eV + p.m_A_eV) * f_J,
        "(m_σ + m_A) / 2": (p.m_sigma_eV + p.m_A_eV) / 2,
        "√(m_σ × m_A)": math.sqrt(p.m_sigma_eV * p.m_A_eV),
        "m_σ × √2 (degen. lift)": p.m_sigma_eV * math.sqrt(2),
        "m_A × √2": p.m_A_eV * math.sqrt(2),
        "m_σ / √α": p.m_sigma_eV / math.sqrt(ALPHA_NWT),
        "m_A / √α": p.m_A_eV / math.sqrt(ALPHA_NWT),
        "m_σ × π² (Derrick κ)": p.m_sigma_eV * math.pi ** 2,
        "m_σ × 1/√(2πα) (cancel α factor)": p.m_sigma_eV / math.sqrt(2 * math.pi * ALPHA_NWT),
        "m_A × 1/√(4πα) (cancel α factor)": p.m_A_eV / math.sqrt(4 * math.pi * ALPHA_NWT),
    }

    print(f"  {'combination':<45} {'mass (keV)':>11}  {'m/m_e':>9}  {'ξ/λ̄_C':>9}")
    print("  " + "-" * 90)
    target_m_e = M_E_EV
    for name, m in candidates_named.items():
        ratio = m / target_m_e
        xi_ratio = (target_m_e / m)  # = m_e/m, which equals ξ/λ̄_C
        marker = " ← matches m_e!" if abs(ratio - 1.0) < 0.05 else ""
        print(f"  {name:<45} {m/1e3:>11.3f}  {ratio:>9.4f}  "
              f"{xi_ratio:>9.4f}{marker}")
    print()

    # ---- 6. Plot ---------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # (a) σ-orbit effective masses
    ax = axes[0, 0]
    orbit_ids = list(range(7))
    schemes_show = ["higgs", "gauge", "polar_weighted", "topological",
                     "wilson_normalized"]
    x = np.arange(len(orbit_ids))
    w = 0.15
    for i, s in enumerate(schemes_show):
        masses = [effective_mass(oid, p, s) / M_E_EV for oid in orbit_ids]
        ax.bar(x + (i - len(schemes_show) / 2) * w, masses, width=w,
               label=s)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.7,
               label='target = m_e')
    ax.set_xticks(x)
    ax.set_xticklabels([f"σ_{i}" for i in orbit_ids], fontsize=9)
    ax.set_ylabel('m_eff / m_e')
    ax.set_title('σ-orbit effective masses by projection scheme')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)
    ax.set_yscale('log')

    # (b) σ-orbit correlation lengths (the ξ/λ̄_C plot)
    ax = axes[0, 1]
    for i, s in enumerate(schemes_show):
        xis = [effective_correlation_length(oid, p, s) / LAMBDA_C_M
                for oid in orbit_ids]
        ax.bar(x + (i - len(schemes_show) / 2) * w, xis, width=w,
               label=s)
    ax.axhline(1.0, color='red', linestyle='--', alpha=0.7,
               label='target = λ̄_C')
    ax.set_xticks(x)
    ax.set_xticklabels([f"σ_{i}" for i in orbit_ids], fontsize=9)
    ax.set_ylabel('ξ / λ̄_C')
    ax.set_title('σ-orbit correlation lengths by projection scheme')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)
    ax.set_yscale('log')

    # (c) σ-orbit topology — polar vs cross edge counts
    ax = axes[1, 0]
    polar_counts = [orbit_invariants(oid).polar_edges for oid in orbit_ids]
    cross_counts = [orbit_invariants(oid).cross_edges for oid in orbit_ids]
    ax.bar(x - 0.2, polar_counts, width=0.4, label='polar edges', color='C0')
    ax.bar(x + 0.2, cross_counts, width=0.4, label='cross edges', color='C3')
    ax.set_xticks(x)
    ax.set_xticklabels([f"σ_{i}\n({SIGMA_ORBITS[i]['name'][:8]})"
                         for i in orbit_ids],
                        fontsize=7, rotation=20, ha='right')
    ax.set_ylabel('# edges')
    ax.set_title('σ-orbit topology (polar P-incident vs E↔F cross)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis='y')

    # (d) candidates_named comparison panel
    ax = axes[1, 1]
    names = list(candidates_named.keys())
    ratios = [m / target_m_e for m in candidates_named.values()]
    ys = np.arange(len(names))
    colors = ['C2' if abs(r - 1.0) < 0.05 else 'C0' for r in ratios]
    ax.barh(ys, ratios, color=colors)
    ax.axvline(1.0, color='red', linestyle='--', alpha=0.7,
                label='target = 1 (m_e)')
    ax.set_yticks(ys)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel('m_combination / m_e')
    ax.set_title('Heuristic combinations of (m_σ, m_A) tested against m_e')
    ax.legend(fontsize=9)
    ax.set_xscale('log')
    ax.grid(alpha=0.3, axis='x')

    fig.suptitle(
        f"Bogoliubov Phase C — K_7 σ-orbit projection test\n"
        f"Does any scheme yield ξ = λ̄_C as a projection invariant?",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_c_sigma_projection.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    # ---- 7. Save numerical outputs --------------------------------------
    schemes_full = ["higgs", "gauge", "polar_weighted", "topological",
                     "BPS_inverse", "wilson_normalized"]
    masses_array = np.array([
        [effective_mass(oid, p, s) / M_E_EV for s in schemes_full]
        for oid in range(7)
    ])
    np.savez(OUT_DIR / "phase_c_sigma_projection.npz",
             orbit_ids=np.array(list(range(7))),
             schemes=np.array(schemes_full),
             masses_over_m_e=masses_array,
             alpha_nwt=ALPHA_NWT,
             m_sigma_over_m_e=p.m_sigma_eV / M_E_EV,
             m_A_over_m_e=p.m_A_eV / M_E_EV,
             f_J=wilson_product_factor(0),
             candidate_names=np.array(list(candidates_named.keys())),
             candidate_masses_over_m_e=np.array(
                [m / target_m_e for m in candidates_named.values()]))
    print(f"  data saved {OUT_DIR / 'phase_c_sigma_projection.npz'}")

    # ---- 8. Headline ----------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase C σ-orbit projection")
    print("=" * 78)
    print()
    print(f"  ★ Best scheme + orbit matching ξ = λ̄_C:")
    print(f"    σ-orbit {best_match[0]} ({SIGMA_ORBITS[best_match[0]]['name']})")
    print(f"    scheme = '{best_match[1]}'")
    print(f"    ξ = {best_match[2]:.4f} × λ̄_C  "
          f"({(best_match[2]-1)*100:+.2f}% deviation)")
    print()
    matching = [(name, m / target_m_e)
                for name, m in candidates_named.items()
                if abs(m / target_m_e - 1.0) < 0.05]
    if matching:
        print(f"  Heuristic combinations matching m_e within 5%:")
        for name, r in matching:
            print(f"    {name:<45} ratio = {r:.4f}")
    else:
        print(f"  No simple heuristic combination matches m_e within 5%.")
        print(f"  The ξ = λ̄_C identification likely requires first-principles")
        print(f"  derivation from σ-orbit Wilson-loop dynamics (Phase D).")


if __name__ == "__main__":
    main()
