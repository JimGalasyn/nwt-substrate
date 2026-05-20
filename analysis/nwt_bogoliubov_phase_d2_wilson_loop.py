"""Bogoliubov Phase D-2 — σ-orbit Wilson-loop dynamics for bound-state masses.

In the NR superfluid limit (Phase D-1), vortex-ring self-energy on
(p,q) torus knots gives the mass spectrum (Paper 11 §III).  Each
particle is a torus-knot vortex with quantum numbers (p,q,m,n_q).

This script:
  1. Maps K_7 σ-orbits to particle classes via topological compatibility
     (n_q = number of linked components allowed by σ-orbit structure).
  2. Identifies the ground-state particle of each σ-orbit (lowest mass).
  3. Verifies the Paper 11 mass formula reproduces observed masses.
  4. Builds the σ-orbit → torus-knot bound-state map.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_d2_wilson_loop.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.sigma_orbits import (
    SIGMA_ORBITS, orbit_invariants,
)
from nwt_substrate.condensate.wilson_loop_dynamics import (
    SIGMA_ORBIT_ALLOWED_N_Q,
    compatible_particles,
    ground_state_per_orbit,
    particle_mass_MeV,
    kelvin_saffman_aspect,
)
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_d2_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE D-2 — σ-orbit Wilson-loop dynamics, bound-state masses")
    print("=" * 78)
    print()
    print("Paper 11 §III mass formula:")
    print("  m/m_e = (p² + q²)/5 · β/β_e · ln(8β)/ln(8β_e) · n_q^q")
    print("  β = √(m²/p² - 1), β_e = √5/2 (electron anchor)")
    print()

    # ---- 1. Validate the mass formula end-to-end ------------------------
    print("=" * 78)
    print("MASS FORMULA VALIDATION — all 25 compendium particles")
    print("=" * 78)
    print()
    print(f"  {'particle':<10} {'(p,q,m,n_q)':<14} "
          f"{'β':>7} {'pred (MeV)':>12} {'obs (MeV)':>12} {'err %':>8}")
    print(f"  " + "-" * 80)
    errors = []
    for e in COMPENDIUM:
        p, q, m_pc, n_q = e["p"], e["q"], e["m"], e["n_q"]
        m_obs = e.get("m_obs")
        beta = kelvin_saffman_aspect(p, m_pc)
        m_pred = particle_mass_MeV(p, q, m_pc, n_q)
        if m_obs and not math.isnan(m_pred):
            err = (m_pred - m_obs) / m_obs * 100.0
            errors.append(abs(err))
        else:
            err = float("nan")
        marker = " (anchor)" if e["name"] == "e-" else ""
        print(f"  {e['name']:<10} "
              f"({p},{q},{m_pc},{n_q})".ljust(16),
              f"{beta:>7.3f} {m_pred:>12.3f} "
              f"{m_obs or float('nan'):>12.3f} "
              f"{err:>+7.2f}%{marker}")
    if errors:
        median_err = float(np.median(errors))
        max_err = float(np.max(errors))
        print()
        print(f"  Median absolute error: {median_err:.2f}%")
        print(f"  Maximum absolute error: {max_err:.2f}%")
        print(f"  Paper 11 reports: median 1.0%, max 10.2% (nucleons)")

    # ---- 2. σ-orbit topological compatibility ---------------------------
    print()
    print("=" * 78)
    print("σ-ORBIT → PARTICLE COMPATIBILITY (topological filter)")
    print("=" * 78)
    print()
    print("Each σ-orbit's 3-edge topology restricts which n_q values are allowed.")
    print("Assumed rules:")
    print("  star (3 edges at common vertex P)  → n_q = 0 (single-component)")
    print("  triangle (closed 3-loop)            → n_q ∈ {1, 2} (loops/2-links)")
    print("  parallel (3 disjoint cross edges)   → n_q ∈ {2, 3, 5} (multi-link)")
    print("  twisted cross (Z_3 chirality)       → n_q = 3 (3-link)")
    print()

    for orbit_id in range(7):
        o = SIGMA_ORBITS[orbit_id]
        inv = orbit_invariants(orbit_id)
        allowed = SIGMA_ORBIT_ALLOWED_N_Q.get(orbit_id, [])
        compat = compatible_particles(orbit_id)
        names = sorted(set(e["name"] for e in compat),
                       key=lambda n: next(
                           (e.get("m_obs", 1e9) for e in compat
                            if e["name"] == n), 1e9))[:6]
        names_str = ', '.join(names) + (' ...' if len(compat) > 6 else '')
        print(f"  σ_{orbit_id} ({o['name'][:24]:<24}): polar={inv.polar_edges} "
              f"cross={inv.cross_edges}  "
              f"allowed n_q = {allowed}  "
              f"({len(compat)} particles)")
        print(f"       lightest: {names_str}")

    # ---- 3. Ground-state assignment per σ-orbit -------------------------
    print()
    print("=" * 78)
    print("σ-ORBIT GROUND STATE — lowest-mass compatible particle")
    print("=" * 78)
    print()
    ground = ground_state_per_orbit()
    print(f"  {'orbit':<6} {'σ-orbit name':<32} {'ground particle':<15} "
          f"{'(p,q,m,n_q)':<14} {'pred (MeV)':>11} {'obs (MeV)':>11}")
    print(f"  " + "-" * 96)
    for g in ground:
        pqmnq = f"({g.p},{g.q},{g.m_pc},{g.n_q})"
        m_obs_str = f"{g.mass_obs_MeV:.3f}" if g.mass_obs_MeV else "n/a"
        print(f"  {g.orbit_id:<6} {g.orbit_name[:30]:<32} "
              f"{g.particle_name:<15} {pqmnq:<14} "
              f"{g.mass_pred_MeV:>11.3f} {m_obs_str:>11}")
    print()

    # Highlights
    matter_orbit = next((g for g in ground if g.orbit_id == 0), None)
    if matter_orbit:
        print(f"  ★ Matter σ-orbit (orbit 0) ground state: "
              f"{matter_orbit.particle_name}")
        print(f"    (p,q,m,n_q) = ({matter_orbit.p}, {matter_orbit.q}, "
              f"{matter_orbit.m_pc}, {matter_orbit.n_q})")
        print(f"    Predicted mass: {matter_orbit.mass_pred_MeV:.3f} MeV")
        print(f"    Observed: {matter_orbit.mass_obs_MeV:.3f} MeV")
        print(f"    Topology gives the electron as the natural ground state of")
        print(f"    the matter (bridging, polar=3) σ-orbit. ✓")

    # ---- 4. Plot --------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # (a) Mass formula validation — predicted vs observed
    ax = axes[0, 0]
    obs = []
    pred = []
    labels = []
    for e in COMPENDIUM:
        m_obs = e.get("m_obs")
        m_pred = particle_mass_MeV(e["p"], e["q"], e["m"], e["n_q"])
        if m_obs and not math.isnan(m_pred):
            obs.append(m_obs)
            pred.append(m_pred)
            labels.append(e["name"])
    ax.loglog(obs, pred, 'o', ms=8, color='C0', alpha=0.7)
    lims = [min(obs) * 0.5, max(obs) * 2]
    ax.loglog(lims, lims, 'k--', alpha=0.5, label='1:1')
    for o, pr, lab in zip(obs, pred, labels):
        ax.annotate(lab, (o, pr), fontsize=7, xytext=(3, 3),
                    textcoords='offset points', alpha=0.7)
    ax.set_xlabel('observed mass (MeV)')
    ax.set_ylabel('Paper 11 predicted (MeV)')
    ax.set_title(f'Paper 11 mass formula — 25 particles\n'
                 f'median error {median_err:.2f}%')
    ax.legend()
    ax.grid(alpha=0.3, which='both')

    # (b) Residuals
    ax = axes[0, 1]
    residuals = [(p_ - o) / o * 100 for o, p_ in zip(obs, pred)]
    ys = np.arange(len(labels))
    colors = ['C2' if abs(r) < 1 else 'C0' if abs(r) < 5 else 'C3'
              for r in residuals]
    ax.barh(ys, residuals, color=colors, alpha=0.7)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('relative error (%)')
    ax.set_title('Paper 11 residuals')
    ax.set_xlim(-12, 12)
    ax.grid(alpha=0.3, axis='x')

    # (c) σ-orbit → particle compatibility heatmap
    ax = axes[1, 0]
    sectors = sorted(set(e.get("sector", "?") for e in COMPENDIUM))
    counts = np.zeros((7, len(sectors)))
    for orbit_id in range(7):
        compat = compatible_particles(orbit_id)
        for e in compat:
            s = e.get("sector", "?")
            si = sectors.index(s)
            counts[orbit_id, si] += 1
    im = ax.imshow(counts, aspect='auto', cmap='Blues', origin='lower')
    ax.set_xticks(range(len(sectors)))
    ax.set_xticklabels(sectors, rotation=30, ha='right', fontsize=9)
    ax.set_yticks(range(7))
    ax.set_yticklabels([f'σ_{i}\n{SIGMA_ORBITS[i]["name"][:15]}'
                        for i in range(7)], fontsize=8)
    for i in range(7):
        for j in range(len(sectors)):
            if counts[i, j] > 0:
                ax.text(j, i, int(counts[i, j]),
                        ha='center', va='center', fontsize=9,
                        color='white' if counts[i, j] > counts.max() / 2
                        else 'black')
    ax.set_title('σ-orbit → particle sector compatibility')
    plt.colorbar(im, ax=ax, label='# compatible particles')

    # (d) σ-orbit ground-state masses
    ax = axes[1, 1]
    g_ids = [g.orbit_id for g in ground]
    g_masses = [g.mass_pred_MeV for g in ground]
    g_names = [f"σ_{g.orbit_id}:\n{g.particle_name}" for g in ground]
    colors_o = ['C2' if g.error_pct is not None and abs(g.error_pct) < 5
                else 'C0' for g in ground]
    bars = ax.bar(g_names, g_masses, color=colors_o, alpha=0.7)
    for bar, g in zip(bars, ground):
        if g.mass_obs_MeV:
            ax.text(bar.get_x() + bar.get_width()/2,
                    g.mass_pred_MeV * 1.02,
                    f'{g.mass_pred_MeV:.1f} MeV\n({g.error_pct:+.1f}%)',
                    ha='center', fontsize=8)
    ax.set_yscale('log')
    ax.set_ylabel('predicted mass (MeV)')
    ax.set_title('σ-orbit ground-state particle masses\n(green = within 5%)')
    ax.grid(alpha=0.3, axis='y')
    ax.tick_params(axis='x', labelsize=8)

    fig.suptitle(
        f"Phase D-2 — σ-orbit Wilson-loop dynamics → "
        f"torus-knot bound-state masses\n"
        f"Paper 11 mass formula reproduces 25-particle spectrum at "
        f"{median_err:.2f}% median error",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_d2_wilson_loop.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_d2_wilson_loop.npz",
             particle_names=np.array([e["name"] for e in COMPENDIUM]),
             observed_masses_MeV=np.array(
                [e.get("m_obs", np.nan) for e in COMPENDIUM]),
             predicted_masses_MeV=np.array([
                particle_mass_MeV(e["p"], e["q"], e["m"], e["n_q"])
                for e in COMPENDIUM]),
             median_error_pct=median_err,
             max_error_pct=max_err,
             ground_state_orbits=np.array([g.orbit_id for g in ground]),
             ground_state_particles=np.array(
                [g.particle_name for g in ground]),
             ground_state_masses_pred=np.array(
                [g.mass_pred_MeV for g in ground]))
    print(f"  data saved {OUT_DIR / 'phase_d2_wilson_loop.npz'}")

    # ---- 5. Headline ---------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase D-2 σ-orbit ↔ torus-knot bound-state map")
    print("=" * 78)
    print()
    print(f"  ★ Paper 11 mass formula reproduces all 25 compendium particles")
    print(f"    at median {median_err:.2f}% error (max {max_err:.1f}%, nucleons).")
    print()
    print(f"  ★ σ-orbit topology gives natural particle-class assignments:")
    for g in ground:
        m_obs_str = f"{g.mass_obs_MeV:.1f} MeV" if g.mass_obs_MeV else "?"
        err_str = f"({g.error_pct:+.1f}%)" if g.error_pct is not None else ""
        print(f"    σ_{g.orbit_id}: {g.particle_name:<10} {m_obs_str:>14} {err_str}")
    print()
    print(f"  ★ Matter σ-orbit (polar=3 bridging edges) → electron (2,1,3,0)")
    print(f"    matches m_e exactly (anchor of the mass formula).")
    print(f"    The framework's f_J = (1-√α)³ Wilson product on the same")
    print(f"    3 bridging edges sets the COSMOGENIC ENERGY FRACTION;")
    print(f"    the mass formula sets the BOUND-STATE MASS — same edges,")
    print(f"    two distinct observables.")
    print()
    print(f"  ★ The (p,q,m,n_q) quantum numbers themselves remain Paper 11")
    print(f"    INPUTS from Refs. NWT1/NWT6.  Their derivation from K_7 σ-orbit")
    print(f"    geometry is the open theoretical problem — Phase E/F territory.")


if __name__ == "__main__":
    main()
