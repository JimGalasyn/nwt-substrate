"""Bogoliubov Phase G-2 — sharpen all 4 substrate predictions.

NWT identified 4 substrate predictions from the 13 walk-realizable
(|p|, |q|) not in compendium:
  (2, 2)  light DM candidate    ~25 MeV
  (2, 3)  QR-dominant resonance ~160 MeV
  (3, 1)  NR-dominant resonance ~160 MeV
  (3, 3)  matter+CP hybrid       ~400 MeV  [SHARPENED in Phase G → n_q=3]

Phase G already refined (3, 3) → n_q = 3 via σ-orbit signature matching.
Phase G-2 applies the same methodology to the other 3 predictions.

Methodology (per Phase F-3 sector signatures):
  Sector       n_q   T/C    σ_4         σ_2-complete
  -------     ---   ----   --------    ------------
  lepton       0    2.25   varies      sometimes
  meson        2    0.97   varies      often
  hyperon      3    0.25   ALWAYS 0    often
  nucleon      5    4.00   ALWAYS 2    no

Each prediction is classified by closest-sector match, refined (n_q, m)
chosen to match NWT's mass estimate, and rendered with comparison.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_g2_all_predictions.py
"""
from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS
from nwt_substrate.condensate.orbit_winding import edge_winding_class
from nwt_substrate.particles.mass import paper6_mass_mev


OUT_DIR = Path(__file__).parent / "phase_g2_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


# Sector signatures from Phase F-3
SECTOR_TC = {
    "lepton (n_q=0)": 2.25,
    "meson (n_q=2)":  0.97,
    "hyperon (n_q=3)": 0.25,
    "nucleon (n_q=5)": 4.00,
}
SECTOR_NQ = {
    "lepton (n_q=0)": 0,
    "meson (n_q=2)":  2,
    "hyperon (n_q=3)": 3,
    "nucleon (n_q=5)": 5,
}


# NWT's 4 substrate predictions with their walks and target masses
PREDICTIONS = {
    (2, 2): {
        "walk": [0, 1, 2, 5, 1, 4, 0],
        "nwt_mass": 25.0,
        "nwt_type": "light DM candidate",
        "comment": "no NR-Hamilton, 33% QR / 67% NR mixed",
    },
    (2, 3): {
        "walk": [0, 1, 2, 3, 4, 0, 1, 4, 0],
        "nwt_mass": 160.0,
        "nwt_type": "QR-dominant resonance",
        "comment": "62% QR (NWT report), no NR-Hamilton",
    },
    (3, 1): {
        "walk": [0, 2, 4, 0, 2, 5, 1, 4, 0],
        "nwt_mass": 160.0,
        "nwt_type": "NR-dominant resonance",
        "comment": "38% QR (NWT report), no NR-Hamilton",
    },
    (3, 3): {
        "walk": [0, 1, 2, 3, 6, 2, 5, 1, 4, 0],
        "nwt_mass": 400.0,
        "nwt_type": "matter+CP hybrid baryon",
        "comment": "33% QR + 67% NR mixed, Hamilton-vertex + Hamilton-edge full",
    },
}


def edge_to_orbit(a: int, b: int) -> int:
    e_norm = (min(a, b), max(a, b))
    for oid, orbit in SIGMA_ORBITS.items():
        if e_norm in orbit["edges"]:
            return oid
    return -1


def walk_invariants(walk: list[int]) -> dict:
    L = len(walk) - 1
    orbit_counts = Counter()
    edges_used = set()
    for i in range(L):
        a, b = walk[i], walk[i + 1]
        oid = edge_to_orbit(a, b)
        orbit_counts[oid] += 1
        edges_used.add((min(a, b), max(a, b)))
    sigma_complete = {}
    for oid in range(7):
        orbit_edges = set(SIGMA_ORBITS[oid]["edges"])
        used = orbit_edges & edges_used
        sigma_complete[oid] = len(used) == 3
    triangle = orbit_counts.get(3, 0) + orbit_counts.get(4, 0)
    cross = (orbit_counts.get(2, 0) + orbit_counts.get(5, 0)
              + orbit_counts.get(6, 0))
    bridging = orbit_counts.get(0, 0) + orbit_counts.get(1, 0)
    tc = triangle / cross if cross > 0 else float('inf')
    n_distinct_verts = len(set(walk))
    return {
        "L": L,
        "orbit_counts": dict(orbit_counts),
        "edges_used": edges_used,
        "sigma_complete": sigma_complete,
        "triangle": triangle,
        "cross": cross,
        "bridging": bridging,
        "tc": tc,
        "sigma_4_count": orbit_counts.get(4, 0),
        "sigma_2_complete": sigma_complete[2],
        "n_distinct_verts": n_distinct_verts,
    }


def classify_sector(inv: dict) -> tuple[str, int]:
    """Return (sector_name, n_q) most consistent with the walk's invariants."""
    # Hard constraints
    has_sigma_4 = inv["sigma_4_count"] > 0
    has_full_sigma_4 = inv["sigma_4_count"] >= 2  # nucleon signature
    # T/C closest sector
    tc = inv["tc"]
    if math.isinf(tc):
        # No cross orbits — likely nucleon-like
        # Nucleon has triangle=4, cross=1 not 0; pure-triangle is extreme
        return "nucleon (n_q=5)", 5
    distances = {s: abs(tc - t) for s, t in SECTOR_TC.items()}
    tc_sector = min(distances, key=distances.get)

    # Apply hard constraints to override T/C if needed:
    if not has_sigma_4:
        # σ_4 = 0 → exclude nucleon (which always has σ_4=2)
        if tc_sector == "nucleon (n_q=5)":
            # downgrade to next-closest non-nucleon sector
            non_nucleon = {s: d for s, d in distances.items()
                           if s != "nucleon (n_q=5)"}
            tc_sector = min(non_nucleon, key=non_nucleon.get)
    if has_full_sigma_4:
        # σ_4 >= 2 strongly suggests nucleon or meson (mesons can have σ_4)
        # If T/C in (1, 5) range, this confirms; if very small T/C, conflict
        pass  # keep T/C-driven classification
    return tc_sector, SECTOR_NQ[tc_sector]


def m_for_target_mass(p: int, q: int, n_q: int,
                       target_mass: float, m_max: int = 30) -> tuple[int, float]:
    """Find m that gives mass closest to target via Paper 11 formula."""
    best_m = -1
    best_mass = 0.0
    best_err = float('inf')
    for m in range(p + 1, m_max):
        mass = paper6_mass_mev(p, q, m, n_q)
        if mass is None:
            continue
        err = abs(mass - target_mass)
        if err < best_err:
            best_err = err
            best_m = m
            best_mass = mass
    return best_m, best_mass


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE G-2 — sharpen all 4 substrate predictions")
    print("=" * 78)
    print()

    results = {}
    for pq, info in PREDICTIONS.items():
        walk = info["walk"]
        inv = walk_invariants(walk)
        sector_name, n_q_refined = classify_sector(inv)

        # Best m for both NWT's original (n_q=5) and Phase G2 refined n_q
        m_orig, mass_orig = m_for_target_mass(pq[0], pq[1], 5, info["nwt_mass"])
        m_refined, mass_refined = m_for_target_mass(
            pq[0], pq[1], n_q_refined, info["nwt_mass"])

        results[pq] = {
            **info,
            "inv": inv,
            "sector": sector_name,
            "n_q_refined": n_q_refined,
            "m_orig": m_orig, "mass_orig": mass_orig,
            "m_refined": m_refined, "mass_refined": mass_refined,
        }

    # ---- Print per-prediction analysis ----------------------------------
    for pq, r in results.items():
        inv = r['inv']
        print("=" * 78)
        print(f"PREDICTION (p, q) = {pq} — {r['nwt_type']}")
        print("=" * 78)
        print(f"  Walk: {'→'.join(str(v) for v in r['walk'])} (L={inv['L']})")
        print(f"  NWT initial: ~{r['nwt_mass']} MeV, {r['comment']}")
        print()
        oc = inv['orbit_counts']
        print(f"  σ-orbit composition: σ_0={oc.get(0,0)}, σ_1={oc.get(1,0)}, "
              f"σ_2={oc.get(2,0)}, σ_3={oc.get(3,0)}, σ_4={oc.get(4,0)}, "
              f"σ_5={oc.get(5,0)}, σ_6={oc.get(6,0)}")
        sig_comp = [f'σ_{o}' for o, c in inv['sigma_complete'].items() if c]
        print(f"  σ-orbits with COMPLETE coverage: "
              f"{', '.join(sig_comp) if sig_comp else 'none'}")
        print(f"  Triangle (σ_3+σ_4) = {inv['triangle']}, "
              f"Cross (σ_2+σ_5+σ_6) = {inv['cross']}, "
              f"T/C = {inv['tc']:.3f}")
        print(f"  Vertex coverage: {inv['n_distinct_verts']}/7")
        print()
        print(f"  ★ Phase G-2 sector match: {r['sector']}")
        print(f"  ★ Refined n_q estimate:    n_q = {r['n_q_refined']}")
        print()
        print(f"  Mass-formula fit to NWT target ~{r['nwt_mass']} MeV:")
        print(f"    with NWT's n_q=5: m = {r['m_orig']}, mass = "
              f"{r['mass_orig']:.1f} MeV")
        print(f"    with refined n_q={r['n_q_refined']}: m = {r['m_refined']}, "
              f"mass = {r['mass_refined']:.1f} MeV")
        print()

    # ---- Summary table --------------------------------------------------
    print("=" * 78)
    print("SUMMARY — Phase G-2 refinement of 4 substrate predictions")
    print("=" * 78)
    print()
    print(f"  {'(p, q)':<8} {'NWT initial':<22} "
          f"{'Phase G-2 refined':<22} {'best (n_q, m)':<14} "
          f"{'mass (MeV)':<11}")
    print("  " + "-" * 80)
    for pq, r in results.items():
        nwt_str = f"~{r['nwt_mass']:.0f} MeV"
        refined_sector = r['sector'].split('(')[0].strip()
        refined_str = f"{refined_sector} (n_q={r['n_q_refined']})"
        nq_m_str = f"({r['n_q_refined']}, {r['m_refined']})"
        mass_str = f"{r['mass_refined']:.1f}"
        print(f"  {str(pq):<8} {nwt_str:<22} {refined_str:<22} "
              f"{nq_m_str:<14} {mass_str:<11}")
    print()

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) T/C ratio of each prediction + sector lines
    ax = axes[0, 0]
    for s, tc in SECTOR_TC.items():
        ax.axvline(tc, ls='--', alpha=0.4)
        ax.text(tc, 0.7, s.split('(')[0].strip(), rotation=90,
                fontsize=9, ha='right', va='center', alpha=0.6)
    for pq, r in results.items():
        ax.scatter(r['inv']['tc'], list(results).index(pq) * 0.2 + 0.1,
                    s=300, alpha=0.8, edgecolor='k')
        ax.annotate(f'({pq[0]},{pq[1]})\n→ n_q={r["n_q_refined"]}',
                    (r['inv']['tc'], list(results).index(pq) * 0.2 + 0.1),
                    xytext=(8, 8), textcoords='offset points',
                    fontsize=9, fontweight='bold')
    ax.set_xscale('log')
    ax.set_xlabel('T/C ratio (log)')
    ax.set_yticks([])
    ax.set_title('Phase G-2: substrate predictions on T/C sector axis')
    ax.set_xlim(0.05, 10)
    ax.grid(alpha=0.3, which='both')

    # (b) σ-orbit composition heatmap (predictions only)
    ax = axes[0, 1]
    pq_list = list(results.keys())
    M = np.zeros((len(pq_list), 7))
    for i, pq in enumerate(pq_list):
        for o in range(7):
            M[i, o] = results[pq]['inv']['orbit_counts'].get(o, 0)
    im = ax.imshow(M, aspect='auto', cmap='YlOrRd', origin='lower')
    ax.set_yticks(range(len(pq_list)))
    ax.set_yticklabels([f'({pq[0]},{pq[1]})' for pq in pq_list])
    ax.set_xticks(range(7))
    ax.set_xticklabels([f'σ_{i}' for i in range(7)])
    for i in range(len(pq_list)):
        for j in range(7):
            v = int(M[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha='center', va='center',
                        fontsize=10,
                        color='white' if v > M.max() / 2 else 'black')
    ax.set_title('σ-orbit composition of substrate predictions')
    plt.colorbar(im, ax=ax, label='edge count')

    # (c) Mass spectra at refined n_q
    ax = axes[1, 0]
    m_range = list(range(3, 30))
    for pq, r in results.items():
        n_q = r['n_q_refined']
        masses = []
        ms_valid = []
        for m_pc in m_range:
            mass = paper6_mass_mev(pq[0], pq[1], m_pc, n_q)
            if mass is not None:
                masses.append(mass)
                ms_valid.append(m_pc)
        ax.plot(ms_valid, masses, 'o-',
                label=f'({pq[0]},{pq[1]}) n_q={n_q}', lw=2, ms=5)
        # mark the best-fit m
        ax.scatter(r['m_refined'], r['mass_refined'],
                    s=200, marker='*', edgecolor='k',
                    linewidth=2, zorder=10)
    # PDG anchors
    pdg_refs = [(105.66, "μ⁻"), (139.57, "π⁺"), (493.68, "K⁺"),
                 (938.27, "p"), (1115.7, "Λ"), (1672.5, "Ω⁻")]
    for m_val, name in pdg_refs:
        ax.axhline(m_val, color='gray', ls=':', alpha=0.4)
        ax.text(29.5, m_val, name, fontsize=8, va='center', color='gray')
    ax.set_yscale('log')
    ax.set_xlabel('m (phase-closure integer)')
    ax.set_ylabel('predicted mass (MeV)')
    ax.set_title('Mass spectra at refined n_q (★ = best fit to NWT mass)')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(alpha=0.3, which='both')

    # (d) Summary table
    ax = axes[1, 1]
    ax.axis('off')
    rows_table = [["(p, q)", "NWT init", "Phase G-2 refined", "mass (MeV)"]]
    for pq, r in results.items():
        rows_table.append([
            f'({pq[0]}, {pq[1]})',
            f"{r['nwt_mass']:.0f} MeV / {r['nwt_type'].split()[-2]}-like",
            f"n_q={r['n_q_refined']} ({r['sector'].split('(')[0].strip()})",
            f"{r['mass_refined']:.1f} (m={r['m_refined']})",
        ])
    table = ax.table(cellText=rows_table, loc='center', cellLoc='left',
                     colWidths=[0.10, 0.32, 0.32, 0.26])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.8)
    for j in range(4):
        table[(0, j)].set_facecolor('#cccccc')
    # highlight (3, 3) since it's the headline
    for j in range(4):
        table[(4, j)].set_facecolor('#d4f4d4')
    ax.set_title('Phase G-2 refinement summary', pad=15)

    fig.suptitle(
        f"Phase G-2 — Substrate prediction refinement for all 4 NWT candidates\n"
        f"(σ-orbit sector matching via Phase F-3 signatures)",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_g2_all_predictions.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_g2_all_predictions.npz",
             pq_list=np.array(list(results.keys())),
             walks=np.array(['-'.join(str(v) for v in r['walk'])
                              for r in results.values()]),
             nwt_masses=np.array([r['nwt_mass'] for r in results.values()]),
             refined_n_q=np.array([r['n_q_refined'] for r in results.values()]),
             refined_m=np.array([r['m_refined'] for r in results.values()]),
             refined_mass=np.array([r['mass_refined']
                                     for r in results.values()]),
             tc_ratios=np.array([r['inv']['tc'] for r in results.values()]))
    print(f"  data saved {OUT_DIR / 'phase_g2_all_predictions.npz'}")

    # ---- Headline -------------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase G-2 refinement of all 4 substrate predictions")
    print("=" * 78)
    print()
    for pq, r in results.items():
        print(f"  ({pq[0]}, {pq[1]}): n_q = {r['n_q_refined']} "
              f"({r['sector'].split('(')[0].strip()})  "
              f"mass = {r['mass_refined']:.1f} MeV (m = {r['m_refined']})")
    print()
    print(f"  ★ Sector distribution of substrate predictions:")
    from collections import defaultdict
    by_sector = defaultdict(list)
    for pq, r in results.items():
        by_sector[r['n_q_refined']].append(pq)
    for n_q in sorted(by_sector.keys()):
        preds = by_sector[n_q]
        sector = {0: 'lepton', 2: 'meson', 3: 'hyperon',
                   5: 'nucleon'}.get(n_q, '?')
        pq_str = ', '.join(f'({p[0]},{p[1]})' for p in preds)
        print(f"    {sector} (n_q={n_q}): {pq_str}  [{len(preds)} prediction(s)]")
    print()


if __name__ == "__main__":
    main()
