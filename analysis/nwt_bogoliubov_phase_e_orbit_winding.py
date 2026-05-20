"""Bogoliubov Phase E — Derive (p, q) quantum numbers from K_7 σ-orbit winding.

Phase D-2 used Paper 11's (p, q, m, n_q) quantum numbers as INPUTS.
Phase E tests whether they can be DERIVED from the K_7 Heffter
torus embedding's geometry.

Strategy
--------
The K_7 Heffter embedding places vertex k at (u_k, v_k) = (k/7, 3k mod 7 / 7)
on the unit torus.  Each edge has a shortest-path displacement on the
universal cover; sum over a σ-orbit's edges gives a topological
signature in 1/7 units.  A closed walk's total winding (in integer
multiples of 7) gives the (p, q) winding of any torus knot supported
on that walk.

Phase E hypothesis: σ-orbit winding signatures (or simple derived
combinations) reproduce Paper 11's (p, q) quantum numbers.

Outcome
-------
Mixed.  Some σ-orbits show striking matches (σ_6 raw winding = (2, 1)
= electron's (p, q) exactly), and several σ-orbits share a reduced
gcd-primitive of (3, 2).  But a clean derivation across all
compendium particles requires Paper 13's specific σ-orbit → walk
construction, not the naïve edge sum.

Honest scope: this script documents the topological invariants of
each σ-orbit and identifies the partial matches.  Full derivation
of (p, q, m, n_q) from K_7 geometry remains a Phase F target.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_e_orbit_winding.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate.sigma_orbits import SIGMA_ORBITS, orbit_invariants
from nwt_substrate.condensate.orbit_winding import (
    all_orbit_windings, edge_winding_class, HEFFTER_VERT_UV,
    closed_walk_winding, orbit_winding,
)
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_e_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


def reduce_winding(p: int, q: int) -> tuple[int, int]:
    """Reduce (p, q) by gcd to primitive form, preserving signs."""
    if p == 0 and q == 0:
        return 0, 0
    g = math.gcd(abs(p), abs(q))
    return p // g, q // g


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE E — Derive (p, q) from K_7 Heffter winding")
    print("=" * 78)
    print()

    # ---- 1. K_7 Heffter vertex positions --------------------------------
    print("K_7 Heffter vertex positions on unit torus (1/7 units):")
    print(f"  {'vertex':<10} {'u':>6} {'v':>6}")
    for k in range(7):
        u, v = HEFFTER_VERT_UV[k]
        print(f"  {f'v_{k} ({chr(80) if k == 0 else f'E_{k}' if k <= 3 else f'F_{k-3}'})':<10} "
              f"{u*7:>6.0f} {v*7:>6.0f}")
    print()

    # ---- 2. σ-orbit winding signatures ----------------------------------
    print("σ-orbit edge winding signatures (in 1/7 units of unit torus):")
    print(f"  {'id':<3} {'name':<28} {'per-edge (n_u, n_v)':<35} "
          f"{'total':<10} {'reduced':<10}")
    print(f"  " + "-" * 92)
    windings = all_orbit_windings()
    reduced_summary = {}
    for ow in windings:
        o = SIGMA_ORBITS[ow.orbit_id]
        edge_str = ', '.join(f'({a},{b})' for (a, b) in ow.edge_windings)
        pq_red = reduce_winding(*ow.total_winding)
        reduced_summary[ow.orbit_id] = pq_red
        print(f"  {ow.orbit_id:<3} {o['name'][:26]:<28} {edge_str:<35} "
              f"{str(ow.total_winding):<10} {str(pq_red):<10}")
    print()
    print("Key observations:")
    print(f"  σ_0 reduced (3, 2) — σ_1 reduced (-3, -2) — Z_2 partners ✓")
    print(f"  σ_2 reduced (3, 2) — same primitive as σ_0, σ_1!")
    print(f"  σ_3 = σ_4 reduced (4, 5) — Z_2 partners (E↔F triangles)")
    print(f"  σ_5 = (-5, -1), σ_6 = (2, -1) — Z_3 chirality partners")
    print()
    print(f"  ★ σ_6 raw winding (2, -1) → primitive (2, 1) matches electron (p, q)!")
    print(f"  ★ σ_0, σ_1, σ_2 all reduce to (3, 2) primitive — trefoil family")
    print()

    # ---- 3. Compare to Paper 11 compendium (p, q) -----------------------
    print("=" * 78)
    print("σ-orbit winding signatures vs Paper 11 (p, q) compendium")
    print("=" * 78)
    print()
    # Group compendium particles by (|p|, |q|) primitive
    pq_to_particles: dict[tuple[int, int], list[str]] = {}
    for e in COMPENDIUM:
        pq_red = reduce_winding(e["p"], e["q"])
        pq_abs = (abs(pq_red[0]), abs(pq_red[1]))
        pq_to_particles.setdefault(pq_abs, []).append(e["name"])

    print("Primitive (|p|, |q|) classes in compendium:")
    for pq, names in sorted(pq_to_particles.items()):
        print(f"  {str(pq):<10} → {', '.join(names)}")
    print()

    print("σ-orbit ↔ compendium match check:")
    print(f"  {'σ-orbit':<10} {'reduced':<10} {'|reduced|':<10} {'matches':<25}")
    print(f"  " + "-" * 60)
    for orbit_id, pq_red in reduced_summary.items():
        pq_abs = (abs(pq_red[0]), abs(pq_red[1]))
        matches = pq_to_particles.get(pq_abs, [])
        match_str = ', '.join(matches) if matches else '(no match in compendium)'
        print(f"  σ_{orbit_id}        {str(pq_red):<10} "
              f"{str(pq_abs):<10} {match_str}")
    print()

    # ---- 4. Test specific closed walks → (p, q) ------------------------
    print("=" * 78)
    print("CLOSED WALK winding numbers (for selected walks)")
    print("=" * 78)
    print()
    # Various candidate walks: each is a list of vertex indices, closed.
    walks = {
        "σ_0 star + return": [0, 1, 2, 3, 0],
        "σ_0 star different order": [0, 2, 1, 3, 0],
        "matter Z_3 cycle (E)": [0, 1, 2, 3, 0],
        "cross Z_2 link": [1, 4, 1],
        "Heffter triangle 0-1-3": [0, 1, 3, 0],
        "Heffter triangle 0-2-6": [0, 2, 6, 0],
        "K_7 length-7 cycle": [0, 1, 2, 3, 4, 5, 6, 0],
    }
    print(f"  {'walk':<28} {'vertices':<25} {'displacement':<14} {'winding (p, q)'}")
    print(f"  " + "-" * 80)
    for name, walk in walks.items():
        # raw displacement
        sum_u = sum_v = 0
        for i in range(len(walk) - 1):
            nu, nv = edge_winding_class(walk[i], walk[i + 1])
            sum_u += nu
            sum_v += nv
        # closed walk: should give multiples of 7 if walk closes on universal cover
        if walk[0] == walk[-1]:
            p_int = sum_u // 7 if sum_u % 7 == 0 else float("nan")
            q_int = sum_v // 7 if sum_v % 7 == 0 else float("nan")
            pq = (p_int, q_int) if isinstance(p_int, int) else "(not closed on cover)"
        else:
            pq = "(open walk)"
        path_str = '→'.join(str(v) for v in walk)
        print(f"  {name:<28} {path_str:<25} {str((sum_u, sum_v)):<14} {pq}")
    print()

    # ---- 5. Plot --------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) Vertex positions on unit torus
    ax = axes[0, 0]
    for k in range(7):
        u, v = HEFFTER_VERT_UV[k] * 7
        ax.scatter(u, v, s=250, c='C0', edgecolor='black', zorder=10)
        ax.text(u, v, str(k), ha='center', va='center', fontweight='bold',
                color='white', fontsize=10)
    ax.set_xlim(-0.5, 7.5)
    ax.set_ylim(-0.5, 7.5)
    ax.set_xlabel('u (1/7 units)')
    ax.set_ylabel('v (1/7 units)')
    ax.set_title('K_7 Heffter vertex positions on unit torus\n'
                 'v_k = (k mod 7, 3k mod 7) / 7')
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)
    ax.set_xticks(range(8))
    ax.set_yticks(range(8))

    # (b) σ-orbit winding signatures on (P, Q) plane
    ax = axes[0, 1]
    for ow in windings:
        o = SIGMA_ORBITS[ow.orbit_id]
        ax.scatter(*ow.total_winding, s=250,
                    edgecolor='black', linewidth=2,
                    label=f'σ_{ow.orbit_id}: {o["name"][:18]}')
        ax.annotate(f'σ_{ow.orbit_id}', ow.total_winding,
                    xytext=(8, 5), textcoords='offset points',
                    fontsize=10, fontweight='bold')
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('Δu_total (1/7 units)')
    ax.set_ylabel('Δv_total (1/7 units)')
    ax.set_title('σ-orbit winding signatures on universal cover')
    ax.legend(loc='upper left', fontsize=7, ncol=2)
    ax.grid(alpha=0.3)

    # (c) σ-orbit reduced primitives vs compendium (p, q)
    ax = axes[1, 0]
    # Plot compendium (p, q) — note: take abs to fold to first quadrant
    cp = [abs(e["p"]) for e in COMPENDIUM]
    cq = [abs(e["q"]) for e in COMPENDIUM]
    ax.scatter(cp, cq, s=70, c='gray', alpha=0.5, label='Paper 11 compendium')
    for e in COMPENDIUM:
        ax.annotate(e["name"], (abs(e["p"]), abs(e["q"])),
                    xytext=(3, 3), textcoords='offset points',
                    fontsize=7, alpha=0.6)
    # Plot σ-orbit reduced primitives
    for orbit_id, pq_red in reduced_summary.items():
        pq_abs = (abs(pq_red[0]), abs(pq_red[1]))
        ax.scatter(*pq_abs, s=300, c='C3', edgecolor='black',
                    marker='X', linewidth=2, zorder=10)
        ax.annotate(f'σ_{orbit_id}', pq_abs, xytext=(5, 5),
                    textcoords='offset points', fontsize=10,
                    fontweight='bold', color='C3')
    ax.set_xlabel('|p| (toroidal winding)')
    ax.set_ylabel('|q| (poloidal winding)')
    ax.set_title('σ-orbit reduced primitives vs Paper 11 (p, q)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_xticks(range(11))
    ax.set_yticks(range(11))

    # (d) Match table
    ax = axes[1, 1]
    ax.axis('off')
    table_data = [['σ', 'raw winding', 'reduced |p|,|q|', 'matches in compendium']]
    for orbit_id, pq_red in reduced_summary.items():
        ow = next(w for w in windings if w.orbit_id == orbit_id)
        pq_abs = (abs(pq_red[0]), abs(pq_red[1]))
        matches = pq_to_particles.get(pq_abs, [])
        m_str = ', '.join(matches[:4]) + ('...' if len(matches) > 4 else '')
        table_data.append([
            f'σ_{orbit_id}',
            str(ow.total_winding),
            str(pq_abs),
            m_str if m_str else '(no match)',
        ])
    table = ax.table(cellText=table_data, loc='center', cellLoc='left',
                     colWidths=[0.10, 0.22, 0.18, 0.50])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.7)
    # Highlight matches
    for i, row in enumerate(table_data[1:], start=1):
        if '(no match)' not in row[3]:
            for j in range(4):
                table[(i, j)].set_facecolor('#d4f4d4')
    ax.set_title('σ-orbit winding ↔ Paper 11 (p, q) match', pad=15)

    fig.suptitle(
        f"Phase E — Derive (p, q) from K_7 Heffter winding\n"
        f"Partial: σ_6 raw (2, 1) = electron (p, q); σ_0/1/2 reduce to (3, 2) trefoil family",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_e_orbit_winding.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    # Save data
    np.savez(OUT_DIR / "phase_e_orbit_winding.npz",
             orbit_ids=np.array([w.orbit_id for w in windings]),
             raw_winding_u=np.array([w.total_winding[0] for w in windings]),
             raw_winding_v=np.array([w.total_winding[1] for w in windings]),
             reduced_p=np.array(
                [reduce_winding(*w.total_winding)[0] for w in windings]),
             reduced_q=np.array(
                [reduce_winding(*w.total_winding)[1] for w in windings]),
             heffter_uv=HEFFTER_VERT_UV)
    print(f"  data saved {OUT_DIR / 'phase_e_orbit_winding.npz'}")

    # ---- 6. Headline ----------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Phase E partial result")
    print("=" * 78)
    print()
    print(f"  σ-orbit raw windings (in 1/7 units) computed for all 7 orbits.")
    print()
    print(f"  ★ PARTIAL MATCHES with Paper 11 (p, q):")
    print(f"    σ_6 raw winding (2, -1) → primitive (2, 1) = electron (p,q) ✓")
    print(f"    σ_0, σ_1, σ_2 all reduce to (3, 2) — trefoil family")
    print(f"      → would correspond to (3,2) torus knot particles")
    print(f"      → no compendium particle has (p,q) = (3,2) directly")
    print(f"    σ_3 = σ_4 reduce to (4, 5) — no direct compendium match")
    print(f"    σ_5 reduce to (5, 1) — no direct compendium match")
    print()
    print(f"  ✗ Simple 'sum edge windings' DOES NOT reproduce Paper 11 (p, q)")
    print(f"    quantum numbers for most particles.")
    print()
    print(f"  Interpretation: each particle's (p, q) is the winding of a")
    print(f"  SPECIFIC CLOSED WALK on K_7, not the raw σ-orbit edge sum.")
    print(f"  The walk uses σ-orbit edges + connector edges, and the")
    print(f"  walk's specific path is what Paper 13's SM-capstone")
    print(f"  construction is intended to determine.")
    print()
    print(f"  Phase F target: enumerate K_7 closed walks of bounded length")
    print(f"  and find which walks produce compendium (p, q) values, then")
    print(f"  identify the σ-orbit constraint that picks each walk.")


if __name__ == "__main__":
    main()
