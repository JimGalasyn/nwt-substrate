"""Bogoliubov Phase K-2 — test NWT's Z_3 chirality jump rule.

NWT 2026-05-20 ~4:30pm structural proposal:

  σ-orbit Z_3 chirality assignment (in our σ_# labels):
    identity: {σ_0, σ_1, σ_2}  (polar bridging × 2 + τ-fixed cross)
    ω:        {σ_3, σ_5}        (E-triangle + cross-twist+)
    ω²:       {σ_4, σ_6}        (F-triangle + cross-twist−)

  Rule: carrier-knot crossings = # ω↔ω² transitions in walk's σ-sequence
  (excluding identity-involving transitions).

This is the substrate-canonical "Z_3 chirality jump" count.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/nwt_bogoliubov_phase_k2_z3_chirality.py
"""
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from nwt_substrate.condensate import (
    bfs_shortest_walks, walk_sigma_sequence,
)
from nwt_substrate.particles.compendium import COMPENDIUM


OUT_DIR = Path(__file__).parent / "phase_k2_outputs"
OUT_DIR.mkdir(exist_ok=True, parents=True)


# NWT chirality assignment (in our σ_# labels)
CHIRALITY = {
    0: "id",   # polar matter (P-E)
    1: "id",   # polar Z_2 partner (P-F)
    2: "id",   # τ-fixed cross (Z_3-fixed)
    3: "ω",    # E-triangle (intra-A)
    4: "ω²",   # F-triangle (intra-B)
    5: "ω",    # cross-twist+ (Z_3 forward)
    6: "ω²",   # cross-twist− (Z_3 backward)
}


def count_chirality_jumps(walk: list[int],
                           variant: str = "strict") -> int:
    """Count chirality transitions in the walk's σ-orbit sequence.

    Variants:
      - "strict":  only ω↔ω² transitions (NWT's canonical proposal)
      - "any":     any chirality change (including id↔nonid)
      - "asym":    only ω→ω² (one direction)
      - "withid":  ω↔ω² PLUS id↔nonid (combined)
    """
    sig_seq = walk_sigma_sequence(walk)
    chir_seq = [CHIRALITY[s] for s in sig_seq]
    n = len(chir_seq)
    jumps = 0
    for i in range(n):
        cur = chir_seq[i]
        nxt = chir_seq[(i + 1) % n]
        if variant == "strict":
            if {cur, nxt} == {"ω", "ω²"}:
                jumps += 1
        elif variant == "any":
            if cur != nxt:
                jumps += 1
        elif variant == "asym":
            if cur == "ω" and nxt == "ω²":
                jumps += 1
        elif variant == "withid":
            if cur != nxt and not (cur == "id" and nxt == "id"):
                jumps += 1
    return jumps


def main() -> None:
    print("=" * 78)
    print("BOGOLIUBOV PHASE K-2 — Z_3 chirality jump rule for n_q")
    print("=" * 78)
    print()
    print("Chirality assignment (our labels):")
    for sid in range(7):
        print(f"  σ_{sid}: {CHIRALITY[sid]}")
    print()

    walks = bfs_shortest_walks(max_length=25)
    print(f"Loaded {len(walks)} (|p|, |q|) shortest walks")
    print()

    # For each compendium walk, compute jump counts under each variant
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        sig_seq = walk_sigma_sequence(walk)
        chir_seq = [CHIRALITY[s] for s in sig_seq]
        rows.append({
            **entry,
            "walk": walk,
            "sig_seq": sig_seq,
            "chir_seq": chir_seq,
            "jumps_strict": count_chirality_jumps(walk, "strict"),
            "jumps_any": count_chirality_jumps(walk, "any"),
            "jumps_asym": count_chirality_jumps(walk, "asym"),
            "jumps_withid": count_chirality_jumps(walk, "withid"),
        })

    print(f"  {'particle':<12} {'n_q':<4} {'walk_L':<7} "
          f"{'strict':<8} {'any':<6} {'asym':<6} {'withid':<8} "
          f"{'σ_seq':<35}")
    print("  " + "-" * 90)
    for r in rows:
        sig_str = ','.join(str(s) for s in r['sig_seq'])
        if len(sig_str) > 33:
            sig_str = sig_str[:30] + '...'
        print(f"  {r['name']:<12} {r['n_q']:<4} {len(r['walk'])-1:<7} "
              f"{r['jumps_strict']:<8} {r['jumps_any']:<6} "
              f"{r['jumps_asym']:<6} {r['jumps_withid']:<8} {sig_str}")
    print()

    # Match analysis
    print("=" * 78)
    print("MATCH ANALYSIS — does any variant reproduce Paper 11 n_q?")
    print("=" * 78)
    print()
    for variant in ["strict", "any", "asym", "withid"]:
        key = f"jumps_{variant}"
        n_q_obs = np.array([r["n_q"] for r in rows])
        n_q_pred = np.array([r[key] for r in rows])
        r_corr = np.corrcoef(n_q_obs, n_q_pred)[0, 1]
        exact_match = sum(1 for r in rows if r[key] == r["n_q"])
        print(f"  variant '{variant}':")
        print(f"    exact matches: {exact_match}/{len(rows)} "
              f"({exact_match/len(rows)*100:.0f}%)")
        print(f"    Pearson r: {r_corr:+.3f}")
        # Stratified by sector
        by_nq = defaultdict(list)
        for r in rows:
            by_nq[r["n_q"]].append(r[key])
        print(f"    by sector (n_q → jumps):")
        for nq in sorted(by_nq.keys()):
            sector = {0:'leptons', 2:'mesons', 3:'hyperons',
                       5:'nucleons'}.get(nq, '?')
            vals = sorted(set(by_nq[nq]))
            print(f"      n_q={nq} ({sector}): jumps = {vals}")
        print()

    # ---- Detailed proton + electron analysis ---------------------------
    print("=" * 78)
    print("DETAILED — proton (n_q=5) + electron (n_q=0) chirality sequences")
    print("=" * 78)
    print()
    for name in ["e-", "p", "Sigma+"]:
        r = next((x for x in rows if x["name"] == name), None)
        if r is None: continue
        print(f"  {name} (n_q={r['n_q']}):")
        print(f"    walk: {'→'.join(str(v) for v in r['walk'])}")
        print(f"    σ-orbits: {r['sig_seq']}")
        print(f"    chirality: {r['chir_seq']}")
        print(f"    transitions (incl wrap):")
        n = len(r['chir_seq'])
        for i in range(n):
            cur = r['chir_seq'][i]
            nxt = r['chir_seq'][(i + 1) % n]
            mark = ""
            if {cur, nxt} == {"ω", "ω²"}:
                mark = " ★ strict jump"
            elif cur != nxt:
                mark = "   (other change)"
            print(f"      σ_{r['sig_seq'][i]}({cur}) → "
                  f"σ_{r['sig_seq'][(i+1)%n]}({nxt}){mark}")
        print()

    # ---- Plot -----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, variant in zip(axes.flatten(),
                            ["strict", "any", "asym", "withid"]):
        key = f"jumps_{variant}"
        sector_colors = {0:'C2', 2:'C0', 3:'C1', 5:'C3'}
        for r in rows:
            c = sector_colors.get(r["n_q"], 'gray')
            ax.scatter(r[key], r["n_q"], c=c, s=80, alpha=0.6,
                       edgecolor='k')
            ax.annotate(r["name"], (r[key], r["n_q"]),
                        xytext=(4, 4), textcoords='offset points',
                        fontsize=7, alpha=0.6)
        lims = [-0.5, max(max(r[key] for r in rows), 6) + 1]
        ax.plot(lims, lims, 'k--', alpha=0.4, label='1:1')
        ax.set_xlabel(f'chirality jumps ({variant})')
        ax.set_ylabel('n_q (Paper 11)')
        r_v = np.corrcoef(
            [r["n_q"] for r in rows],
            [r[key] for r in rows])[0, 1]
        ax.set_title(f'{variant}: r = {r_v:+.3f}')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle(
        "Phase K-2 — NWT's Z_3 chirality jump rule\n"
        "(4 counting variants tested vs Paper 11 n_q)",
        fontsize=11, y=1.005,
    )
    fig.tight_layout()
    fig_path = OUT_DIR / "phase_k2_z3_chirality.png"
    fig.savefig(fig_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure saved {fig_path}")

    np.savez(OUT_DIR / "phase_k2_z3_chirality.npz",
             names=np.array([r["name"] for r in rows]),
             n_q=np.array([r["n_q"] for r in rows]),
             jumps_strict=np.array([r["jumps_strict"] for r in rows]),
             jumps_any=np.array([r["jumps_any"] for r in rows]),
             jumps_asym=np.array([r["jumps_asym"] for r in rows]),
             jumps_withid=np.array([r["jumps_withid"] for r in rows]))
    print(f"  data saved {OUT_DIR / 'phase_k2_z3_chirality.npz'}")

    # ---- Headline -------------------------------------------------------
    print()
    print("=" * 78)
    print("HEADLINE — Z_3 chirality jump rule test")
    print("=" * 78)
    print()
    for variant in ["strict", "any", "asym", "withid"]:
        key = f"jumps_{variant}"
        exact = sum(1 for r in rows if r[key] == r["n_q"])
        print(f"  {variant:<8} variant: {exact}/{len(rows)} exact match "
              f"({exact/len(rows)*100:.0f}%)")
    print()
    # Best variant
    best = max(["strict", "any", "asym", "withid"],
                key=lambda v: sum(1 for r in rows
                                   if r[f"jumps_{v}"] == r["n_q"]))
    print(f"  Best variant: '{best}'")


if __name__ == "__main__":
    main()
