"""Steane-code reframe — Check 3: direction-sensitive walk → Pauli map.

For each walk edge (v_i → v_{i+1}) with signed step d = (v_{i+1} - v_i) mod 7:
  - If d ∈ QR = {1, 2, 4}:  contribute X_{v_{i+1}} to the walk's Pauli
  - If d ∈ NR = {3, 5, 6}:  contribute Z_{v_{i+1}} to the walk's Pauli

Walk Pauli = product over edges.  Then check centralizer + reduce to logical
coset.  Compare matter walk (forward) vs antimatter walk (reverse).

Hypothesis (framework's strongest):
  - Pure-QR walks (e.g. proton)        → X⊗7 → X_L coset
  - Pure-NR walks (e.g. antiproton)    → Z⊗7 → Z_L coset
  - Walk reversal ↔ H⊗7 conjugation ↔ X_L ↔ Z_L coset swap

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_check3_direction_sensitive_pauli.py
"""
from __future__ import annotations

from collections import defaultdict, Counter
from pathlib import Path

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.particles.compendium import COMPENDIUM


VERTEX_TO_FANO = {
    0: (1, 1, 1), 1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1),
    4: (0, 1, 1), 5: (1, 0, 1), 6: (1, 1, 0),
}
VERTEX_LABEL = {0: 'P', 1: 'E1', 2: 'E2', 3: 'E3', 4: 'F1', 5: 'F2', 6: 'F3'}


def std_position(point: tuple) -> int:
    return 4 * point[0] + 2 * point[1] + 1 * point[2]


VERTEX_TO_QUBIT = {v: std_position(VERTEX_TO_FANO[v]) - 1 for v in range(7)}
QUBIT_TO_VERTEX = {q: v for v, q in VERTEX_TO_QUBIT.items()}

QR_DIRECTIONS = {1, 2, 4}
NR_DIRECTIONS = {3, 5, 6}


# ---------------------------------------------------------------------------
# Steane [[7,1,3]] stabilizers + logical operators
# ---------------------------------------------------------------------------

def hamming_parity_check():
    H = np.zeros((3, 7), dtype=int)
    for k in range(3):
        for q in range(7):
            pos = q + 1
            H[k, q] = (pos >> k) & 1
    return H


HAMMING = hamming_parity_check()
X_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]
Z_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]


def in_x_stab_subspace(v: tuple) -> bool:
    target = np.array(v, dtype=int)
    for mask in range(8):
        x = np.zeros(7, dtype=int)
        for k in range(3):
            if (mask >> k) & 1: x ^= np.array(X_STAB_GENS[k])
        if np.array_equal(x, target): return True
    return False


def in_z_stab_subspace(v: tuple) -> bool:
    return in_x_stab_subspace(v)  # CSS symmetric


def commutes_with_all_z_stabs(x_support: tuple) -> bool:
    """X-only Pauli on `x_support` commutes with all Z-stab generators iff
    |x_support ∩ each Z-stab support| is even.
    """
    for g in Z_STAB_GENS:
        if sum(x_support[i]*g[i] for i in range(7)) % 2 != 0:
            return False
    return True


def commutes_with_all_x_stabs(z_support: tuple) -> bool:
    for g in X_STAB_GENS:
        if sum(z_support[i]*g[i] for i in range(7)) % 2 != 0:
            return False
    return True


def in_centralizer(x_part: tuple, z_part: tuple) -> bool:
    """Pauli (X-part, Z-part) commutes with all 6 stabilizers."""
    return (commutes_with_all_z_stabs(x_part)
             and commutes_with_all_x_stabs(z_part))


def logical_coset(x_part: tuple, z_part: tuple) -> str:
    """Reduce in-centralizer Pauli to one of {I_L, X_L, Y_L, Z_L}."""
    x_bit = 0 if in_x_stab_subspace(x_part) else 1
    z_bit = 0 if in_z_stab_subspace(z_part) else 1
    return {(0,0):"I_L", (1,0):"X_L", (0,1):"Z_L", (1,1):"Y_L"}[(x_bit, z_bit)]


# ---------------------------------------------------------------------------
# Direction-sensitive walk → Pauli
# ---------------------------------------------------------------------------

def walk_direction_pauli(walk: list[int]) -> tuple[tuple, tuple]:
    """Returns (X-part, Z-part) as 7-bit tuples indexed by qubit.

    For each step (v_i → v_{i+1}):
      d = (v_{i+1} - v_i) mod 7
      if d ∈ QR:  toggle X-bit on qubit(v_{i+1})
      if d ∈ NR:  toggle Z-bit on qubit(v_{i+1})
    """
    x = [0]*7
    z = [0]*7
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i+1]
        d = (b - a) % 7
        qb = VERTEX_TO_QUBIT[b]
        if d in QR_DIRECTIONS:
            x[qb] ^= 1
        elif d in NR_DIRECTIONS:
            z[qb] ^= 1
        else:
            raise ValueError(f"step {a}→{b} has d=0 (self-loop)")
    return tuple(x), tuple(z)


def reverse_walk(walk: list[int]) -> list[int]:
    """Reverse a closed walk (CPT analog at walk level)."""
    return list(reversed(walk))


def hadamard_conjugate(x_part: tuple, z_part: tuple) -> tuple:
    """H⊗7 conjugation: swap X-part and Z-part."""
    return z_part, x_part


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("Steane Check 3 — direction-sensitive walk → Pauli")
    print("=" * 78)
    print()
    print("Rule: QR-step (d ∈ {1,2,4}) → X on destination; NR-step (d ∈ {3,5,6}) → Z on destination")
    print()
    print("Hypothesis:")
    print("  matter walks  → X-rich Paulis (predominantly X_L coset)")
    print("  antimatter walks → Z-rich Paulis (predominantly Z_L coset)")
    print("  walk reversal ↔ H⊗7 conjugation (swap X-part ↔ Z-part)")
    print()

    walks = bfs_shortest_walks(max_length=25)

    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks: continue
        walk = walks[key]
        x_part, z_part = walk_direction_pauli(walk)
        cent = in_centralizer(x_part, z_part)
        coset = logical_coset(x_part, z_part) if cent else "—"

        # Reverse walk
        rwalk = reverse_walk(walk)
        rx_part, rz_part = walk_direction_pauli(rwalk)
        rcent = in_centralizer(rx_part, rz_part)
        rcoset = logical_coset(rx_part, rz_part) if rcent else "—"

        # H⊗7 conjugate of forward Pauli
        hx, hz = hadamard_conjugate(x_part, z_part)
        # Compare: does reverse Pauli match H⊗7 conjugate (mod stabilizer)?
        # (Without stabilizer reduction, we compare directly first.)
        cpt_exact_match = (rx_part == hx and rz_part == hz)

        # Modulo stabilizer match: are H⊗7(P_fwd) and P_rev in same coset?
        hcoset = logical_coset(hx, hz) if in_centralizer(hx, hz) else "—"
        cpt_coset_match = (cent and rcent and rcoset == hcoset)

        rows.append({**entry, "walk": walk,
                      "x_part": x_part, "z_part": z_part,
                      "cent": cent, "coset": coset,
                      "rx_part": rx_part, "rz_part": rz_part,
                      "rcent": rcent, "rcoset": rcoset,
                      "hcoset": hcoset,
                      "cpt_exact_match": cpt_exact_match,
                      "cpt_coset_match": cpt_coset_match})

    # ---- Per-particle table -------------------------------------------
    print(f"{'particle':<10} {'n_q':<3} {'(p,q)':<8} {'L':<3} {'X-part':<10} "
           f"{'Z-part':<10} {'cent':<5} {'coset':<5}    "
           f"{'rev X':<10} {'rev Z':<10} {'rev cos':<8} {'CPT=H⊗7?'}")
    print("-" * 130)
    for r in rows:
        x_str = ''.join(str(b) for b in r['x_part'])
        z_str = ''.join(str(b) for b in r['z_part'])
        rx_str = ''.join(str(b) for b in r['rx_part'])
        rz_str = ''.join(str(b) for b in r['rz_part'])
        cpt = "exact" if r['cpt_exact_match'] else ("coset" if r['cpt_coset_match'] else "✗")
        print(f"{r['name']:<10} {r['n_q']:<3} "
               f"({r['p']:>2},{r['q']:>2})  {len(r['walk'])-1:<3} {x_str:<10} "
               f"{z_str:<10} {str(r['cent']):<5} {r['coset']:<5}    "
               f"{rx_str:<10} {rz_str:<10} {r['rcoset']:<8} {cpt}")
    print()

    # ---- Sector consistency --------------------------------------------
    print("=" * 78)
    print("SECTOR CONSISTENCY by carrier-knot")
    print("=" * 78)
    print()
    sector_name = {0: 'leptons (unknot)', 2: 'mesons (Hopf)',
                    3: 'hyperons (trefoil)', 5: 'nucleons (cinquefoil)'}
    by_nq = defaultdict(list)
    for r in rows:
        by_nq[r['n_q']].append(r)
    for nq in sorted(by_nq.keys()):
        sec = sector_name.get(nq, f'n_q={nq}')
        entries = by_nq[nq]
        cent_count = sum(1 for r in entries if r['cent'])
        cosets = Counter(r['coset'] for r in entries if r['cent'])
        out_cnt = len(entries) - cent_count
        print(f"  {sec}: {len(entries)} walks, {cent_count} in centralizer, "
               f"{out_cnt} outside")
        if cosets:
            print(f"    cosets: {dict(cosets)}")
        # forward / reverse coset pairs
        print(f"    (forward → reverse) coset pairs:")
        pairs = Counter((r['coset'], r['rcoset']) for r in entries)
        for (f_c, r_c), cnt in sorted(pairs.items()):
            print(f"      {f_c} → {r_c}: {cnt}")
    print()

    # ---- CPT = H⊗7 check ------------------------------------------------
    print("=" * 78)
    print("CPT-CHECK — does walk reversal act as H⊗7 (X ↔ Z swap)?")
    print("=" * 78)
    print()
    n_exact = sum(1 for r in rows if r['cpt_exact_match'])
    n_coset = sum(1 for r in rows if r['cpt_coset_match'])
    n_cent = sum(1 for r in rows if r['cent'] and r['rcent'])
    n_total = len(rows)
    print(f"  Exact equality:      {n_exact}/{n_total} walks satisfy P_reverse = H⊗7 P_forward exactly")
    print(f"  Coset equality:      {n_coset}/{n_total} walks satisfy [P_reverse] = [H⊗7 P_forward] mod stab")
    print(f"  Both in centralizer: {n_cent}/{n_total}")
    print()

    # ---- Symmetry check: do matter walks land in X_L and antimatter in Z_L?
    print("=" * 78)
    print("MATTER vs ANTIMATTER coset asymmetry")
    print("=" * 78)
    print()
    fwd_cosets = Counter(r['coset'] for r in rows if r['cent'])
    rev_cosets = Counter(r['rcoset'] for r in rows if r['rcent'])
    print(f"  Forward (matter)   coset distribution: {dict(fwd_cosets)}")
    print(f"  Reverse (antimatter) coset distribution: {dict(rev_cosets)}")
    print()
    print("  Per-pair check (forward, reverse):")
    pair_counts = Counter((r['coset'], r['rcoset']) for r in rows
                          if r['cent'] and r['rcent'])
    for (fc, rc), n in sorted(pair_counts.items()):
        print(f"    ({fc}, {rc}): {n}")
    print()

    # ---- Headline ------------------------------------------------------
    print("=" * 78)
    print("HEADLINE")
    print("=" * 78)
    print()
    print(f"  Walks in centralizer:        {sum(1 for r in rows if r['cent'])}/{n_total}")
    print(f"  Reverse walks in centralizer: {sum(1 for r in rows if r['rcent'])}/{n_total}")
    print(f"  Walks where reverse = H⊗7(forward) exactly: {n_exact}/{n_total}")
    print(f"  Walks where [reverse] = [H⊗7(forward)] mod stab: {n_coset}/{n_total}")
    print()
    print("  Proton walk (all-QR): expect X-part = 1111111, Z-part = 0000000.")
    r_p = next((r for r in rows if r['name'] == 'p'), None)
    if r_p:
        print(f"    actual: X-part = {''.join(str(b) for b in r_p['x_part'])}, "
               f"Z-part = {''.join(str(b) for b in r_p['z_part'])}, "
               f"coset = {r_p['coset']}")
        print(f"    reverse: X-part = {''.join(str(b) for b in r_p['rx_part'])}, "
               f"Z-part = {''.join(str(b) for b in r_p['rz_part'])}, "
               f"coset = {r_p['rcoset']}")
    print()


if __name__ == "__main__":
    main()
