"""Steane-code reframe — Check 5: syndrome → Fano-point map (substrate reading).

For each outside-centralizer walk W:
  1. Compute Check 3 direction-sensitive Pauli (x_part, z_part).
  2. Compute X-syndrome (z-part anti-commutation with X-stabs, 3 bits) and
     Z-syndrome (x-part anti-commutation with Z-stabs, 3 bits).
  3. Each 3-bit syndrome (s_0, s_1, s_2) is interpreted as a binary index
     (s_0 + 2 s_1 + 4 s_2) ∈ {0..7}; nonzero indices identify a Hamming-
     position qubit = Fano point.
  4. Each particle thus has an (X-error qubit, Z-error qubit) pair —
     two Fano points marking "where the matter component is concentrated"
     and "where the antimatter component is concentrated."
  5. Tabulate by (p, q) class, by sector, and check matter/antimatter
     swap under walk reversal.

NWT's hypothesis: the syndrome → Fano point map identifies the qubit
on which the particle's wavefunction is concentrated.  If true,
particles in the same (p, q) class share their syndrome Fano points,
and matter/antimatter walks should have swapped X-error/Z-error
Fano points under H⊗7 = walk reversal.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_check5_syndrome_fano_map.py
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


def walk_direction_pauli(walk: list[int]) -> tuple[tuple, tuple]:
    x = [0]*7
    z = [0]*7
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i+1]
        d = (b - a) % 7
        qb = VERTEX_TO_QUBIT[b]
        if d in QR_DIRECTIONS: x[qb] ^= 1
        elif d in NR_DIRECTIONS: z[qb] ^= 1
    return tuple(x), tuple(z)


def x_syndrome(z_part: tuple) -> tuple:
    """3-bit syndrome from Z-part anti-commutation with X-stabilizers."""
    return tuple(
        sum(z_part[i]*g[i] for i in range(7)) % 2 for g in X_STAB_GENS
    )


def z_syndrome(x_part: tuple) -> tuple:
    """3-bit syndrome from X-part anti-commutation with Z-stabilizers."""
    return tuple(
        sum(x_part[i]*g[i] for i in range(7)) % 2 for g in Z_STAB_GENS
    )


def syndrome_to_qubit(syn: tuple) -> int | None:
    """3-bit Hamming syndrome (s_0, s_1, s_2) → qubit index 0..6 or None."""
    idx = syn[0] + 2*syn[1] + 4*syn[2]
    return idx - 1 if idx > 0 else None


def qubit_to_label(q: int | None) -> str:
    if q is None: return "I"
    return VERTEX_LABEL[QUBIT_TO_VERTEX[q]]


def reverse_walk(walk: list[int]) -> list[int]:
    return list(reversed(walk))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("Steane Check 5 — syndrome → Fano-point map")
    print("=" * 78)
    print()
    print("Each walk's direction-sensitive Pauli has 2 syndromes:")
    print("  X-error qubit = Fano point where 'matter component' is concentrated")
    print("  Z-error qubit = Fano point where 'antimatter component' is concentrated")
    print()
    print("Hypothesis (NWT): syndrome → Fano-point identifies the qubit on which")
    print("the particle's wavefunction is concentrated.")
    print()

    walks = bfs_shortest_walks(max_length=25)

    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks: continue
        walk = walks[key]
        x_part, z_part = walk_direction_pauli(walk)
        synd_X = x_syndrome(z_part)
        synd_Z = z_syndrome(x_part)
        q_X = syndrome_to_qubit(synd_X)
        q_Z = syndrome_to_qubit(synd_Z)

        # Reverse walk
        rwalk = reverse_walk(walk)
        rx, rz = walk_direction_pauli(rwalk)
        rsynd_X = x_syndrome(rz)
        rsynd_Z = z_syndrome(rx)
        rq_X = syndrome_to_qubit(rsynd_X)
        rq_Z = syndrome_to_qubit(rsynd_Z)

        rows.append({**entry, "walk": walk,
                      "syn_X": synd_X, "syn_Z": synd_Z,
                      "q_X": q_X, "q_Z": q_Z,
                      "rsyn_X": rsynd_X, "rsyn_Z": rsynd_Z,
                      "rq_X": rq_X, "rq_Z": rq_Z})

    # ---- Per-particle table -------------------------------------------
    print(f"{'particle':<10} {'n_q':<3} {'(p,q)':<8} {'X-syn':<5} {'Z-syn':<5} "
           f"{'X-err':<7} {'Z-err':<7}  | "
           f"{'rev X-syn':<10} {'rev Z-syn':<10} {'rev X-err':<10} {'rev Z-err'}")
    print("-" * 120)
    for r in rows:
        sx = ''.join(str(b) for b in r['syn_X'])
        sz = ''.join(str(b) for b in r['syn_Z'])
        rsx = ''.join(str(b) for b in r['rsyn_X'])
        rsz = ''.join(str(b) for b in r['rsyn_Z'])
        print(f"{r['name']:<10} {r['n_q']:<3} "
               f"({r['p']:>2},{r['q']:>2})  {sx:<5} {sz:<5} "
               f"{qubit_to_label(r['q_X']):<7} {qubit_to_label(r['q_Z']):<7}  | "
               f"{rsx:<10} {rsz:<10} "
               f"{qubit_to_label(r['rq_X']):<10} {qubit_to_label(r['rq_Z'])}")
    print()

    # ---- Hypothesis 1: matter/antimatter syndromes swap under reversal --
    print("=" * 78)
    print("HYPOTHESIS 1: walk reversal swaps X-error ↔ Z-error syndromes")
    print("=" * 78)
    print()
    swap_count = 0
    same_count = 0
    for r in rows:
        # Reverse walk's X-error should equal forward walk's Z-error
        # and Reverse walk's Z-error should equal forward walk's X-error
        swap_match = (r['rq_X'] == r['q_Z']) and (r['rq_Z'] == r['q_X'])
        if swap_match: swap_count += 1
        else: same_count += 1
        if not swap_match:
            print(f"  ✗ {r['name']:<10} fwd=({qubit_to_label(r['q_X'])},{qubit_to_label(r['q_Z'])}) "
                   f"rev=({qubit_to_label(r['rq_X'])},{qubit_to_label(r['rq_Z'])})")
    print()
    print(f"  Swap-match (rev syndrome = swap of fwd): {swap_count}/{len(rows)}")
    print()

    # ---- Hypothesis 2: same (p,q) → same syndrome ---------------------
    print("=" * 78)
    print("HYPOTHESIS 2: particles in same (|p|,|q|) class share syndrome")
    print("=" * 78)
    print()
    by_pq = defaultdict(list)
    for r in rows:
        by_pq[(abs(r['p']), abs(r['q']))].append(r)
    print(f"  (|p|,|q|) class breakdown:")
    print(f"  {'class':<10} {'particles':<40} {'X-errs':<10} {'Z-errs':<10}")
    print("  " + "-"*70)
    for key in sorted(by_pq.keys()):
        entries = by_pq[key]
        names = ','.join(e['name'] for e in entries)
        x_errs = set(qubit_to_label(e['q_X']) for e in entries)
        z_errs = set(qubit_to_label(e['q_Z']) for e in entries)
        print(f"  {str(key):<10} {names:<40} {str(x_errs):<10} {str(z_errs):<10}")
    print()

    # ---- Hypothesis 3: (p,q) → Fano-point function exists --------------
    print("=" * 78)
    print("HYPOTHESIS 3: (|p|,|q|) → Fano-point function?")
    print("=" * 78)
    print()
    pq_to_xerr = {}
    pq_to_zerr = {}
    for key, entries in by_pq.items():
        x_errs = set(e['q_X'] for e in entries)
        z_errs = set(e['q_Z'] for e in entries)
        # Function if a single value per key
        if len(x_errs) == 1:
            pq_to_xerr[key] = next(iter(x_errs))
        if len(z_errs) == 1:
            pq_to_zerr[key] = next(iter(z_errs))
    print(f"  classes with unique X-error: {len(pq_to_xerr)}/{len(by_pq)}")
    print(f"  classes with unique Z-error: {len(pq_to_zerr)}/{len(by_pq)}")
    print()

    # All particles share syndrome by (p,q)? If yes, function exists.
    function_x = len(pq_to_xerr) == len(by_pq)
    function_z = len(pq_to_zerr) == len(by_pq)
    print(f"  Pure function (p,q) → X-err Fano point exists: {function_x}")
    print(f"  Pure function (p,q) → Z-err Fano point exists: {function_z}")
    print()

    # Print the function
    if function_x:
        print("  (|p|,|q|) → X-error Fano point:")
        for key in sorted(by_pq.keys()):
            print(f"    ({key[0]},{key[1]}) → {qubit_to_label(pq_to_xerr[key])}")
    print()

    # ---- Hypothesis 4: Fano-point distribution by sector --------------
    print("=" * 78)
    print("HYPOTHESIS 4: Fano-point distribution by carrier-knot sector")
    print("=" * 78)
    print()
    sector_name = {0: 'leptons (unknot)', 2: 'mesons (Hopf)',
                    3: 'hyperons (trefoil)', 5: 'nucleons (cinquefoil)'}
    for nq in sorted(set(r['n_q'] for r in rows)):
        sec = sector_name.get(nq, f'n_q={nq}')
        entries = [r for r in rows if r['n_q'] == nq]
        x_errs = Counter(qubit_to_label(r['q_X']) for r in entries)
        z_errs = Counter(qubit_to_label(r['q_Z']) for r in entries)
        print(f"  {sec}: {len(entries)} walks")
        print(f"    X-err distribution: {dict(x_errs)}")
        print(f"    Z-err distribution: {dict(z_errs)}")
    print()

    # ---- Hypothesis 5: syndrome-pair structure ------------------------
    print("=" * 78)
    print("HYPOTHESIS 5: syndrome-pair (X-err, Z-err) substrate structure")
    print("=" * 78)
    print()
    pair_counts = Counter((qubit_to_label(r['q_X']), qubit_to_label(r['q_Z']))
                          for r in rows)
    print(f"  All distinct (X-err, Z-err) pairs in compendium:")
    for pair, cnt in sorted(pair_counts.items()):
        print(f"    {pair}: {cnt} walks")
    print()

    # Group by symmetric pair {a, b}
    print(f"  Grouped by unordered set {{X-err, Z-err}}:")
    unordered = defaultdict(int)
    for r in rows:
        s = frozenset([qubit_to_label(r['q_X']), qubit_to_label(r['q_Z'])])
        unordered[s] += 1
    for s in sorted(unordered.keys(), key=lambda x: sorted(x)):
        print(f"    {{ {', '.join(sorted(s))} }}: {unordered[s]} walks")
    print()


if __name__ == "__main__":
    main()
