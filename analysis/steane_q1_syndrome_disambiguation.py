"""Steane-code reframe — Q1: syndrome augmentation for full (p, q) encoding.

Check 5 established that the (X-err Fano, Z-err Fano) syndrome map on 16
compendium (|p|, |q|) classes is many-to-one — 7 cross-sector equivalences
collapse multiple distinct (p, q) classes onto a single syndrome pair.

Question (Jim, 2026-05-21 morning): does augmenting the syndrome with
additional substrate invariants — walk length L, σ-orbit composition
signature — restore an injective map (|p|, |q|) → invariant tuple?

This script:
  1. Computes the augmented invariant tuple for every compendium walk:
       (X-err Fano, Z-err Fano, L, σ-orbit signature)
     where σ-orbit signature is the multiset {orbit_id : edge_count}
     describing how the walk's edges distribute across the 7 σ-orbits.
  2. Tests which subsets of {syndrome_pair, L, σ-orbit signature} suffice
     to disambiguate every (|p|, |q|) cross-sector pair.
  3. Reports the MINIMAL sufficient invariant set per cross-sector pair.

Honest framing: the test is per-(|p|, |q|)-pair disambiguation, not full
injectivity of the map; multiple particles share a (|p|, |q|) class and so
share all substrate invariants by construction (e.g. p, n, Σ⁺ all have
(1,3)).

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q1_syndrome_disambiguation.py
"""
from __future__ import annotations

from collections import defaultdict, Counter
from itertools import combinations

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.condensate.walks import edge_to_orbit
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
    import numpy as np
    H = np.zeros((3, 7), dtype=int)
    for k in range(3):
        for q in range(7):
            H[k, q] = ((q + 1) >> k) & 1
    return H


HAMMING = hamming_parity_check()
X_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]
Z_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]


def walk_direction_pauli(walk: list[int]) -> tuple[tuple, tuple]:
    x = [0] * 7
    z = [0] * 7
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        d = (b - a) % 7
        qb = VERTEX_TO_QUBIT[b]
        if d in QR_DIRECTIONS:
            x[qb] ^= 1
        elif d in NR_DIRECTIONS:
            z[qb] ^= 1
    return tuple(x), tuple(z)


def x_syndrome(z_part: tuple) -> tuple:
    return tuple(sum(z_part[i] * g[i] for i in range(7)) % 2
                 for g in X_STAB_GENS)


def z_syndrome(x_part: tuple) -> tuple:
    return tuple(sum(x_part[i] * g[i] for i in range(7)) % 2
                 for g in Z_STAB_GENS)


def syndrome_to_qubit(syn: tuple) -> int | None:
    idx = syn[0] + 2 * syn[1] + 4 * syn[2]
    return idx - 1 if idx > 0 else None


def qubit_to_label(q: int | None) -> str:
    if q is None:
        return "I"
    return VERTEX_LABEL[QUBIT_TO_VERTEX[q]]


def sigma_orbit_signature(walk: list[int]) -> tuple[int, ...]:
    """Multiset of σ-orbit IDs visited by the walk, as a sorted-count vector
    of length 7 (one entry per σ-orbit). The signature is direction-
    insensitive: edge (a, b) and edge (b, a) hit the same σ-orbit.
    """
    counts = [0] * 7
    for i in range(len(walk) - 1):
        a, b = walk[i], walk[i + 1]
        oid = edge_to_orbit(a, b)
        if 0 <= oid < 7:
            counts[oid] += 1
    return tuple(counts)


def main() -> None:
    print("=" * 78)
    print("Steane Q1 — syndrome augmentation for full (p, q) encoding")
    print("=" * 78)
    print()

    walks = bfs_shortest_walks(max_length=25)

    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        x_part, z_part = walk_direction_pauli(walk)
        synd_X = x_syndrome(z_part)
        synd_Z = z_syndrome(x_part)
        q_X = syndrome_to_qubit(synd_X)
        q_Z = syndrome_to_qubit(synd_Z)
        L = len(walk) - 1
        sig = sigma_orbit_signature(walk)
        rows.append({
            **entry,
            "walk": walk,
            "L": L,
            "synd_pair": (qubit_to_label(q_X), qubit_to_label(q_Z)),
            "sigma_sig": sig,
        })

    # -- Per-particle dump -------------------------------------------------
    print(f"{'particle':<10} {'(p,q)':<8} {'L':<3} "
          f"{'syndrome':<12} {'σ-signature (o0..o6)'}")
    print("-" * 78)
    for r in rows:
        sig_str = ' '.join(str(c) for c in r['sigma_sig'])
        sp = f"({r['synd_pair'][0]},{r['synd_pair'][1]})"
        print(f"{r['name']:<10} ({r['p']:>2},{r['q']:>2})  {r['L']:<3} "
              f"{sp:<12} [{sig_str}]")
    print()

    # -- Group particles by (|p|,|q|) class --------------------------------
    by_pq: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for r in rows:
        by_pq[(abs(r['p']), abs(r['q']))].append(r)

    # Build the (|p|,|q|) → invariant-tuple table using one representative
    # per class (all members of a class have the same walk so same invariants).
    pq_invariants = {}
    for key, entries in by_pq.items():
        r = entries[0]
        pq_invariants[key] = {
            "synd_pair": r['synd_pair'],
            "L": r['L'],
            "sigma_sig": r['sigma_sig'],
        }

    # -- Build cross-sector equivalence groups (Check 5's "twins") ---------
    by_syndrome: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for key, inv in pq_invariants.items():
        by_syndrome[inv['synd_pair']].append(key)

    print("=" * 78)
    print("Check 5 cross-sector twins (syndrome-degenerate (|p|,|q|) classes)")
    print("=" * 78)
    print()
    twin_groups = [grp for grp in by_syndrome.values() if len(grp) > 1]
    for grp in twin_groups:
        sp = pq_invariants[grp[0]]['synd_pair']
        names_per_class = []
        for k in grp:
            names = ','.join(e['name'] for e in by_pq[k])
            names_per_class.append(f"{k}={{{names}}}")
        print(f"  syndrome {sp}: " + ' / '.join(names_per_class))
    print(f"\n  {len(twin_groups)} cross-sector twin groups, "
          f"covering {sum(len(g) for g in twin_groups)} (|p|,|q|) classes")
    print()

    # -- For each twin group, test which augmentations disambiguate --------
    print("=" * 78)
    print("Disambiguation test per cross-sector twin group")
    print("=" * 78)
    print()
    print("For each augmentation candidate, count distinct values across")
    print("the twin (|p|,|q|) classes. A value equal to group size means")
    print("the augmentation disambiguates them; less than group size means")
    print("at least two (|p|,|q|) classes share that augmentation value too.")
    print()

    candidates = [
        ("L",              lambda inv: inv['L']),
        ("σ-sig",          lambda inv: inv['sigma_sig']),
        ("L + σ-sig",      lambda inv: (inv['L'], inv['sigma_sig'])),
    ]

    full_disambig_count = {name: 0 for name, _ in candidates}
    for grp in twin_groups:
        sp = pq_invariants[grp[0]]['synd_pair']
        print(f"  Twin {sp}: classes {grp} (size {len(grp)})")
        for name, extract in candidates:
            vals = [extract(pq_invariants[k]) for k in grp]
            n_distinct = len(set(vals))
            ok = "✓" if n_distinct == len(grp) else "✗"
            if n_distinct == len(grp):
                full_disambig_count[name] += 1
            print(f"    {name:<14} → {n_distinct}/{len(grp)} distinct  {ok}")
        print()

    n_twins = len(twin_groups)
    print(f"\n  SUMMARY across {n_twins} twin groups:")
    for name, _ in candidates:
        cnt = full_disambig_count[name]
        print(f"    syndrome + {name:<14} disambiguates {cnt}/{n_twins} groups")
    print()

    # -- Test global injectivity: does (syndrome, L, σ-sig) distinguish all
    #    16 (|p|,|q|) classes?
    print("=" * 78)
    print("Global injectivity test on all 16 (|p|,|q|) classes")
    print("=" * 78)
    print()
    invariant_groups = [
        ("synd",                   lambda i: i['synd_pair']),
        ("synd + L",               lambda i: (i['synd_pair'], i['L'])),
        ("synd + σ-sig",           lambda i: (i['synd_pair'], i['sigma_sig'])),
        ("synd + L + σ-sig",       lambda i: (i['synd_pair'], i['L'], i['sigma_sig'])),
        ("L + σ-sig only",         lambda i: (i['L'], i['sigma_sig'])),
        ("σ-sig only",             lambda i: i['sigma_sig']),
        ("L only",                 lambda i: i['L']),
    ]
    n_classes = len(pq_invariants)
    for name, extract in invariant_groups:
        values = [extract(inv) for inv in pq_invariants.values()]
        n_distinct = len(set(values))
        verdict = "INJECTIVE ✓" if n_distinct == n_classes else f"loses {n_classes - n_distinct} bit(s)"
        print(f"  {name:<24}  {n_distinct}/{n_classes} distinct   {verdict}")
    print()

    # -- For each non-injective invariant set, show the collisions --------
    print("=" * 78)
    print("Collision detail for non-injective invariant sets")
    print("=" * 78)
    print()
    for name, extract in invariant_groups:
        coll = defaultdict(list)
        for k, inv in pq_invariants.items():
            coll[extract(inv)].append(k)
        collisions = {v: ks for v, ks in coll.items() if len(ks) > 1}
        if not collisions:
            continue
        print(f"  {name}:")
        for v, ks in collisions.items():
            print(f"    {v} → classes {ks}")
        print()


if __name__ == "__main__":
    main()
