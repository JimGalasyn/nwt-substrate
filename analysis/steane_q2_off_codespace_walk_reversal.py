"""Steane-code reframe — Q2: characterize off-codespace walk reversal.

Check 3 established that for the 6 in-centralizer compendium walks
(nucleons + D⁰), walk-reversal acts as H⊗7 EXACTLY at the Pauli-bit
level. For the 19 out-of-centralizer walks, walk-reversal ≠ H⊗7 even
bit-level.

Question (Jim, 2026-05-21 morning): is the off-codespace walk-reversal
*some* closed-form modification of H⊗7? Specifically, is the discrepancy
  Δ = P_reverse - H⊗7(P_forward)        (Pauli-word XOR)
always a stabilizer element, or always a stabilizer × syndrome-dependent
correction?

This script:
  1. Computes (P_forward, P_reverse, H⊗7(P_forward)) for every compendium walk.
  2. Computes the Pauli-word XOR Δ_x = P_rev.x ⊕ (H⊗7 P_fwd).x = P_rev.x ⊕ P_fwd.z
     and Δ_z = P_rev.z ⊕ (H⊗7 P_fwd).z = P_rev.z ⊕ P_fwd.x.
  3. Tests four hypotheses about Δ:
     (a) Δ is identity (CPT = H⊗7 holds exactly) — already known to hold on 6/25.
     (b) Δ_x and Δ_z each lie in the X-stabilizer / Z-stabilizer subspace
         (CPT = H⊗7 mod stabilizer).
     (c) Δ is determined by the forward Pauli's *syndrome* (function of synd_X, synd_Z).
     (d) Δ is determined by σ-orbit signature mismatch under reversal.
  4. Looks for a closed-form Δ(syndrome) rule by tabulating per-syndrome
     Δ vectors.

Honest framing: if Δ has a simple syndrome-dependent description, then
"substrate-CPT" = H⊗7 ∘ Δ(syndrome) and H⊗7 is fundamentally CPT modulo
the substrate's syndrome correction. If Δ has no clean structure, then
walk-reversal is fundamentally CPT and H⊗7 is just a special-case image
on the codespace.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q2_off_codespace_walk_reversal.py
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.particles.compendium import COMPENDIUM


VERTEX_TO_FANO = {
    0: (1, 1, 1), 1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1),
    4: (0, 1, 1), 5: (1, 0, 1), 6: (1, 1, 0),
}
VERTEX_LABEL = {0: 'P', 1: 'E1', 2: 'E2', 3: 'E3', 4: 'F1', 5: 'F2', 6: 'F3'}


def std_position(point):
    return 4 * point[0] + 2 * point[1] + 1 * point[2]


VERTEX_TO_QUBIT = {v: std_position(VERTEX_TO_FANO[v]) - 1 for v in range(7)}
QUBIT_TO_VERTEX = {q: v for v, q in VERTEX_TO_QUBIT.items()}

QR_DIRECTIONS = {1, 2, 4}
NR_DIRECTIONS = {3, 5, 6}


def hamming_parity_check():
    H = np.zeros((3, 7), dtype=int)
    for k in range(3):
        for q in range(7):
            H[k, q] = ((q + 1) >> k) & 1
    return H


HAMMING = hamming_parity_check()
X_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]
Z_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]


def walk_direction_pauli(walk):
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


def x_syndrome(z_part):
    return tuple(sum(z_part[i] * g[i] for i in range(7)) % 2 for g in X_STAB_GENS)


def z_syndrome(x_part):
    return tuple(sum(x_part[i] * g[i] for i in range(7)) % 2 for g in Z_STAB_GENS)


def syndrome_to_qubit(syn):
    idx = syn[0] + 2 * syn[1] + 4 * syn[2]
    return idx - 1 if idx > 0 else None


def qubit_label(q):
    return "I" if q is None else VERTEX_LABEL[QUBIT_TO_VERTEX[q]]


def x_stab_subspace():
    """Enumerate all 8 elements of the X-stabilizer subspace (Z₂-span of 3 generators)."""
    sub = []
    for mask in range(8):
        v = np.zeros(7, dtype=int)
        for k in range(3):
            if (mask >> k) & 1:
                v ^= np.array(X_STAB_GENS[k])
        sub.append(tuple(int(x) for x in v))
    return sub


X_STAB_SUBSPACE = frozenset(x_stab_subspace())
Z_STAB_SUBSPACE = X_STAB_SUBSPACE  # CSS symmetric


def xor(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def main():
    print("=" * 78)
    print("Steane Q2 — off-codespace walk-reversal structure")
    print("=" * 78)
    print()
    print("For each walk W:")
    print("  P_fwd  = direction-sensitive Pauli of W            (x_fwd, z_fwd)")
    print("  P_rev  = direction-sensitive Pauli of reverse W    (x_rev, z_rev)")
    print("  H⊗7 P_fwd = (z_fwd, x_fwd)        (transversal Hadamard swap)")
    print("  Δ_x = x_rev ⊕ z_fwd               (X-part discrepancy)")
    print("  Δ_z = z_rev ⊕ x_fwd               (Z-part discrepancy)")
    print()
    print("CPT = H⊗7 holds bit-level iff Δ_x = Δ_z = 0.")
    print("CPT = H⊗7 holds mod stabilizer iff Δ_x ∈ X-stab subspace AND")
    print("Δ_z ∈ Z-stab subspace.")
    print()

    walks = bfs_shortest_walks(max_length=25)

    rows = []
    seen_pq = set()
    for entry in COMPENDIUM:
        key = (abs(entry['p']), abs(entry['q']))
        if key in seen_pq or key not in walks:
            continue
        seen_pq.add(key)
        walk = walks[key]
        rwalk = list(reversed(walk))
        x_fwd, z_fwd = walk_direction_pauli(walk)
        x_rev, z_rev = walk_direction_pauli(rwalk)
        delta_x = xor(x_rev, z_fwd)
        delta_z = xor(z_rev, x_fwd)
        synd_X = x_syndrome(z_fwd)
        synd_Z = z_syndrome(x_fwd)
        rows.append({
            'name': entry['name'], 'p': key[0], 'q': key[1],
            'x_fwd': x_fwd, 'z_fwd': z_fwd,
            'x_rev': x_rev, 'z_rev': z_rev,
            'delta_x': delta_x, 'delta_z': delta_z,
            'synd_pair': (qubit_label(syndrome_to_qubit(synd_X)),
                          qubit_label(syndrome_to_qubit(synd_Z))),
        })

    # -- Per-walk table ----------------------------------------------------
    print(f"{'(p,q)':<8} {'name':<8} {'syndrome':<12} "
          f"{'Δ_x':<10} {'Δ_z':<10} "
          f"{'Δ_x∈⟨X⟩':<10} {'Δ_z∈⟨Z⟩':<10} {'CPT=H⊗7?'}")
    print("-" * 90)
    n_exact = n_modstab = 0
    for r in rows:
        dx_in = r['delta_x'] in X_STAB_SUBSPACE
        dz_in = r['delta_z'] in Z_STAB_SUBSPACE
        if r['delta_x'] == (0,) * 7 and r['delta_z'] == (0,) * 7:
            verdict = "EXACT"
            n_exact += 1
            n_modstab += 1
        elif dx_in and dz_in:
            verdict = "mod stab"
            n_modstab += 1
        else:
            verdict = "✗"
        dx_str = ''.join(str(b) for b in r['delta_x'])
        dz_str = ''.join(str(b) for b in r['delta_z'])
        sp = f"({r['synd_pair'][0]},{r['synd_pair'][1]})"
        print(f"({r['p']:>2},{r['q']:>2})  {r['name']:<8} {sp:<12} "
              f"{dx_str:<10} {dz_str:<10} "
              f"{str(dx_in):<10} {str(dz_in):<10} {verdict}")
    print()
    print(f"  EXACT CPT=H⊗7:       {n_exact}/{len(rows)}")
    print(f"  CPT=H⊗7 mod stab:    {n_modstab}/{len(rows)}")
    print(f"  Genuinely off:       {len(rows) - n_modstab}/{len(rows)}")
    print()

    # -- Hypothesis: Δ depends only on syndrome ----------------------------
    print("=" * 78)
    print("HYPOTHESIS A: Δ = function of syndrome pair only")
    print("=" * 78)
    print()
    by_syn = defaultdict(list)
    for r in rows:
        by_syn[r['synd_pair']].append(r)
    consistent = True
    for sp, group in sorted(by_syn.items()):
        dxs = set(r['delta_x'] for r in group)
        dzs = set(r['delta_z'] for r in group)
        ok = len(dxs) == 1 and len(dzs) == 1
        marker = "✓" if ok else "✗"
        if not ok:
            consistent = False
        names = ','.join(r['name'] for r in group)
        sp_str = f"({sp[0]},{sp[1]})"
        print(f"  syndrome {sp_str:<14} [{names}]: "
              f"{len(dxs)} distinct Δ_x, {len(dzs)} distinct Δ_z  {marker}")
    print()
    print(f"  HYPOTHESIS A {'CONFIRMED' if consistent else 'FALSIFIED'}: "
          f"Δ is{'' if consistent else ' NOT'} a function of syndrome alone")
    print()

    # -- Hypothesis: Δ_x = z_fwd XOR something_structural ------------------
    print("=" * 78)
    print("HYPOTHESIS B: Δ_x relates to z_fwd via a syndrome-keyed transform")
    print("=" * 78)
    print()
    print("Δ_x = z_rev ⊕ x_fwd? (Δ_z = x_rev ⊕ z_fwd was the by-definition.")
    print("Now test the cross-pairing.)")
    print()
    for r in rows:
        cross_x = xor(r['x_rev'], r['x_fwd'])
        cross_z = xor(r['z_rev'], r['z_fwd'])
        cx_str = ''.join(str(b) for b in cross_x)
        cz_str = ''.join(str(b) for b in cross_z)
        sp = f"({r['synd_pair'][0]},{r['synd_pair'][1]})"
        print(f"  ({r['p']:>2},{r['q']:>2}) {r['name']:<8} {sp:<14} "
              f"x_rev⊕x_fwd={cx_str}  z_rev⊕z_fwd={cz_str}")
    print()

    # -- Hypothesis: walk reversal acts as H⊗7 mod additional stabilizer
    # -- equivalence — count classes -----------------------------------
    print("=" * 78)
    print("HYPOTHESIS C: how many walks satisfy CPT=H⊗7 mod stabilizer?")
    print("=" * 78)
    print()
    print(f"  CPT = H⊗7 mod stabilizer requires Δ_x ∈ X-stab subspace")
    print(f"  and Δ_z ∈ Z-stab subspace (8 elements each).")
    print()
    print(f"  Result: {n_modstab}/{len(rows)} walks have CPT=H⊗7 mod stab.")
    print()
    if n_modstab > n_exact:
        added = n_modstab - n_exact
        print(f"  → {added} additional walks beyond the centralizer satisfy")
        print(f"    CPT=H⊗7 modulo stabilizer correction. These are walks")
        print(f"    where Δ is a stabilizer element (nontrivial syndrome,")
        print(f"    same logical-coset image, walk-reversal up to gauge).")
    print()

    # -- Hypothesis: the off-codespace walks split by Δ-coset --------------
    print("=" * 78)
    print("HYPOTHESIS D: Δ-cosets of out-of-centralizer walks (substrate structure)")
    print("=" * 78)
    print()
    delta_classes = defaultdict(list)
    for r in rows:
        delta_classes[(r['delta_x'], r['delta_z'])].append(r)
    print(f"  {len(delta_classes)} distinct Δ-pairs across {len(rows)} walks")
    print()
    for (dx, dz), group in sorted(delta_classes.items(),
                                    key=lambda kv: -len(kv[1])):
        names = ','.join(f"{r['name']}({r['p']},{r['q']})" for r in group)
        dx_str = ''.join(str(b) for b in dx)
        dz_str = ''.join(str(b) for b in dz)
        dx_in = dx in X_STAB_SUBSPACE
        dz_in = dz in Z_STAB_SUBSPACE
        kind = ("identity" if (dx == (0,) * 7 and dz == (0,) * 7)
                else "stab" if (dx_in and dz_in)
                else "non-stab")
        print(f"  Δ=({dx_str},{dz_str}) [{kind}] — {len(group)} walks: {names}")
    print()


if __name__ == "__main__":
    main()
