"""Steane-code reframe — Check 2: walk → logical Pauli coset (visit-parity).

For each compendium walk W:
  1. Compute visit-count parity on each K_7 vertex.
  2. Form three candidate Pauli words from that parity:
       X_W = X on each odd-parity vertex
       Z_W = Z on each odd-parity vertex
       Y_W = Y on each odd-parity vertex   (= X_W · Z_W up to phase)
  3. Check centralizer condition (commutes with all 6 Steane stabilizers).
  4. If in centralizer, reduce modulo stabilizer to one of {I_L, X_L, Y_L, Z_L}.
  5. Tabulate by Paper-11 carrier-knot sector n_q ∈ {0, 2, 3, 5}.

Limitation flagged in advance (the "wrinkle"): visit-parity is the same for
matter and antimatter walks of the same particle (the reverse of a closed
walk visits the same vertex multiset).  So this attempt cannot directly
distinguish e^- from e^+.  But it might still cleanly partition the four
carrier-knot sectors {unknot, Hopf, trefoil, cinquefoil}.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_check2_walk_to_logical_pauli.py
"""
from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.particles.compendium import COMPENDIUM


# ---------------------------------------------------------------------------
# Fano-point labeling (matches Check 1 / 1b)
#   K_7 vertex k (k=0..6) is mapped to a Fano point.
#   Visualization convention (Paper 21 §4 figures):
#       0 = P,  1 = E_1,  2 = E_2,  3 = E_3,  4 = F_1,  5 = F_2,  6 = F_3
# ---------------------------------------------------------------------------

VERTEX_TO_FANO = {
    0: (1, 1, 1),    # P (polar)
    1: (1, 0, 0),    # E_1
    2: (0, 1, 0),    # E_2
    3: (0, 0, 1),    # E_3
    4: (0, 1, 1),    # F_1
    5: (1, 0, 1),    # F_2
    6: (1, 1, 0),    # F_3
}
VERTEX_LABEL = {0: 'P', 1: 'E1', 2: 'E2', 3: 'E3', 4: 'F1', 5: 'F2', 6: 'F3'}


def std_position(point: tuple) -> int:
    """Map Fano point (b2, b1, b0) → integer 1..7 (standard Steane indexing)."""
    return 4 * point[0] + 2 * point[1] + 1 * point[2]


# Reindex vertices 0..6 by std_position - 1, so qubit indices match Hamming positions
VERTEX_TO_QUBIT = {v: std_position(VERTEX_TO_FANO[v]) - 1 for v in range(7)}
QUBIT_TO_VERTEX = {q: v for v, q in VERTEX_TO_QUBIT.items()}


# ---------------------------------------------------------------------------
# Steane [[7,1,3]] stabilizers and logical operators (symplectic rep on F_2^14)
# ---------------------------------------------------------------------------

def hamming_parity_check():
    """3 × 7 Hamming (7,4) parity-check rows; entry i,j = bit i of (j+1)."""
    H = np.zeros((3, 7), dtype=int)
    for k in range(3):
        for q in range(7):
            pos = q + 1
            H[k, q] = (pos >> k) & 1
    return H


HAMMING = hamming_parity_check()
# 3 X-stabilizer generators (supports = rows of H)
X_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]
Z_STAB_GENS = [tuple(HAMMING[k]) for k in range(3)]

# Logical operators: weight-3 codewords of the (7,4) Hamming code.
# We pick X_L_support = qubits {0, 1, 2} = positions {1, 2, 3} = Fano line.
# (001 + 010 + 011 = 000 ✓.)  Z_L on the same support.
X_L_SUPPORT = (1, 1, 1, 0, 0, 0, 0)  # qubits 0, 1, 2
Z_L_SUPPORT = (1, 1, 1, 0, 0, 0, 0)


def in_x_stab_subspace(v: tuple) -> bool:
    """Is v ∈ (F_2)^7 in the 3-dim X-stabilizer subspace?"""
    target = np.array(v, dtype=int)
    for mask in range(8):
        x = np.zeros(7, dtype=int)
        for k in range(3):
            if (mask >> k) & 1: x ^= np.array(X_STAB_GENS[k])
        if np.array_equal(x, target): return True
    return False


def in_z_stab_subspace(v: tuple) -> bool:
    target = np.array(v, dtype=int)
    for mask in range(8):
        z = np.zeros(7, dtype=int)
        for k in range(3):
            if (mask >> k) & 1: z ^= np.array(Z_STAB_GENS[k])
        if np.array_equal(z, target): return True
    return False


def x_only_syndrome(v: tuple) -> tuple:
    """For an X-only Pauli with support v, return Z-stabilizer syndrome."""
    s = []
    for g in Z_STAB_GENS:
        s.append(sum(v[i]*g[i] for i in range(7)) % 2)
    return tuple(s)


# ---------------------------------------------------------------------------
# Symplectic Pauli helpers
# ---------------------------------------------------------------------------

def pauli_commute(p1, p2) -> bool:
    """Commutation in symplectic rep: [(a,b), (c,d)] = 0 iff a·d + b·c even."""
    a, b = p1
    c, d = p2
    s = sum(a[i]*d[i] for i in range(7)) + sum(b[i]*c[i] for i in range(7))
    return s % 2 == 0


def in_centralizer(pauli) -> bool:
    """Pauli commutes with all 6 Steane stabilizer generators?"""
    for s in X_STAB_GENS:
        if not pauli_commute(pauli, (s, (0,)*7)):
            return False
    for s in Z_STAB_GENS:
        if not pauli_commute(pauli, ((0,)*7, s)):
            return False
    return True


# Stabilizer group elements (8 X-stabs × 8 Z-stabs = 64 elements)
def all_stab_elements():
    out = []
    for x_mask in range(8):
        x = np.zeros(7, dtype=int)
        for k in range(3):
            if (x_mask >> k) & 1: x ^= np.array(X_STAB_GENS[k])
        for z_mask in range(8):
            z = np.zeros(7, dtype=int)
            for k in range(3):
                if (z_mask >> k) & 1: z ^= np.array(Z_STAB_GENS[k])
            out.append((tuple(x), tuple(z)))
    return out


STAB_ELEMENTS = all_stab_elements()
STAB_SET = set(STAB_ELEMENTS)


def reduce_modulo_stabilizer(pauli) -> str:
    """Identify the logical-Pauli coset of `pauli` ∈ centralizer.

    Returns one of "I_L", "X_L", "Z_L", "Y_L".

    Cosets determined by 2 bits: (X-part mod X-stab subspace, Z-part mod Z-stab).
    Each bit is 0 if the part is in the stabilizer subspace, 1 if it is in
    the logical-X/Z coset.
    """
    x_bit = 0 if in_x_stab_subspace(pauli[0]) else 1
    z_bit = 0 if in_z_stab_subspace(pauli[1]) else 1
    return {(0, 0): "I_L", (1, 0): "X_L",
             (0, 1): "Z_L", (1, 1): "Y_L"}[(x_bit, z_bit)]


# ---------------------------------------------------------------------------
# Walk → visit-parity Pauli
# ---------------------------------------------------------------------------

def walk_visit_parity(walk: list[int]) -> tuple:
    """Per-K_7-vertex parity of visit count.

    Convention: visit list = walk[:-1] (all vertices traversed before closing).
    Returns a 7-tuple of 0/1 indexed by *qubit* position (not K_7 vertex).
    """
    parity = [0] * 7
    for v in walk[:-1]:
        q = VERTEX_TO_QUBIT[v]
        parity[q] ^= 1
    return tuple(parity)


def pauli_X(parity: tuple) -> tuple:
    return (parity, (0,)*7)


def pauli_Z(parity: tuple) -> tuple:
    return ((0,)*7, parity)


def pauli_Y(parity: tuple) -> tuple:
    return (parity, parity)


def pauli_fmt(p) -> str:
    return ''.join(str(b) for b in p[0]) + '|' + ''.join(str(b) for b in p[1])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("Steane Check 2 — walk → logical Pauli coset via visit parity")
    print("=" * 78)
    print()
    print("K_7 vertex ↔ Steane qubit mapping:")
    for v in range(7):
        q = VERTEX_TO_QUBIT[v]
        print(f"  vertex {v} ({VERTEX_LABEL[v]}, Fano {VERTEX_TO_FANO[v]}) "
               f"→ qubit {q} (Hamming pos {q+1})")
    print()
    print(f"Logical X_L support = qubits {[i for i, b in enumerate(X_L_SUPPORT) if b]}")
    print(f"Logical Z_L support = qubits {[i for i, b in enumerate(Z_L_SUPPORT) if b]}")
    print()

    walks = bfs_shortest_walks(max_length=25)
    print(f"Loaded {len(walks)} (|p|, |q|) shortest walks")
    print()

    # ---- Per-walk computation -----------------------------------------
    rows = []
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        walk = walks[key]
        parity = walk_visit_parity(walk)

        pX = pauli_X(parity)
        pY = pauli_Y(parity)
        pZ = pauli_Z(parity)

        # Centralizer membership
        cx = in_centralizer(pX)
        cy = in_centralizer(pY)
        cz = in_centralizer(pZ)

        # Reduce mod stabilizer if in centralizer; otherwise compute syndrome
        rx = reduce_modulo_stabilizer(pX) if cx else "—"
        ry = reduce_modulo_stabilizer(pY) if cy else "—"
        rz = reduce_modulo_stabilizer(pZ) if cz else "—"
        synd = x_only_syndrome(parity)  # syndrome from Z-stab anti-commutation
        synd_int = synd[0] + 2*synd[1] + 4*synd[2]  # 0..7 (= corrected qubit + 1)
        corrected_qubit = synd_int - 1 if synd_int > 0 else None

        rows.append({**entry, "walk": walk, "parity": parity,
                      "cx": cx, "cy": cy, "cz": cz,
                      "rx": rx, "ry": ry, "rz": rz,
                      "syndrome": synd, "corrected_qubit": corrected_qubit})

    print(f"{'particle':<10} {'n_q':<3} {'(p,q)':<8} {'L':<3} {'parity':<10} "
           f"{'X-cent':<7} {'coset':<6} {'syndrome':<10} {'corrected qubit'}")
    print("-" * 80)
    for r in rows:
        par_str = ''.join(str(b) for b in r['parity'])
        syn_str = ''.join(str(b) for b in r['syndrome'])
        cq = r['corrected_qubit']
        cq_str = (f"q{cq} ({VERTEX_LABEL[QUBIT_TO_VERTEX[cq]]})"
                   if cq is not None else "—")
        print(f"{r['name']:<10} {r['n_q']:<3} "
               f"({r['p']:>2},{r['q']:>2})  {len(r['walk'])-1:<3} {par_str:<10} "
               f"{str(r['cx']):<7} {r['rx']:<6} {syn_str:<10} {cq_str}")
    print()

    # ---- Sector consistency analysis -----------------------------------
    print("=" * 78)
    print("SECTOR CONSISTENCY — does each n_q sector map to a unique coset?")
    print("=" * 78)
    print()

    sector_name = {0: 'leptons (unknot)', 2: 'mesons (Hopf)',
                    3: 'hyperons (trefoil)', 5: 'nucleons (cinquefoil)'}
    for variant, key_cent, key_cos in [("X", "cx", "rx"), ("Y", "cy", "ry"),
                                         ("Z", "cz", "rz")]:
        print(f"  --- {variant}-variant ---")
        by_nq = defaultdict(list)
        for r in rows:
            by_nq[r["n_q"]].append((r["name"], r[key_cent], r[key_cos]))
        for nq in sorted(by_nq.keys()):
            sec = sector_name.get(nq, f'n_q={nq}')
            entries = by_nq[nq]
            in_cent_count = sum(1 for _, c, _ in entries if c)
            cosets = set(co for _, c, co in entries if c)
            outside = [n for n, c, _ in entries if not c]
            print(f"    n_q={nq} ({sec}): {len(entries)} walks, "
                   f"{in_cent_count} in centralizer; cosets seen = {cosets}")
            if outside:
                print(f"      outside centralizer: {outside}")
        print()

    # ---- Syndrome distribution by sector ------------------------------
    print("=" * 78)
    print("SYNDROME DISTRIBUTION (X-only Pauli, outside-centralizer walks)")
    print("=" * 78)
    print()
    print("Each X-only Pauli that's not in the centralizer has a 3-bit")
    print("Hamming syndrome from anti-commutation with Z-stabilizers.")
    print("Syndrome (s_1, s_2, s_3) → integer s_1 + 2 s_2 + 4 s_3 = qubit index + 1")
    print("This identifies the qubit (= Fano point) that the Hamming code")
    print("would 'correct' to bring the Pauli into the centralizer.")
    print()
    print(f"  {'sector':<10} {'walks':<6} {'syndromes seen':<60}")
    print("  " + "-" * 75)
    by_nq_synd = defaultdict(list)
    for r in rows:
        synd_int = r['syndrome'][0] + 2*r['syndrome'][1] + 4*r['syndrome'][2]
        if synd_int > 0:
            q = synd_int - 1
            v = QUBIT_TO_VERTEX[q]
            by_nq_synd[r['n_q']].append(f"{r['name']}→{VERTEX_LABEL[v]}")
    for nq in sorted(by_nq_synd.keys()):
        sec = sector_name.get(nq, f'n_q={nq}')
        entries = by_nq_synd[nq]
        print(f"  {sec:<10} {len(entries):<6} {', '.join(entries)}")
    print()
    # Tabulate which Fano point each syndrome identifies, per sector
    print("  Syndrome-identified Fano point tally by sector:")
    for nq in sorted(by_nq_synd.keys()):
        sec = sector_name.get(nq, f'n_q={nq}')
        fano_pts = [e.split('→')[1] for e in by_nq_synd[nq]]
        from collections import Counter
        c = Counter(fano_pts)
        print(f"    n_q={nq} ({sec}): {dict(c)}")
    print()

    # ---- Inspect particularly: matter walks vs reverse walks ----------
    print("=" * 78)
    print("Matter/antimatter check (proton vs antiproton, etc.)")
    print("=" * 78)
    print()
    # Reverse compendium walks and re-run
    print("  Visit parity is invariant under walk reversal (by construction).")
    print("  So matter/antimatter walks have IDENTICAL Pauli words under this map.")
    print("  → visit-parity alone cannot distinguish matter from antimatter.")
    print()

    # ---- Headline ------------------------------------------------------
    print("=" * 78)
    print("HEADLINE")
    print("=" * 78)
    print()
    n_total = len(rows)
    for variant, key_cent in [("X", "cx"), ("Y", "cy"), ("Z", "cz")]:
        in_cent = sum(1 for r in rows if r[key_cent])
        print(f"  {variant}-variant: {in_cent}/{n_total} walks in centralizer")
    print()


if __name__ == "__main__":
    main()
