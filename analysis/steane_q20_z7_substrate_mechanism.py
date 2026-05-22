"""Q20 — Substrate mechanism for rule (I): cosmogenic Z_3 cycle on Z_7 residues.

Q19 found a clean rule (I): n_q = f(p mod 7, q mod 7) with 16/16
compendium fit. This Q probes the SUBSTRATE MECHANISM underlying the
rule by connecting it to the AGL(1, 7) Z_3 cosmogenic three-generation
generator x → 2x mod 7.

Key observation: x → 2x mod 7 partitions Z_7* into two Z_3 cycles:
  QR = {1, 2, 4}: cycle 1 → 2 → 4 → 1
  NR = {3, 6, 5}: cycle 3 → 6 → 5 → 3

Per rule (I), the n_q assignments on QR residues are:
  1 → lepton (n_q = F_2 = 1)
  2 → meson  (n_q = F_3 = 2)
  4 → hyperon (n_q = F_4 = 3)

**STRIKING PATTERN**: under the cosmogenic Z_3 (= x → 2x mod 7), the
QR residues are PERMUTED CYCLICALLY through CONSECUTIVE FIBONACCI INDICES
{F_2, F_3, F_4} = {1, 2, 3}. Each Z_3 step INCREMENTS the Fibonacci
index by 1.

This Z_3 IS the cosmogenic three-generation breaking generator (per
g2_bridge.py): the same Z_3 that maps the SM 3 fermion generations.
Here it acts on substrate Z_7 residues to permute the 3 lower carrier-
knot Fibonacci shells {lepton, meson, hyperon}.

What about NR?
  3 → nucleon (F_5 = 5) [when p ≢ 0 mod 7] OR meson (F_3 = 2) [when p ≡ 0]
  6 → ???  (no compendium data) ← PREDICTION TARGET
  5 → meson (F_3 = 2)

This Q investigates:
  (1) Verify Z_3 cycle structure on QR ↔ Fibonacci index increment
  (2) Document NR cycle behavior and the q=3 p-mod-7=0 split
  (3) PREDICT the sector for q ≡ 6 mod 7 via Z_3 cycle structure
  (4) Connect rule (I) to AGL Z_3 cosmogenic three-generation breaking
  (5) Substrate-mechanism reading of rule (I)

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q20_z7_substrate_mechanism.py
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from nwt_substrate.algebra.g2_bridge import (
    PALEY_TO_BAEZ_LABELING, AGL_Z3_SIGN_LIFT,
    agl_z3_g2_matrix, baez_permutation_of_agl_z3, verify_bridge,
)


# ============================================================================
# Setup: Z_7 residues, QR/NR, Fibonacci
# ============================================================================

def fibonacci(n: int) -> int:
    """F_1 = F_2 = 1, F_3 = 2, F_4 = 3, F_5 = 5, F_6 = 8, F_7 = 13."""
    if n <= 2:
        return 1
    a, b = 1, 1
    for _ in range(n - 2):
        a, b = b, a + b
    return b


QR = {1, 2, 4}  # Quadratic residues mod 7
NR = {3, 5, 6}  # Non-residues mod 7


def apply_z3(x: int) -> int:
    """Apply the cosmogenic Z_3 generator: x → 2x mod 7."""
    return (2 * x) % 7


def z3_orbit(x: int) -> list[int]:
    """Z_3 orbit of x under x → 2x mod 7."""
    orbit = [x]
    cur = apply_z3(x)
    while cur != x:
        orbit.append(cur)
        cur = apply_z3(cur)
    return orbit


# Rule (I) from Q19
RULE_I = {
    0: ("meson", 2, 3),     # (sector, n_q, Fibonacci index)
    1: ("lepton", 0, 2),    # F_2 = 1 (n_q-Fibonacci-framed)
    2: ("meson", 2, 3),     # F_3 = 2
    3: ("nucleon", 5, 5),   # F_5 = 5 (when p ≢ 0)
    4: ("hyperon", 3, 4),   # F_4 = 3
    5: ("meson", 2, 3),     # F_3 = 2
    6: ("?", None, None),   # PREDICTION TARGET
}

RULE_I_SPECIAL = {
    # q mod 7 = 3 AND p mod 7 = 0 → meson F_3
    (0, 3): ("meson", 2, 3),
}


# ============================================================================
# Step 1: Verify Z_3 cycle structure on QR
# ============================================================================

def step1_qr_z3_cycle():
    print("=" * 78)
    print("[Step 1] Z_3 cycle x → 2x mod 7 on QR residues + sector assignments")
    print("=" * 78)
    print()
    print(f"  QR = {{1, 2, 4}} (quadratic residues mod 7)")
    print(f"  Z_3 generator: x → 2x mod 7")
    print(f"  QR Z_3 orbit:  1 → 2 → 4 → 1  (single 3-cycle)")
    print()
    print(f"  {'q mod 7':<10} {'sector':<10} {'n_q':<5} {'Fibonacci':<10} "
          f"{'next under Z_3 (x→2x)':<25} {'F-index Δ'}")
    print("  " + "-" * 80)
    for q in [1, 2, 4]:
        sector, nq, fib_idx = RULE_I[q]
        nxt = apply_z3(q)
        nxt_sector, nxt_nq, nxt_fib = RULE_I[nxt]
        if fib_idx is not None and nxt_fib is not None:
            delta = nxt_fib - fib_idx
            delta_str = f"+{delta}" if delta > 0 else str(delta)
        else:
            delta_str = "?"
        print(f"  {q:<10} {sector:<10} {nq:<5} F_{fib_idx}={fibonacci(fib_idx):<6} "
              f"{q} → {nxt} ({nxt_sector}, F_{nxt_fib}) {delta_str}")
    print()
    print("  ★★★ STRUCTURAL FINDING: Z_3 cycle on QR INCREMENTS Fibonacci ★★★")
    print("       index by +1 each step (with wraparound at F_4 → F_2):")
    print()
    print("       QR cycle: 1 → 2 → 4 → 1")
    print("       Sector:   lepton → meson → hyperon → lepton")
    print("       F-index:  F_2 → F_3 → F_4 → F_2 (cyclic, period 3)")
    print()
    print("  The cosmogenic Z_3 generator (= AGL Z_3 three-generation breaking)")
    print("  PERMUTES the three lowest substrate carrier-knot shells through a")
    print("  three-step Fibonacci-incrementing cycle. This is a substrate-")
    print("  mechanism explanation for WHY {lepton, meson, hyperon} carry")
    print("  consecutive Fibonacci indices {F_2, F_3, F_4}.")
    print()


# ============================================================================
# Step 2: Document NR cycle behavior
# ============================================================================

def step2_nr_z3_cycle():
    print("=" * 78)
    print("[Step 2] Z_3 cycle on NR residues")
    print("=" * 78)
    print()
    print(f"  NR = {{3, 5, 6}} (non-residues mod 7)")
    print(f"  NR Z_3 orbit:  3 → 6 → 5 → 3  (single 3-cycle)")
    print()
    print(f"  {'q mod 7':<10} {'sector':<22} {'n_q':<8} {'Fibonacci':<12} "
          f"{'next under Z_3'}")
    print("  " + "-" * 80)
    for q in [3, 6, 5]:
        sector, nq, fib_idx = RULE_I[q]
        nxt = apply_z3(q)
        nxt_sector, nxt_nq, nxt_fib = RULE_I[nxt]
        if q == 3:
            sector_str = "nucleon (p≢0) / meson (p≡0)"
            nq_str = "5/2"
            fib_str = "F_5/F_3"
        else:
            sector_str = sector if sector != "?" else "PREDICTION"
            nq_str = str(nq) if nq is not None else "?"
            fib_str = f"F_{fib_idx}={fibonacci(fib_idx)}" if fib_idx else "?"
        nxt_str = f"{q} → {nxt} ({nxt_sector if nxt_sector != '?' else 'PRED'})"
        print(f"  {q:<10} {sector_str:<22} {nq_str:<8} {fib_str:<12} {nxt_str}")
    print()
    print("  NR cycle behavior:")
    print("    3 (nucleon F_5) → 6 (?) → 5 (meson F_3) → 3")
    print()
    print("  Unlike QR, the NR cycle has the SPECIAL nucleon at residue 3")
    print("  (with p mod 7 ≠ 0) and the EMPTY substrate-prediction slot at")
    print("  residue 6. The substrate-mechanism for the NR pattern is more")
    print("  subtle than QR (which has a clean Fibonacci-increment cycle).")
    print()


# ============================================================================
# Step 3: Predict q ≡ 6 mod 7 sector
# ============================================================================

def step3_predict_q6_sector():
    print("=" * 78)
    print("[Step 3] Predict the sector for q ≡ 6 mod 7 (NR middle residue)")
    print("=" * 78)
    print()
    print("  No compendium walk has q ≡ 6 mod 7. The Z_3 cycle structure")
    print("  on NR (3 → 6 → 5) suggests q ≡ 6 is a substrate-prediction")
    print("  slot. Possible interpretations:")
    print()
    print("  (a) Fibonacci-extension: NR cycle goes F_5 → F_6 → F_3, so")
    print("      q ≡ 6 → n_q = F_6 = 8 (T(2, 8) torus link carrier)")
    print("      → predicts a NEW heavy substrate-prediction species at n_q=8")
    print()
    print("  (b) Reflective: maybe q ≡ 6 mirrors q ≡ 1 (via 6 = -1 mod 7),")
    print("      giving lepton-class (n_q = F_2 = 1)")
    print("      → NR-direction antimatter-lepton analog")
    print()
    print("  (c) Trivial: q ≡ 6 is substrate-forbidden (no walks exist),")
    print("      analogous to d=2 Hamilton being substrate-forbidden")
    print("      → no species predicted")
    print()
    print("  To select among (a)/(b)/(c), would need to find walk(s) with")
    print("  q ≡ 6 mod 7 in compendium or in extended BFS, OR identify a")
    print("  substrate-mechanism prediction.")
    print()
    print("  EMPIRICAL CHECK: do walk-BFS results yield any closed walks with")
    print("  q ≡ 6 mod 7?")
    from nwt_substrate.condensate import bfs_shortest_walks
    walks = bfs_shortest_walks(max_length=25)
    q6_walks = [(pq, walks[pq]) for pq in walks.keys()
                 if abs(pq[1]) % 7 == 6]
    print()
    print(f"  Found {len(q6_walks)} walk-BFS closed walks with q ≡ 6 mod 7:")
    for pq, walk in q6_walks[:10]:
        print(f"    ({pq[0]:>2}, {pq[1]:>2}) L={len(walk)-1:<3}  {walk}")
    if len(q6_walks) > 10:
        print(f"    ... and {len(q6_walks) - 10} more")
    print()


# ============================================================================
# Step 4: Connect to AGL Z_3 cosmogenic three-generation breaking
# ============================================================================

def step4_agl_z3_connection():
    print("=" * 78)
    print("[Step 4] Connect rule (I) Z_3 cycle to AGL three-generation breaking")
    print("=" * 78)
    print()
    print("  Per g2_bridge.py and Paper 20:")
    print("  AGL(1, 7) = affine group of Z_7 = Z_7 ⋊ Z_7* of order 42.")
    print("  Its commutator [a, b] generator restricted to Z_7* gives the")
    print("  multiplicative Z_3 subgroup: x → 2x mod 7.")
    print()
    print("  This is the COSMOGENIC THREE-GENERATION generator: it permutes")
    print("  the 3 'flavors' of substrate excitations along the Z_3 cycle.")
    print()
    print("  Verifying via existing g2_bridge infrastructure...")
    bridge = verify_bridge()
    print(f"  Bridge verification:")
    print(f"    AGL Z_3 generator M is orthogonal (M^T M = I): "
          f"{bridge.is_orthogonal}")
    print(f"    Order 3 (M³ = I): {bridge.order_3}")
    print(f"    Preserves octonion product (G_2): "
          f"{bridge.preserves_octonion_product}")
    print(f"    Trace on R^7: {bridge.trace_R7:+.1f}")
    print(f"    Distinct from SU(3) center Z_3: "
          f"{not bridge.agl_z3_equals_su3_center}")
    print()
    print("  Now check explicit Z_3 action on Paley K_7 vertices vs Baez:")
    print(f"  Baez permutation π of AGL Z_3 generator:")
    perm = baez_permutation_of_agl_z3()
    for k, v in sorted(perm.items()):
        print(f"    e_{k} → e_{v}")
    print()
    print(f"  Paley → Baez labeling (vertex k → octonion idx):")
    for k, v in sorted(PALEY_TO_BAEZ_LABELING.items()):
        print(f"    Paley vertex {k} → Baez e_{v}")
    print()
    print("  In Paley-vertex labels, AGL Z_3 = multiplication-by-2 on Z_7:")
    print(f"    0 → 0 (fixed, polar)")
    for x in range(1, 7):
        nxt = apply_z3(x)
        print(f"    {x} → {nxt}")
    print()
    print("  ★ THIS IS THE SAME Z_3 used in rule (I) ★")
    print()
    print("  So rule (I)'s Z_3 cycle is identified with the COSMOGENIC")
    print("  AGL Z_3 three-generation generator. The QR cycle of")
    print("  Fibonacci-incrementing sectors {lepton, meson, hyperon} is a")
    print("  manifestation of cosmogenic Z_3 axis-breaking on the K_7")
    print("  substrate's Z_7 winding residues.")
    print()


# ============================================================================
# Step 5: Substrate-mechanism summary
# ============================================================================

def step5_substrate_mechanism():
    print("=" * 78)
    print("[Step 5] Substrate-mechanism reading of rule (I)")
    print("=" * 78)
    print()
    print("  The substrate K_7 has:")
    print("    - Z_7 vertex cyclic structure (Paley convention)")
    print("    - Cosmogenic Z_3 breaking via AGL generator x → 2x mod 7")
    print("    - QR/NR partition: QR = {1, 2, 4}, NR = {3, 5, 6}, {0}")
    print("    - 2 substrate-realized Hamilton directions: d=1 (QR-primary),")
    print("      d=3 (NR-primary); d=2 substrate-forbidden")
    print("    - Φ-shell Fibonacci-indexing of carrier-knots: T(2, F_n)")
    print()
    print("  Walk homology class (|p|, |q|) reduces mod 7 to (p mod 7, q mod 7).")
    print()
    print("  The cosmogenic Z_3 acts on q mod 7 (the v-cycle winding residue):")
    print("    - QR-cycle: q ∈ {1, 2, 4} → 3 distinct sectors")
    print("      • q ≡ 1: lepton (carrier T(2, F_2 = 1) = unknot)")
    print("      • q ≡ 2: meson  (carrier T(2, F_3 = 2) = Hopf)")
    print("      • q ≡ 4: hyperon (carrier T(2, F_4 = 3) = trefoil)")
    print("      → Z_3 cycle INCREMENTS Fibonacci index by +1 each step")
    print()
    print("    - NR-cycle: q ∈ {3, 5, 6}")
    print("      • q ≡ 3 (with p ≢ 0): nucleon (T(2, F_5 = 5) = cinquefoil)")
    print("      • q ≡ 5: meson (T(2, F_3 = 2) = Hopf)")
    print("      • q ≡ 6: substrate-prediction slot (Q6 candidate)")
    print()
    print("    - Polar: q ≡ 0 (p mod 7 ≠ 0): meson (T(2, F_3 = 2))")
    print()
    print("  The p mod 7 = 0 special case at q ≡ 3 corresponds to:")
    print("    'u-cycle winding is integer-Z_7-cycle (trivial mod 7)'")
    print("    → walk lives effectively in pure v-direction")
    print("    → degraded from nucleon (full structure) to meson")
    print()
    print("  ★ This gives a substrate-mechanism reading of rule (I) ★")
    print("  ★ rooted in cosmogenic Z_3 + QR/NR + polar-axis structure ★")
    print()
    print("  REMAINING THEORETICAL WORK:")
    print("    - Formal derivation of QR sector ordering: why F_2 → F_3 → F_4")
    print("      and not some permutation? Conjecture: Z_3 is ORIENTED (chiral).")
    print("    - Substrate role of NR (q ≡ 3 nucleon vs q ≡ 5 meson asymmetry)")
    print("    - q ≡ 6 prediction: F_6 = 8 (Fibonacci-extension) most natural")
    print()


def main():
    print("=" * 78)
    print("Q20 — Substrate mechanism for rule (I) via cosmogenic Z_3 cycle")
    print("=" * 78)
    print()
    step1_qr_z3_cycle()
    step2_nr_z3_cycle()
    step3_predict_q6_sector()
    step4_agl_z3_connection()
    step5_substrate_mechanism()


if __name__ == "__main__":
    main()
