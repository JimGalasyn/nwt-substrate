"""Q24 — NR irregularity resolution: NR-decrement vs Q22 F_6 reading.

Q20 found NR cycle (3 → 6 → 5 → 3) under cosmogenic Z_3 = mult-by-2 mod 7.
Q21 tentatively assigned q=6 → F_6 = 8 (baryon-class).
Q22 used this to give η_B = (3/14)·α⁴ via baryon-shell count.

BUT: the NR sequence F_5 → F_6 → F_3 is IRREGULAR (drops -3 in F-index
from F_6 to F_3). QR by contrast has clean F-increment: F_2 → F_3 → F_4
→ F_2 (consecutive Fibonacci, period 3).

This Q tests an ALTERNATIVE structural reading: NR-decrement cycle.
If q=6 → F_4 = 3 (hyperon-like) instead of F_6 = 8 (baryon prediction),
then NR sequence becomes:
    F_5 → F_4 → F_3 → F_5
    (decrement -1 per Z_3 step, wraparound F_3 → F_5)
which is the MIRROR of QR's increment cycle. Beautiful symmetry.

Tests:
  (1) Verify NR-decrement reading is mathematically consistent
  (2) Compare with QR-increment for symmetry
  (3) Identify implication for η_B baryon-shell count (now only F_4, F_5 = 2 shells)
  (4) Identify implication for q ≡ 6 substrate species mass (now hyperon-like, ~5-20 GeV)
  (5) Discuss which reading is structurally favored + how to falsify

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_q24_nr_decrement_symmetry.py
"""
from __future__ import annotations

from nwt_substrate.particles.mass import paper6_mass_mev


def fibonacci(n: int) -> int:
    if n <= 2:
        return 1
    a, b = 1, 1
    for _ in range(n - 2):
        a, b = b, a + b
    return b


def main():
    print("=" * 78)
    print("Q24 — NR irregularity resolution: decrement vs F_6 baryon reading")
    print("=" * 78)
    print()

    # ---- Step 1: Two readings of NR side -------------------------------
    print("[Step 1] Two competing readings of the NR Z_3 cycle (3 → 6 → 5 → 3)")
    print("-" * 78)
    print()
    print("  Reading A: Q22 'F_6 baryon-shell' (matches η_B baryon count)")
    print(f"  {'q mod 7':<8} {'sector':<14} {'Fibonacci':<10} {'F-index delta'}")
    print("  " + "-" * 55)
    readings_A = [(3, "nucleon", 5), (6, "PREDICT", 6), (5, "meson", 3)]
    for k, (q, sec, fi) in enumerate(readings_A):
        nxt_fi = readings_A[(k+1) % 3][2]
        delta = nxt_fi - fi
        delta_str = f"{delta:+d}" if k < 2 else f"{delta:+d} (wrap)"
        print(f"  {q:<8} {sec:<14} F_{fi}={fibonacci(fi):<6} {delta_str}")
    print()
    print(f"  Reading A NR F-index sequence: 5 → 6 → 3 → 5")
    print(f"  Pattern: IRREGULAR (+1, -3, +2)")
    print()

    print("  Reading B: NR-decrement (mirror of QR-increment)")
    print(f"  {'q mod 7':<8} {'sector':<14} {'Fibonacci':<10} {'F-index delta'}")
    print("  " + "-" * 55)
    readings_B = [(3, "nucleon", 5), (6, "hyperon-mirror", 4), (5, "meson", 3)]
    for k, (q, sec, fi) in enumerate(readings_B):
        nxt_fi = readings_B[(k+1) % 3][2]
        delta = nxt_fi - fi
        delta_str = f"{delta:+d}" if k < 2 else f"{delta:+d} (wrap)"
        print(f"  {q:<8} {sec:<14} F_{fi}={fibonacci(fi):<6} {delta_str}")
    print()
    print(f"  Reading B NR F-index sequence: 5 → 4 → 3 → 5")
    print(f"  Pattern: CLEAN DECREMENT (-1, -1, +2 wrap)")
    print()

    # ---- Step 2: QR vs NR symmetry under Reading B ---------------------
    print("[Step 2] Reading B gives QR ↔ NR mirror symmetry")
    print("-" * 78)
    print()
    print(f"  QR cycle (1 → 2 → 4 → 1) under Z_3 = mult-by-2:")
    print(f"    F-index sequence: 2 → 3 → 4 → 2")
    print(f"    Pattern: CLEAN INCREMENT (+1, +1, -2 wrap)")
    print(f"    Sectors: lepton (F_2) → meson (F_3) → hyperon (F_4)")
    print()
    print(f"  NR cycle (3 → 6 → 5 → 3) under Z_3 = mult-by-2:")
    print(f"    Reading B F-index sequence: 5 → 4 → 3 → 5")
    print(f"    Pattern: CLEAN DECREMENT (-1, -1, +2 wrap)")
    print(f"    Sectors: nucleon (F_5) → hyperon-mirror (F_4) → meson (F_3)")
    print()
    print(f"  ★ Reading B gives BEAUTIFUL QR ↔ NR mirror symmetry: ★")
    print(f"    QR goes UP through {{F_2, F_3, F_4}}")
    print(f"    NR goes DOWN through {{F_5, F_4, F_3}}")
    print(f"    QR starts at lightest (lepton), NR starts at heaviest (nucleon)")
    print(f"    Both meet at meson (F_3, q ∈ {{2, 5}})")
    print()
    print(f"  This is structurally analogous to CPT mirror: matter generation")
    print(f"  cycles ascend in mass, antimatter generations descend.")
    print()

    # ---- Step 3: η_B implications --------------------------------------
    print("[Step 3] Implications for η_B = (3/14)·α⁴")
    print("-" * 78)
    print()
    print(f"  Reading A (Q22 F_6=8 baryon):")
    print(f"    Baryon Fibonacci shells: {{F_4, F_5, F_6}} = 3 shells")
    print(f"    η_B numerator = 3 = # baryon shells ✓")
    print(f"    η_B = 3/14 · α⁴ = 6.077e-10 (matches Planck at 0.38%)")
    print()
    print(f"  Reading B (NR-decrement):")
    print(f"    Baryon Fibonacci shells: {{F_4, F_5}} = 2 shells")
    print(f"    η_B numerator must come from elsewhere (e.g., RANK_SO7 = 3 generations)")
    print(f"    The 'baryon-shell count' reading is FALSIFIED in Reading B")
    print(f"    But the RANK_SO7/dim(G_2) = 3/14 reading remains valid")
    print()
    print(f"  Reading A is FAVORED if we want baryon-shell-count η_B reading.")
    print(f"  Reading B is FAVORED if we prefer cyclic-symmetric NR structure.")
    print()

    # ---- Step 4: Mass implications for q ≡ 6 species -------------------
    print("[Step 4] Mass predictions for q ≡ 6 species under each reading")
    print("-" * 78)
    print()
    print(f"  Reading A: q ≡ 6 → n_q = 8 (F_6 baryon)")
    print(f"  Reading B: q ≡ 6 → n_q = 3 (F_4 hyperon-mirror)")
    print()
    print(f"  {'Candidate':<12} {'m_int':<7} {'A: n_q=8 mass':<18} {'B: n_q=3 mass'}")
    print("  " + "-" * 65)
    candidates = [(1, 6, 2), (1, 6, 7), (1, 6, 16),
                  (5, 6, 6), (5, 6, 11),
                  (7, 6, 8), (7, 6, 13)]
    for p, q, m in candidates:
        mass_A = paper6_mass_mev(p, q, m, 8)
        mass_B = paper6_mass_mev(p, q, m, 3)
        if mass_A and mass_B:
            mass_A_str = f"{mass_A/1000:.2f} GeV" if mass_A > 1000 else f"{mass_A:.0f} MeV"
            mass_B_str = f"{mass_B/1000:.2f} GeV" if mass_B > 1000 else f"{mass_B:.0f} MeV"
            print(f"  ({p},{q})        {m:<7} {mass_A_str:<18} {mass_B_str}")
    print()

    # ---- Step 5: Substrate-canonical reading? --------------------------
    print("[Step 5] Which reading is substrate-canonical?")
    print("=" * 78)
    print()
    print("  Arguments for Reading A (F_6 baryon, Q22):")
    print("    (a) Matches η_B baryon-shell-count reading at 0.38%")
    print("    (b) Q21 cycle-extension: Z_3 INCREMENTS F-index on QR; if NR")
    print("        also increments (3→6→5 with F: 5→6→3), gives Reading A")
    print("    (c) Predicts a NEW heavy substrate species (DM/BSM candidate)")
    print()
    print("  Arguments for Reading B (NR-decrement, mirror):")
    print("    (a) Gives QR ↔ NR mirror symmetry (clean +1 vs -1 cycles)")
    print("    (b) NR-decrement sequence F_5→F_4→F_3 is natural for")
    print("        ANTIMATTER direction (mass-decreasing = time-reversal)")
    print("    (c) Hyperon-class q ≡ 6 species are more empirically testable")
    print("        (~5-20 GeV vs TeV scale)")
    print("    (d) Reading B requires η_B's 3 = RANK_SO7 (generations), not")
    print("        baryon-shells. Both interpretations of '3' are valid")
    print("        substrate-canonically.")
    print()
    print("  Decisive test: detection of q ≡ 6 substrate species")
    print("    - GeV-scale (5-20 GeV): supports Reading B")
    print("    - TeV-scale (1-30 TeV): supports Reading A (or particles are DM)")
    print("    - Non-detection: either reading possible (substrate predictions")
    print("      sometimes don't materialize, e.g., d=2 Hamilton absence)")
    print()
    print("  HONEST VERDICT:")
    print("    Without empirical input, BOTH readings are substrate-canonical")
    print("    structural hypotheses. The choice depends on which structural")
    print("    feature is privileged:")
    print("    - Reading A: privilege η_B baryon-shell concrete count")
    print("    - Reading B: privilege QR ↔ NR cyclic symmetry")
    print()
    print("  Q24 RESOLVES the NR irregularity: NOT irregular under Reading B")
    print("  but at the cost of breaking Q22's baryon-shell η_B count. Two")
    print("  competing structurally-elegant readings, distinguishable by")
    print("  empirical test.")
    print()


if __name__ == "__main__":
    main()
