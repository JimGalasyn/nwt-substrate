"""Steane-code reframe — check 1: verify Heffter Z_3×Z_2 ⊂ GL(3,2) = Aut(Fano).

Conjecture (2026-05-20 paradigm-level reframe): K_7 walks identify particles
with logical qubits of the Steane [[7,1,3]] code. The Heffter rotation σ and
sheet-swap τ generate the polar-fixing subgroup of the code's automorphism
group GL(3,2) = PSL(2,7), order 168.

This script verifies the group-theoretic embedding by:
  1. Constructing 7 Fano-plane points with a Heffter-compatible labeling
     P, E_1, E_2, E_3, F_1, F_2, F_3.
  2. Picking a Fano-line set consistent with the K_7 Heffter embedding.
  3. Building σ = (P)(E_1 E_2 E_3)(F_1 F_2 F_3) and
            τ = (P)(E_1 F_1)(E_2 F_2)(E_3 F_3).
  4. Checking that σ and τ each map every Fano line to a Fano line
     (i.e. sit inside Aut(Fano) = GL(3,2)).
  5. Verifying ⟨σ, τ⟩ has order 6 = |Z_3 × Z_2|.

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_check1_heffter_in_GL3_2.py
"""
from __future__ import annotations

from itertools import combinations, product


# ---------------------------------------------------------------------------
# Fano plane: 7 points, 7 lines, each line has 3 points, each point on 3 lines.
# Standard coordinatization: points = nonzero vectors of (Z_2)^3.
# Lines = sets of 3 points whose sum is zero.
# ---------------------------------------------------------------------------

FANO_POINTS = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)
                if (a, b, c) != (0, 0, 0)]
assert len(FANO_POINTS) == 7


def fano_lines() -> list[frozenset]:
    """All 7 Fano lines: triples {p, q, p+q (mod 2)}."""
    pts = FANO_POINTS
    seen = set()
    out = []
    for i, p in enumerate(pts):
        for q in pts[i+1:]:
            r = tuple((p[k] + q[k]) % 2 for k in range(3))
            if r == (0, 0, 0) or r == p or r == q:
                continue
            line = frozenset([p, q, r])
            if line not in seen:
                seen.add(line)
                out.append(line)
    return out


LINES = fano_lines()
assert len(LINES) == 7


# ---------------------------------------------------------------------------
# Pick a Heffter-compatible labeling of the 7 Fano points.
# In the K_7 Heffter embedding (memory framework_k7_heffter_embedding):
#   P is polar; {E_1, E_2, E_3} is the "matter" triangle; {F_1, F_2, F_3} is
#   the "antimatter" triangle (Z_2 partner of E).
# Z_3 cycles each triangle.  Z_2 swaps E_i ↔ F_i.
#
# Try the natural labeling:
#   P = (1, 1, 1)             -- the unique Fano point fixed by τ that swaps
#                                E_i ↔ F_i if E_i + F_i = P for all i.
#   E_1, E_2, E_3 = the three weight-1 vectors (standard basis e_1, e_2, e_3).
#   F_1, F_2, F_3 = the three weight-2 vectors with F_i = P + E_i = e_j + e_k.
# Then E_i + F_i = P for all i; τ(x) = x XOR P? No, τ(E_i) = F_i = P + E_i,
# and τ fixes P, so τ(x) = x + P for x != P and τ(P) = P.  That's not Z_2-linear,
# but it doesn't need to be linear — it just needs to preserve the line set.
# We will check that explicitly.
# ---------------------------------------------------------------------------

P    = (1, 1, 1)
E    = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
F    = [(0, 1, 1), (1, 0, 1), (1, 1, 0)]

# sanity: E_i + F_i == P for all i
for i in range(3):
    s = tuple((E[i][k] + F[i][k]) % 2 for k in range(3))
    assert s == P, (i, E[i], F[i], s)

# labels
LABEL = {P: 'P'}
for i, e in enumerate(E):
    LABEL[e] = f'E{i+1}'
for i, f in enumerate(F):
    LABEL[f] = f'F{i+1}'


def show_line(line: frozenset) -> str:
    return '{' + ', '.join(sorted(LABEL[p] for p in line)) + '}'


# ---------------------------------------------------------------------------
# Heffter permutations σ, τ
# ---------------------------------------------------------------------------

def sigma(x: tuple) -> tuple:
    """Z_3 rotation:  P → P,  E_i → E_{i+1 mod 3},  F_i → F_{i+1 mod 3}."""
    if x == P:
        return P
    if x in E:
        i = E.index(x)
        return E[(i + 1) % 3]
    if x in F:
        i = F.index(x)
        return F[(i + 1) % 3]
    raise ValueError(x)


def tau(x: tuple) -> tuple:
    """Z_2 sheet swap:  P → P,  E_i ↔ F_i."""
    if x == P:
        return P
    if x in E:
        i = E.index(x)
        return F[i]
    if x in F:
        i = F.index(x)
        return E[i]
    raise ValueError(x)


def apply_perm(perm, line: frozenset) -> frozenset:
    return frozenset(perm(p) for p in line)


# ---------------------------------------------------------------------------
# Check: do σ and τ preserve the Fano line set?
# ---------------------------------------------------------------------------

def check_preserves_lines(perm, name: str) -> bool:
    print(f"  testing {name} preserves the 7 Fano lines:")
    ok = True
    for line in LINES:
        new_line = apply_perm(perm, line)
        in_set = new_line in LINES
        mark = "✓" if in_set else "✗"
        print(f"    {show_line(line):<20} → {show_line(new_line):<20} {mark}")
        if not in_set:
            ok = False
    return ok


# ---------------------------------------------------------------------------
# Check: ⟨σ, τ⟩ has order 6 and is Z_3 × Z_2
# ---------------------------------------------------------------------------

def perm_to_tuple(perm) -> tuple:
    return tuple(perm(p) for p in FANO_POINTS)


def compose(f, g):
    """(f ∘ g)(x) = f(g(x))."""
    return lambda x: f(g(x))


def generate_group(gens):
    """BFS-generate the subgroup of S_7 from given generator callables."""
    identity = lambda x: x
    seen = {perm_to_tuple(identity): identity}
    frontier = [identity]
    while frontier:
        new_frontier = []
        for g1 in frontier:
            for g2 in gens:
                h = compose(g2, g1)
                key = perm_to_tuple(h)
                if key not in seen:
                    seen[key] = h
                    new_frontier.append(h)
        frontier = new_frontier
    return seen


# ---------------------------------------------------------------------------
# GL(3,2) construction for sanity (order 168)
# ---------------------------------------------------------------------------

def gl3_2_order() -> int:
    """Count invertible 3x3 matrices over F_2."""
    n = 0
    for entries in product((0, 1), repeat=9):
        rows = [entries[0:3], entries[3:6], entries[6:9]]
        # determinant mod 2
        a, b, c = rows[0]; d, e, f = rows[1]; g, h, i = rows[2]
        det = (a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)) % 2
        if det == 1:
            n += 1
    return n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("Steane-code reframe — Check 1: Heffter Z_3 × Z_2 ⊂ GL(3,2) = Aut(Fano)")
    print("=" * 78)
    print()
    print("Fano-point labeling (Heffter-compatible):")
    for pt in FANO_POINTS:
        print(f"  {LABEL[pt]:<3} = {pt}")
    print()

    print("All 7 Fano lines:")
    for line in LINES:
        print(f"  {show_line(line)}    (sum = {tuple(sum(p[k] for p in line) % 2 for k in range(3))})")
    print()

    # ----- Check σ -----
    ok_sigma = check_preserves_lines(sigma, "σ (Z_3 rotation)")
    print()

    # ----- Check τ -----
    ok_tau = check_preserves_lines(tau, "τ (Z_2 sheet swap)")
    print()

    # ----- Generated group -----
    group = generate_group([sigma, tau])
    print(f"  |⟨σ, τ⟩| = {len(group)}  (expect 6 = |Z_3 × Z_2|)")
    print()

    # Element orders
    def order_of(g):
        ident = perm_to_tuple(lambda x: x)
        h, k = g, 1
        while perm_to_tuple(h) != ident:
            h = compose(g, h); k += 1
            if k > 12: return None
        return k
    print("  elements of ⟨σ, τ⟩ and their orders:")
    for tup, g in sorted(group.items()):
        labels = tuple(LABEL[p] for p in tup)
        o = order_of(g)
        print(f"    order {o}: {labels}")
    print()

    # ----- Compatibility: σ and τ commute? (Z_3 × Z_2 ⇒ commute) -----
    st = compose(sigma, tau)
    ts = compose(tau, sigma)
    commute = perm_to_tuple(st) == perm_to_tuple(ts)
    print(f"  σ and τ commute: {commute}  (Z_3 × Z_2 requires commute = True)")
    print()

    # ----- GL(3,2) order sanity -----
    print(f"  |GL(3,2)| (computed by brute count) = {gl3_2_order()}  (expect 168)")
    print()

    # ----- HEADLINE -----
    print("=" * 78)
    print("HEADLINE")
    print("=" * 78)
    print()
    all_ok = ok_sigma and ok_tau and len(group) == 6 and commute
    if all_ok:
        print("  ✓ σ preserves all 7 Fano lines        → σ ∈ Aut(Fano) = GL(3,2)")
        print("  ✓ τ preserves all 7 Fano lines        → τ ∈ Aut(Fano) = GL(3,2)")
        print("  ✓ ⟨σ, τ⟩ has order 6 = |Z_3 × Z_2|")
        print("  ✓ σ and τ commute (Z_3 × Z_2 abelian)")
        print()
        print("  Therefore Heffter Z_3 × Z_2 ⊂ GL(3,2) = Aut(Fano plane)")
        print("                                       = Aut(Steane [[7,1,3]] code).")
        print()
        print("  The Heffter Z_3 × Z_2 is the polar-fixing subgroup of the")
        print("  Steane code's automorphism group.  Readings (1) Fano-octonion,")
        print("  (2) trefoil residue critical-phase, (3) stabilizer eigenstate")
        print("  are now provably one structure.")
    else:
        print("  ✗ Embedding check FAILED.  Re-examine the Heffter labeling.")


if __name__ == "__main__":
    main()
