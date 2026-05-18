"""
B_9 K_7-edge investigation — follow-up to Tier B.5 bonus finding.

The Wade's-rules audit (FT B.5) found that the tricapped trigonal prism
(TTP) deltahedron — the structure of the closo borane B_9 H_9^{2-} —
has E = 21 edges, exactly matching N_EDGES_K7 = dim(so(7)) = 21.

Question: is this a STRUCTURAL realization (TTP edge graph encodes K_7 /
so(7) algebra), or just a numerical coincidence at the edge count?

Audit: compare TTP and K_7 on (a) vertex count, (b) edge count, (c)
subgraph containment K_7 ⊂ TTP, (d) chromatic number, (e) automorphism
group / edge orbit decomposition, (f) adjacency-matrix spectrum.

If any structural property matches → suggestive substrate signal.
If only edge counts match → numerical coincidence only.

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_b9_k7_edge_investigation.py
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Build adjacency matrices
# ---------------------------------------------------------------------------

def build_k7_adjacency():
    """K_7: complete graph on 7 vertices, all 21 edges present."""
    n = 7
    A = np.ones((n, n), dtype=int) - np.eye(n, dtype=int)
    return A


def build_ttp_adjacency():
    """Tricapped trigonal prism (B_9 closo deltahedron).

    Vertex labels:
        0, 1, 2 — top triangle
        3, 4, 5 — bottom triangle
        6, 7, 8 — three caps (each capping a quadrilateral face of the prism)

    Edges:
        Top triangle: 0-1, 1-2, 0-2  (3 edges)
        Bottom triangle: 3-4, 4-5, 3-5  (3 edges)
        Vertical prism: 0-3, 1-4, 2-5  (3 edges)
        Cap-to-top + cap-to-bottom: each cap caps one of the three
          quadrilateral faces (e.g., 0-1-4-3, 1-2-5-4, 0-2-5-3) and
          connects to all four corners. So:
            Cap 6 caps face 0-1-4-3: edges 6-0, 6-1, 6-3, 6-4
            Cap 7 caps face 1-2-5-4: edges 7-1, 7-2, 7-4, 7-5
            Cap 8 caps face 0-2-5-3: edges 8-0, 8-2, 8-3, 8-5
          12 edges total

    Total: 3 + 3 + 3 + 12 = 21 edges. χ = 9 − 21 + 14 = 2 (genus 0). ✓
    """
    edges = []
    # Top triangle
    edges += [(0, 1), (1, 2), (0, 2)]
    # Bottom triangle
    edges += [(3, 4), (4, 5), (3, 5)]
    # Vertical prism
    edges += [(0, 3), (1, 4), (2, 5)]
    # Cap 6 caps face (0,1,4,3)
    edges += [(6, 0), (6, 1), (6, 3), (6, 4)]
    # Cap 7 caps face (1,2,5,4)
    edges += [(7, 1), (7, 2), (7, 4), (7, 5)]
    # Cap 8 caps face (0,2,5,3)
    edges += [(8, 0), (8, 2), (8, 3), (8, 5)]

    n = 9
    A = np.zeros((n, n), dtype=int)
    for u, v in edges:
        A[u, v] = 1
        A[v, u] = 1
    return A, edges


# ---------------------------------------------------------------------------
# Graph invariants
# ---------------------------------------------------------------------------

def edge_count(A):
    return int(A.sum()) // 2


def euler_characteristic(V, E, F):
    return V - E + F


def chromatic_number_brute(A):
    """Find chromatic number by brute-force coloring."""
    n = A.shape[0]
    for k in range(1, n + 1):
        if _is_k_colorable(A, k):
            return k
    return n


def _is_k_colorable(A, k):
    """Greedy + DFS test for k-colorability."""
    n = A.shape[0]
    coloring = [-1] * n

    def is_safe(v, color):
        for u in range(n):
            if A[v, u] == 1 and coloring[u] == color:
                return False
        return True

    def graph_color(v):
        if v == n:
            return True
        for color in range(k):
            if is_safe(v, color):
                coloring[v] = color
                if graph_color(v + 1):
                    return True
                coloring[v] = -1
        return False

    return graph_color(0)


def has_k7_subgraph(A):
    """Does A contain K_7 as a subgraph? (All 21 edges of some 7-vertex set
    are present in A.)"""
    n = A.shape[0]
    if n < 7:
        return False
    from itertools import combinations
    for subset in combinations(range(n), 7):
        # Check all C(7,2) = 21 pairs are edges
        all_edges = True
        for i, j in combinations(subset, 2):
            if A[i, j] != 1:
                all_edges = False
                break
        if all_edges:
            return True
    return False


def adjacency_spectrum(A):
    """Sorted eigenvalues of adjacency matrix (decreasing)."""
    eigs = np.linalg.eigvalsh(A.astype(float))
    return np.sort(eigs)[::-1]


def degree_sequence(A):
    """Sorted vertex degrees (decreasing)."""
    return sorted(A.sum(axis=1).tolist(), reverse=True)


# ---------------------------------------------------------------------------
# Edge-orbit decomposition under TTP D_3h symmetry
# ---------------------------------------------------------------------------

def ttp_edge_orbits():
    """TTP has D_3h symmetry. Edges decompose into orbits under the
    3-fold rotation + horizontal mirror.

    By vertex labeling above:
      - Top triangle: {(0,1), (1,2), (0,2)}
      - Bottom triangle: {(3,4), (4,5), (3,5)}
      - Vertical prism: {(0,3), (1,4), (2,5)}
      - Cap-to-top: 6 edges (each cap connects to 2 top vertices)
        {(6,0),(6,1)}, {(7,1),(7,2)}, {(8,0),(8,2)}
      - Cap-to-bottom: 6 edges
        {(6,3),(6,4)}, {(7,4),(7,5)}, {(8,3),(8,5)}

    Under horizontal mirror, top↔bottom triangles swap, and cap-to-top
    ↔ cap-to-bottom swap. So:
      - {top triangle, bottom triangle}: 2-orbit of 3 each, but merge
        under horizontal mirror into single orbit of 6 edges
      - vertical prism: 3-orbit (fixed under mirror)
      - {cap-to-top, cap-to-bottom}: merge under mirror → 12-edge orbit

    Total orbits under D_3h: 3 orbits, sizes (6, 3, 12) = 21
    Total orbits under C_3v only (no horizontal mirror): 5 orbits (3,3,3,6,6) = 21
    """
    return {
        "D_3h (full TTP symmetry)": [
            ("triangle (top+bottom)", 6),
            ("vertical prism", 3),
            ("cap (top+bottom)", 12),
        ],
        "C_3v (no horizontal mirror)": [
            ("top triangle", 3),
            ("bottom triangle", 3),
            ("vertical prism", 3),
            ("cap-to-top", 6),
            ("cap-to-bottom", 6),
        ],
    }


# ---------------------------------------------------------------------------
# Main audit
# ---------------------------------------------------------------------------

def run_audit():
    K7 = build_k7_adjacency()
    TTP, ttp_edges = build_ttp_adjacency()

    return {
        "K7": {
            "V": K7.shape[0],
            "E": edge_count(K7),
            "degree_sequence": degree_sequence(K7),
            "chromatic_number": chromatic_number_brute(K7),
            "spectrum": adjacency_spectrum(K7).tolist(),
            "has_K7_subgraph_in_self": has_k7_subgraph(K7),
        },
        "TTP": {
            "V": TTP.shape[0],
            "E": edge_count(TTP),
            "degree_sequence": degree_sequence(TTP),
            "chromatic_number": chromatic_number_brute(TTP),
            "spectrum": adjacency_spectrum(TTP).tolist(),
            "has_K7_subgraph": has_k7_subgraph(TTP),
            "edge_orbits": ttp_edge_orbits(),
        },
    }


def render_report(audit):
    out = []
    out.append("=" * 78)
    out.append("B_9 K_7-edge investigation (bonus finding from Tier B.5)")
    out.append("=" * 78)
    out.append("")
    out.append("Question: does B_9 tricapped trigonal prism (TTP) deltahedron realize")
    out.append("the K_7 edge graph, or is the 21-edge match a numerical coincidence?")
    out.append("")

    k7 = audit["K7"]
    ttp = audit["TTP"]

    # Comparison table
    out.append("-" * 78)
    out.append("GRAPH-INVARIANT COMPARISON")
    out.append("-" * 78)
    out.append("")
    out.append(f"  {'invariant':<28} {'K_7':>15} {'TTP (B_9)':>15} {'match?':>10}")
    out.append("  " + "-" * 70)

    out.append(f"  {'vertices V':<28} {k7['V']:>15} {ttp['V']:>15} {'NO' if k7['V'] != ttp['V'] else 'YES':>10}")
    out.append(f"  {'edges E':<28} {k7['E']:>15} {ttp['E']:>15} {'YES' if k7['E'] == ttp['E'] else 'NO':>10}")
    out.append(f"  {'chromatic number χ(G)':<28} {k7['chromatic_number']:>15} {ttp['chromatic_number']:>15} {'NO' if k7['chromatic_number'] != ttp['chromatic_number'] else 'YES':>10}")
    out.append(f"  {'max degree':<28} {max(k7['degree_sequence']):>15} {max(ttp['degree_sequence']):>15} {'NO' if max(k7['degree_sequence']) != max(ttp['degree_sequence']) else 'YES':>10}")
    out.append(f"  {'min degree':<28} {min(k7['degree_sequence']):>15} {min(ttp['degree_sequence']):>15} {'NO' if min(k7['degree_sequence']) != min(ttp['degree_sequence']) else 'YES':>10}")
    out.append(f"  {'K_7 subgraph?':<28} {'(self)':>15} {str(ttp['has_K7_subgraph']):>15}")
    out.append("")

    out.append(f"  K_7 degree sequence:  {k7['degree_sequence']}")
    out.append(f"  TTP degree sequence:  {ttp['degree_sequence']}")
    out.append("")
    out.append(f"  K_7 spectrum:  {[round(x, 4) for x in k7['spectrum']]}")
    out.append(f"  TTP spectrum:  {[round(x, 4) for x in ttp['spectrum']]}")
    out.append("")

    # TTP edge orbits
    out.append("-" * 78)
    out.append("TTP EDGE ORBIT DECOMPOSITION")
    out.append("-" * 78)
    for sym_group, orbits in ttp["edge_orbits"].items():
        out.append(f"  Under {sym_group}:")
        for name, size in orbits:
            out.append(f"    {name:<30}: {size} edges")
        total = sum(s for _, s in orbits)
        out.append(f"    {'TOTAL':<30}: {total} edges")
        out.append("")

    out.append("  K_7 edge orbits under S_7 (full automorphism group):")
    out.append("    transitive: 1 orbit of 21 edges (all edges equivalent)")
    out.append("")

    # Verdict
    out.append("=" * 78)
    out.append("VERDICT")
    out.append("=" * 78)
    out.append("")
    out.append("  Vertex count mismatch:   K_7 has 7, TTP has 9. Cannot be graph-isomorphic.")
    out.append(f"  K_7 subgraph in TTP?     {ttp['has_K7_subgraph']}.")
    out.append(f"  Chromatic number:        K_7 = {k7['chromatic_number']} vs TTP = {ttp['chromatic_number']}. NOT MATCHING.")
    out.append(f"  Edge orbit structure:    K_7 = 1 orbit (S_7-transitive)")
    out.append(f"                            TTP = 3 orbits (D_3h) or 5 orbits (C_3v). NOT MATCHING.")
    out.append(f"  Spectrum:                K_7 = {{6 (×1), -1 (×6)}}, TTP has different spectrum.")
    out.append("")
    out.append("  CONCLUSION:")
    out.append("    The 21-edge match between TTP and K_7 is NOT a structural realization.")
    out.append("    Every graph invariant beyond edge count DIFFERS:")
    out.append("      - K_7 is 7-vertex complete (regular degree 6, 1 transitive edge orbit)")
    out.append("      - TTP is 9-vertex sparse deltahedron (degree 4 or 5, 3-5 edge orbits)")
    out.append("")
    out.append("    The 21 = 21 match is at the INTEGER level only — both happen to")
    out.append("    equal N_EDGES_K7 = dim(so(7)) = 21, but the graphs encode different")
    out.append("    algebraic structure. The substrate identifies B_9 via PRIMITIVE")
    out.append("    integers (n=9=N_POS_ROOTS_SO7 and E=21=N_EDGES_K7) but NOT via")
    out.append("    graph realization.")
    out.append("")
    out.append("  IMPLICATION:")
    out.append("    Bonus B.5 finding 'B_9 deltahedron has 21 edges' is HONESTLY a cross-arc")
    out.append("    INTEGER REUSE (substrate-canonical integer 21 appears in both K_7 and")
    out.append("    TTP edge counts), not a graph-theoretic substrate realization. Same")
    out.append("    pattern as A.3 shell-32 = K_7_TRIANGLES − RANK_SO7: integer coincidence")
    out.append("    at the substrate-canonical level, not structural identity.")
    out.append("")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
