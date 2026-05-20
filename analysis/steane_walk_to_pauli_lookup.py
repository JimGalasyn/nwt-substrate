"""Generate the walk-to-Pauli lookup table for Exp 12 circuit construction.

For each compendium (|p|, |q|) class, output:
  - BFS-shortest K_7 walk
  - Direction-sensitive Pauli word (X-part, Z-part) per Check 3
  - Predicted (X-error, Z-error) Fano point pair per Check 5
  - Centralizer membership (True iff Pauli is logical, all-zeros syndrome)
  - K_7-vertex-to-Steane-qubit mapping for Heron implementation

Output: `steane_exp12_walk_pauli_table.json` (machine-readable for circuit
construction in `nwt_substrate.heron.steane_pair_synthesis`).

Run:
    cd /home/jim/repos/nwt-substrate
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 \\
        analysis/steane_walk_to_pauli_lookup.py
"""
from __future__ import annotations

import json
from pathlib import Path

from nwt_substrate.condensate import bfs_shortest_walks
from nwt_substrate.particles.compendium import COMPENDIUM


# K_7 vertex ↔ Fano point ↔ Steane qubit mapping (matches Check 3 / 5)
VERTEX_TO_FANO = {
    0: (1, 1, 1), 1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1),
    4: (0, 1, 1), 5: (1, 0, 1), 6: (1, 1, 0),
}
VERTEX_LABEL = {0: 'P', 1: 'E1', 2: 'E2', 3: 'E3', 4: 'F1', 5: 'F2', 6: 'F3'}


def std_position(point: tuple) -> int:
    """Map Fano point (b2, b1, b0) → integer Hamming position 1..7."""
    return 4 * point[0] + 2 * point[1] + 1 * point[2]


VERTEX_TO_QUBIT = {v: std_position(VERTEX_TO_FANO[v]) - 1 for v in range(7)}
QUBIT_TO_VERTEX = {q: v for v, q in VERTEX_TO_QUBIT.items()}

QR_DIRECTIONS = {1, 2, 4}
NR_DIRECTIONS = {3, 5, 6}


# Hamming (7,4) parity-check matrix rows = Steane X-stab supports
HAMMING_ROWS = [
    tuple((q + 1) >> k & 1 for q in range(7))
    for k in range(3)
]
X_STAB_GENS = HAMMING_ROWS
Z_STAB_GENS = HAMMING_ROWS


def walk_direction_pauli(walk):
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
    return tuple(x), tuple(z)


def syndrome(x_part, z_part):
    """Compute (X-syndrome, Z-syndrome) — 3 bits each."""
    s_X = tuple(sum(z_part[i]*g[i] for i in range(7)) % 2 for g in X_STAB_GENS)
    s_Z = tuple(sum(x_part[i]*g[i] for i in range(7)) % 2 for g in Z_STAB_GENS)
    return s_X, s_Z


def syndrome_to_qubit(syn):
    idx = syn[0] + 2*syn[1] + 4*syn[2]
    return idx - 1 if idx > 0 else None


def qubit_to_label(q):
    if q is None: return "I"
    return VERTEX_LABEL[QUBIT_TO_VERTEX[q]]


def pauli_word_string(x_part, z_part):
    """Render Pauli word as 'X_qZ_q' notation for the operator string."""
    parts = []
    for q in range(7):
        if x_part[q] and z_part[q]:
            parts.append(f"Y_{q}")
        elif x_part[q]:
            parts.append(f"X_{q}")
        elif z_part[q]:
            parts.append(f"Z_{q}")
    return " ".join(parts) if parts else "I"


# Logical-Z support choice (matches steane_pair_synthesis helper):
# weight-3 Z on Fano line {qubit 0, 1, 2} = positions {1, 2, 3}
LOGICAL_Z_SUPPORT = (0, 1, 2)


def predicted_logical_z(x_part, z_part) -> int:
    """Predicted logical-Z eigenvalue (+1 or -1) for the post-Pauli state.

    Z_L acting on the post-Pauli state |ψ⟩ = P_W|0_L⟩ has eigenvalue
    determined by whether P_W commutes with Z_L:
      - Commutes: eigenvalue = +1 (|0_L⟩-like)
      - Anti-commutes: eigenvalue = -1 (|1_L⟩-like)
      - Mixed (P_W contains Y components, or X-part overlap with Z_L_support
        is odd while Z-part on Z_L_support qubits matter): need to check both
        parts.

    For our purposes:
      Z_L = Z_0 Z_1 Z_2.
      [Z_L, X_q] = 0 iff q ∉ {0, 1, 2}; anti-commutes iff q ∈ {0, 1, 2}.
      [Z_L, Z_q] = 0 always.
      [Z_L, Y_q] inherits from [Z_L, X_q] (Y is Pauli X with sign).
    So Z_L (anti-)commutes with P_W iff
      (X-part ∩ Z_L_support) is even-weighted.
    """
    overlap = sum(x_part[i] for i in LOGICAL_Z_SUPPORT)
    return +1 if overlap % 2 == 0 else -1


def main():
    walks = bfs_shortest_walks(max_length=25)

    table = {
        "metadata": {
            "source": "VV Check 3 + Check 5 (steane_check3, steane_check5)",
            "date": "2026-05-20",
            "convention": {
                "qubit_indexing": "Steane qubits 0..6 correspond to Hamming positions 1..7",
                "vertex_to_qubit": {str(v): VERTEX_TO_QUBIT[v]
                                       for v in range(7)},
                "vertex_to_fano_label": VERTEX_LABEL,
                "qubit_to_fano_label": {
                    str(q): VERTEX_LABEL[QUBIT_TO_VERTEX[q]] for q in range(7)
                },
                "QR_directions": sorted(QR_DIRECTIONS),
                "NR_directions": sorted(NR_DIRECTIONS),
                "pauli_rule": "X on destination qubit for QR-direction step; Z on destination for NR-direction step",
                "x_stab_supports_qubits": [
                    [q for q in range(7) if g[q]] for g in X_STAB_GENS
                ],
                "z_stab_supports_qubits": [
                    [q for q in range(7) if g[q]] for g in Z_STAB_GENS
                ],
                "logical_z_support_qubits": list(LOGICAL_Z_SUPPORT),
                "logical_z_eigenvalue_rule": "+1 if X-part overlap with Z_L support {0,1,2} is even-weighted, else -1",
            }
        },
        "particles": [],
    }

    # Group by (|p|, |q|) and gather all particles in each class
    by_pq = {}
    for entry in COMPENDIUM:
        key = (abs(entry["p"]), abs(entry["q"]))
        if key not in walks:
            continue
        by_pq.setdefault(key, []).append(entry)

    for (p, q) in sorted(by_pq.keys()):
        walk = walks[(p, q)]
        x_part, z_part = walk_direction_pauli(walk)
        s_X, s_Z = syndrome(x_part, z_part)
        q_X = syndrome_to_qubit(s_X)
        q_Z = syndrome_to_qubit(s_Z)
        in_centralizer = (q_X is None and q_Z is None)

        # Reverse walk
        rwalk = list(reversed(walk))
        rx_part, rz_part = walk_direction_pauli(rwalk)
        rs_X, rs_Z = syndrome(rx_part, rz_part)

        names = [e["name"] for e in by_pq[(p, q)]]

        table["particles"].append({
            "p": p, "q": q,
            "particles": names,
            "n_q": by_pq[(p, q)][0]["n_q"],
            "walk_K7_vertices": walk,
            "walk_steane_qubits": [VERTEX_TO_QUBIT[v] for v in walk],
            "walk_length": len(walk) - 1,
            "forward_pauli": {
                "x_part_bits": list(x_part),
                "z_part_bits": list(z_part),
                "pauli_string": pauli_word_string(x_part, z_part),
                "in_centralizer": in_centralizer,
                "x_syndrome_bits": list(s_X),
                "z_syndrome_bits": list(s_Z),
                "x_err_qubit": q_X,
                "z_err_qubit": q_Z,
                "x_err_fano": qubit_to_label(q_X),
                "z_err_fano": qubit_to_label(q_Z),
                "predicted_syndrome_fano_pair": [qubit_to_label(q_X),
                                                   qubit_to_label(q_Z)],
                "predicted_logical_z": predicted_logical_z(x_part, z_part),
            },
            "reverse_pauli": {
                "x_part_bits": list(rx_part),
                "z_part_bits": list(rz_part),
                "pauli_string": pauli_word_string(rx_part, rz_part),
                "x_syndrome_bits": list(rs_X),
                "z_syndrome_bits": list(rs_Z),
                "x_err_fano": qubit_to_label(syndrome_to_qubit(rs_X)),
                "z_err_fano": qubit_to_label(syndrome_to_qubit(rs_Z)),
                "predicted_logical_z": predicted_logical_z(rx_part, rz_part),
            },
        })

    # Cross-sector equivalence classes (groups of particle-classes sharing syndrome pair)
    equiv = {}
    for p in table["particles"]:
        sf = tuple(p["forward_pauli"]["predicted_syndrome_fano_pair"])
        equiv.setdefault(sf, []).append({
            "pq": [p["p"], p["q"]],
            "particles": p["particles"],
        })

    table["cross_sector_equivalences"] = [
        {
            "syndrome_pair": list(sf),
            "classes": cl,
            "n_particles": sum(len(c["particles"]) for c in cl),
        }
        for sf, cl in sorted(equiv.items())
    ]

    out_path = Path(__file__).parent / "steane_exp12_walk_pauli_table.json"
    with open(out_path, "w") as f:
        json.dump(table, f, indent=2)
    print(f"Wrote lookup table → {out_path}")
    print(f"  {len(table['particles'])} (p, q) classes")
    print(f"  {len(table['cross_sector_equivalences'])} distinct syndrome pairs")
    print()

    # Brief summary
    in_cent = sum(1 for p in table["particles"] if p["forward_pauli"]["in_centralizer"])
    print(f"Centralizer particles (logical operators, syndrome=(I,I)): {in_cent}")
    print(f"  These will be indistinguishable from one another and from identity")
    print(f"  via stabilizer syndromes alone. Logical readout needed for separation.")
    print()
    print("Cross-sector equivalences (syndrome pairs shared across (p,q) classes):")
    for eq in table["cross_sector_equivalences"]:
        all_particles = [n for c in eq["classes"] for n in c["particles"]]
        if len(eq["classes"]) > 1:
            print(f"  {tuple(eq['syndrome_pair'])}: {all_particles}")


if __name__ == "__main__":
    main()
