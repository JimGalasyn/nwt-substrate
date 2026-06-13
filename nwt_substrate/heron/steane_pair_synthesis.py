"""Heron Exp 12 — Steane pair-synthesis circuits for substrate-monism's
syndrome-classification protocol.

Builds a Qiskit `QuantumCircuit` for each compendium (p, q) class that
implements:
  1. Steane [[7, 1, 3]] |0_L> preparation (3 H + 9 CX standard encoder)
  2. Walk's direction-sensitive Pauli word from the lookup table
     (X-on-destination for QR-steps; Z-on-destination for NR-steps)
  3. 6 stabilizer-syndrome measurements (3 X-type + 3 Z-type) via ancillas
  4. Logical-Z transversal readout (destructive by default; ancilla-coupled
     option via flag)

Lookup table: `analysis/steane_exp12_walk_pauli_table.json` (VV Check 3 + 5).

Circuit layout (default):
  - 7 system qubits (Steane block): indices 0..6 (Hamming positions 1..7)
  - 6 syndrome ancillas: indices 7..12 (one per stabilizer generator)
  - 1 optional logical-Z ancilla: index 13 (ancilla-coupled mode only)
  - Classical bits: 6 (syndromes) + 7 (logical Z destructive) or
                    6 + 1 (logical-Z ancilla-coupled)

Usage:
    from nwt_substrate.heron.steane_pair_synthesis import (
        build_exp12_circuit, build_control_circuit, COMPENDIUM_LOOKUP,
    )
    # Build a circuit for the proton's (1, 3) walk
    qc = build_exp12_circuit(p=1, q=3, logical_z_mode='destructive')
    # Build an identity control
    qc = build_control_circuit('identity')

Run as script for a sanity print:
    PYTHONPATH=/home/jim/repos/nwt-substrate python3 -m \\
        nwt_substrate.heron.steane_pair_synthesis
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

from nwt_substrate.isa.constants import N_VERTICES_K7

try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    HAS_QISKIT = True
except ImportError:
    HAS_QISKIT = False


# Path to the VV-generated lookup table
_LOOKUP_PATH = (Path(__file__).parents[2] / "analysis"
                / "steane_exp12_walk_pauli_table.json")


def _load_lookup() -> dict:
    """Load the Steane Exp 12 walk-to-Pauli lookup table from disk.

    Returns
    -------
    dict
        Lookup table mapping walk coordinates to Pauli-word bit arrays.
    """
    with open(_LOOKUP_PATH) as f:
        return json.load(f)


COMPENDIUM_LOOKUP = _load_lookup() if _LOOKUP_PATH.exists() else None


# ---------------------------------------------------------------------------
# Steane code structural constants
# ---------------------------------------------------------------------------
# Hamming (7, 4) parity-check matrix rows = Steane stabilizer supports
# Row k = positions q (0-indexed) where bit-k of (q+1) is set
STEANE_STAB_SUPPORTS = [
    [q for q in range(N_VERTICES_K7) if ((q + 1) >> k) & 1]
    for k in range(3)  # 3 = Hamming(7,4) parity bits (n − k); no substrate-constant alias
]
# Each row has weight 4 (anti-Fano-line). Stabilizers are X⊗⊗⊗ on
# supports[0..2] (X-type) and Z⊗⊗⊗ on supports[0..2] (Z-type).

# Logical X representative: weight-3 X on a Fano line
# Standard choice: qubits {0, 1, 2} (positions 1, 2, 3) — Fano line summing to 0
STEANE_LOGICAL_X_SUPPORT = [0, 1, 2]
STEANE_LOGICAL_Z_SUPPORT = [0, 1, 2]


# ---------------------------------------------------------------------------
# Steane |0_L> preparation
# ---------------------------------------------------------------------------

def _add_steane_zero_prep(qc: "QuantumCircuit", qubits: list[int]) -> None:
    """Add Steane [[7, 1, 3]] |0_L> preparation to circuit.

    Convention: qubits 0..6 are the Steane block (Hamming positions 1..7).
    The encoder uses 3 H gates on the "indicator" qubits (0, 1, 3 =
    positions 1, 2, 4) and 9 CNOTs to propagate.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to add the |0_L> preparation to (mutated in place).
    qubits : list[int]
        The 7 qubit indices of the Steane block.
    """
    # Step 1: H on the 3 indicator qubits (positions 1, 2, 4)
    qc.h(qubits[0])  # pos 1
    qc.h(qubits[1])  # pos 2
    qc.h(qubits[3])  # pos 4
    # Step 2: propagate via CNOTs. The pattern is determined by Hamming
    # parity-check positions: each indicator qubit drives the positions
    # that have its bit set.
    # bit-0 indicator (qubit 0) drives positions with bit-0 set in (pos - 1):
    #   targets are q2 (pos 3), q4 (pos 5), q6 (pos 7)
    qc.cx(qubits[0], qubits[2])
    qc.cx(qubits[0], qubits[4])
    qc.cx(qubits[0], qubits[6])
    # bit-1 indicator (qubit 1) drives positions with bit-1 set in (pos - 1):
    #   targets are q2 (pos 3), q5 (pos 6), q6 (pos 7)
    qc.cx(qubits[1], qubits[2])
    qc.cx(qubits[1], qubits[5])
    qc.cx(qubits[1], qubits[6])
    # bit-2 indicator (qubit 3) drives positions with bit-2 set in (pos - 1):
    #   targets are q4 (pos 5), q5 (pos 6), q6 (pos 7)
    qc.cx(qubits[3], qubits[4])
    qc.cx(qubits[3], qubits[5])
    qc.cx(qubits[3], qubits[6])


# ---------------------------------------------------------------------------
# Pauli-word application from lookup table
# ---------------------------------------------------------------------------

def _apply_pauli_word(qc: "QuantumCircuit", qubits: list[int],
                      x_bits: list[int], z_bits: list[int]) -> None:
    """Apply the walk's Pauli word: X_q if x_bits[q]; Z_q if z_bits[q].

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to apply the Pauli word to (mutated in place).
    qubits : list[int]
        Qubit indices of the Steane block.
    x_bits, z_bits : list[int]
        Per-qubit indicator bits for X and Z gates respectively.
    """
    for q in range(N_VERTICES_K7):
        if x_bits[q] and z_bits[q]:
            qc.y(qubits[q])  # Y up to phase (sufficient for stabilizer purposes)
        elif x_bits[q]:
            qc.x(qubits[q])
        elif z_bits[q]:
            qc.z(qubits[q])


# ---------------------------------------------------------------------------
# Stabilizer-syndrome measurement
# ---------------------------------------------------------------------------

def _add_syndrome_measurement(qc: "QuantumCircuit", sys_qubits: list[int],
                              anc_qubits: list[int], cbits: "ClassicalRegister") -> None:
    """Measure the 6 Steane stabilizer generators using 6 ancillas.

    For each X-type stabilizer X_a X_b X_c X_d on support S:
        H ancilla
        CNOT ancilla -> each q in S
        H ancilla
        measure ancilla

    For each Z-type stabilizer Z_a Z_b Z_c Z_d on support S:
        CNOT each q in S -> ancilla
        measure ancilla

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to add the syndrome measurement to (mutated in place).
    sys_qubits : list[int]
        The 7 system qubit indices.
    anc_qubits : list[int]
        The 6 ancilla qubit indices (X-type first, then Z-type).
    cbits : ClassicalRegister
        Classical register receiving the 6 syndrome bits.
    """
    n_stab = len(STEANE_STAB_SUPPORTS)   # stabilizer generators per type (X, then Z)
    # X-type stabilizers (ancillas 0 .. n_stab-1)
    for k in range(n_stab):
        anc = anc_qubits[k]
        cb = cbits[k]
        support = STEANE_STAB_SUPPORTS[k]
        qc.h(anc)
        for q_idx in support:
            qc.cx(anc, sys_qubits[q_idx])
        qc.h(anc)
        qc.measure(anc, cb)
    # Z-type stabilizers (ancillas n_stab .. 2*n_stab-1)
    for k in range(n_stab):
        anc = anc_qubits[n_stab + k]
        cb = cbits[n_stab + k]
        support = STEANE_STAB_SUPPORTS[k]
        for q_idx in support:
            qc.cx(sys_qubits[q_idx], anc)
        qc.measure(anc, cb)


# ---------------------------------------------------------------------------
# Logical-Z readout
# ---------------------------------------------------------------------------

def _add_logical_z_destructive(qc: "QuantumCircuit", sys_qubits: list[int],
                               cbits: "ClassicalRegister") -> None:
    """Destructive logical-Z readout: measure all 7 system qubits in Z basis,
    then post-process (compute parity on a Z_L representative support) on
    classical side.

    Writes all 7 single-qubit Z outcomes to cbits[0..6].

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to add the measurement to (mutated in place).
    sys_qubits : list[int]
        The 7 system qubit indices.
    cbits : ClassicalRegister
        Classical register receiving the 7 readout bits.
    """
    for q in range(N_VERTICES_K7):
        qc.measure(sys_qubits[q], cbits[q])


def _add_logical_z_ancilla(qc: "QuantumCircuit", sys_qubits: list[int],
                           anc_qubit: "Qubit", cbit: "Clbit") -> None:
    """Ancilla-coupled logical-Z: CNOT from each Z_L-support qubit to ancilla,
    measure ancilla. Non-destructive for the system qubits.

    Z_L = Z on qubits {0, 1, 2} (Fano line through positions 1, 2, 3).

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to add the measurement to (mutated in place).
    sys_qubits : list[int]
        The 7 system qubit indices.
    anc_qubit : Qubit
        Single ancilla qubit for the non-destructive Z_L measurement.
    cbit : Clbit
        Classical bit receiving the measurement outcome.
    """
    for q_idx in STEANE_LOGICAL_Z_SUPPORT:
        qc.cx(sys_qubits[q_idx], anc_qubit)
    qc.measure(anc_qubit, cbit)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_exp12_circuit(p: int, q: int,
                         logical_z_mode: Literal['destructive', 'ancilla'] = 'destructive',
                         ) -> "QuantumCircuit":
    """Build the Exp 12 syndrome + logical-Z circuit for the (|p|, |q|) walk.

    Args:
        p, q: walk winding (absolute values; lookup table is keyed on |p|, |q|)
        logical_z_mode: 'destructive' (default; measure all 7 Z) or
                        'ancilla' (single ancilla-coupled measurement)

    Returns:
        QuantumCircuit ready to be transpiled and submitted to Heron.
    """
    if not HAS_QISKIT:
        raise ImportError("qiskit is required for build_exp12_circuit")

    # Find the entry in the lookup table
    entry = next(
        (e for e in COMPENDIUM_LOOKUP["particles"]
          if e["p"] == abs(p) and e["q"] == abs(q)),
        None,
    )
    if entry is None:
        raise ValueError(f"No compendium entry for (|p|, |q|) = ({p}, {q})")

    return _build_circuit_from_pauli(
        entry["forward_pauli"]["x_part_bits"],
        entry["forward_pauli"]["z_part_bits"],
        label=f"({p},{q}):{','.join(entry['particles'])}",
        logical_z_mode=logical_z_mode,
    )


def build_control_circuit(kind: Literal['identity', 'x7', 'z7', 'h7'],
                           logical_z_mode: Literal['destructive', 'ancilla'] = 'destructive',
                           ) -> "QuantumCircuit":
    """Build a control circuit for pre-reg §5.5 validation runs.

    Args:
        kind:
            'identity' = no Pauli applied (expect (I, I), logical Z = +1)
            'x7'       = X⊗7 applied (logical X; expect (I, I), logical Z = -1)
            'z7'       = Z⊗7 applied (logical Z; expect (I, I), logical Z = +1)
            'h7'       = H⊗7 applied (logical Hadamard; expect (I, I), logical Z ≈ 0)
    """
    if not HAS_QISKIT:
        raise ImportError("qiskit is required for build_control_circuit")
    return _build_circuit_from_control(kind, logical_z_mode)


def _build_circuit_from_pauli(x_bits, z_bits, label: str,
                                logical_z_mode: str):
    sys = QuantumRegister(7, "sys")
    anc_synd = QuantumRegister(6, "anc_syn")
    cbits_syn = ClassicalRegister(6, "c_syn")
    if logical_z_mode == 'destructive':
        cbits_lz = ClassicalRegister(7, "c_lz")
        qc = QuantumCircuit(sys, anc_synd, cbits_syn, cbits_lz, name=label)
    else:
        anc_lz = QuantumRegister(1, "anc_lz")
        cbits_lz = ClassicalRegister(1, "c_lz")
        qc = QuantumCircuit(sys, anc_synd, anc_lz, cbits_syn, cbits_lz,
                             name=label)

    _add_steane_zero_prep(qc, sys)
    qc.barrier(label="pauli")
    _apply_pauli_word(qc, sys, x_bits, z_bits)
    qc.barrier(label="syn")
    _add_syndrome_measurement(qc, sys, anc_synd, cbits_syn)
    qc.barrier(label="lz")
    if logical_z_mode == 'destructive':
        _add_logical_z_destructive(qc, sys, cbits_lz)
    else:
        _add_logical_z_ancilla(qc, sys, anc_lz[0], cbits_lz[0])
    return qc


def _build_circuit_from_control(kind, logical_z_mode):
    sys = QuantumRegister(7, "sys")
    anc_synd = QuantumRegister(6, "anc_syn")
    cbits_syn = ClassicalRegister(6, "c_syn")
    if logical_z_mode == 'destructive':
        cbits_lz = ClassicalRegister(7, "c_lz")
        qc = QuantumCircuit(sys, anc_synd, cbits_syn, cbits_lz,
                             name=f"control:{kind}")
    else:
        anc_lz = QuantumRegister(1, "anc_lz")
        cbits_lz = ClassicalRegister(1, "c_lz")
        qc = QuantumCircuit(sys, anc_synd, anc_lz, cbits_syn, cbits_lz,
                             name=f"control:{kind}")

    _add_steane_zero_prep(qc, sys)
    qc.barrier(label="pauli")
    if kind == 'identity':
        pass  # no Pauli
    elif kind == 'x7':
        for q in range(N_VERTICES_K7): qc.x(sys[q])
    elif kind == 'z7':
        for q in range(N_VERTICES_K7): qc.z(sys[q])
    elif kind == 'h7':
        for q in range(N_VERTICES_K7): qc.h(sys[q])
    else:
        raise ValueError(f"unknown control: {kind}")
    qc.barrier(label="syn")
    _add_syndrome_measurement(qc, sys, anc_synd, cbits_syn)
    qc.barrier(label="lz")
    if logical_z_mode == 'destructive':
        _add_logical_z_destructive(qc, sys, cbits_lz)
    else:
        _add_logical_z_ancilla(qc, sys, anc_lz[0], cbits_lz[0])
    return qc


# ---------------------------------------------------------------------------
# Post-processing helpers
# ---------------------------------------------------------------------------

def decode_syndrome(syndrome_bits: list[int]) -> tuple[int | None, int | None]:
    """Decode 6-bit syndrome (s_X^1, s_X^2, s_X^3, s_Z^1, s_Z^2, s_Z^3) →
    (X-error qubit, Z-error qubit). Returns None for any zero 3-bit
    sub-syndrome (= no error on that rail).

    Caller convention: cbits[0..2] = X-syndrome bits, cbits[3..5] = Z-syndrome.
    Hamming decoding: 3-bit syndrome s ↦ qubit index s_0 + 2 s_1 + 4 s_2 - 1
    (None if all zero).
    """
    sx = syndrome_bits[0] + 2 * syndrome_bits[1] + 4 * syndrome_bits[2]
    sz = syndrome_bits[3] + 2 * syndrome_bits[4] + 4 * syndrome_bits[5]
    return (sx - 1 if sx > 0 else None, sz - 1 if sz > 0 else None)


def decode_logical_z_destructive(z_bits: list[int]) -> int:
    """Compute logical-Z eigenvalue from 7 destructive single-qubit Z
    measurements. Logical-Z eigenvalue ∈ {+1, -1}.

    Use Z_L = Z on qubits {0, 1, 2} (Fano line through positions 1, 2, 3):
    eigenvalue = (-1)^(z_bits[0] + z_bits[1] + z_bits[2]).
    """
    parity = sum(z_bits[i] for i in STEANE_LOGICAL_Z_SUPPORT) % 2
    return -1 if parity else +1


def decode_logical_z_ancilla(anc_bit: int) -> int:
    return -1 if anc_bit else +1


# ---------------------------------------------------------------------------
# Script-mode sanity print
# ---------------------------------------------------------------------------

def main():
    if not HAS_QISKIT:
        print("qiskit not installed; cannot build circuits.")
        return
    if COMPENDIUM_LOOKUP is None:
        print(f"Lookup table not found at {_LOOKUP_PATH}.")
        return

    print(f"Loaded {len(COMPENDIUM_LOOKUP['particles'])} (p, q) classes "
           f"from lookup table.")
    print()

    # Sanity: build proton, electron, mu- circuits + identity control
    samples = [
        ("proton", 1, 3),
        ("electron", 2, 1),
        ("muon", 1, 8),
        ("pi+", 3, 5),
        ("D0", 3, 7),
    ]
    for name, p, q in samples:
        qc = build_exp12_circuit(p, q, logical_z_mode='destructive')
        depth = qc.depth()
        cx_count = qc.count_ops().get("cx", 0)
        h_count = qc.count_ops().get("h", 0)
        sg_count = sum(qc.count_ops().get(g, 0) for g in ("x", "y", "z"))
        meas_count = qc.count_ops().get("measure", 0)
        print(f"  {name:<10} (p,q)=({p},{q}): "
               f"depth={depth}, CX={cx_count}, H={h_count}, "
               f"single-Pauli={sg_count}, measure={meas_count}")

    print()
    for ctrl in ('identity', 'x7', 'z7', 'h7'):
        qc = build_control_circuit(ctrl, logical_z_mode='destructive')
        print(f"  control {ctrl:<8}: depth={qc.depth()}, "
               f"CX={qc.count_ops().get('cx', 0)}, "
               f"H={qc.count_ops().get('h', 0)}, "
               f"single-Pauli={sum(qc.count_ops().get(g, 0) for g in ('x','y','z'))}")

    # Spot-check: build proton with ancilla-coupled logical-Z
    qc = build_exp12_circuit(1, 3, logical_z_mode='ancilla')
    print()
    print(f"  proton (ancilla LZ): depth={qc.depth()}, "
           f"CX={qc.count_ops().get('cx', 0)}, "
           f"total qubits={qc.num_qubits}")


if __name__ == "__main__":
    main()
