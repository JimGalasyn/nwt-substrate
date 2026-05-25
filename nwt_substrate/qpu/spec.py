"""Thin circuit IR + Steane builders + measurement-scheme strategies.

Our quantum vocabulary is small (Clifford prep + Pauli words + basis-change +
measurement), so we own a ~tiny IR and emit Qiskit or Braket from it. This is
option (a) of the design doc: zero new dependencies and it hands us the
canonical measurement layout for both SDKs for free.

`to_qiskit`/`to_braket` import their SDK lazily, so importing this module never
requires qiskit or braket to be installed.
"""
from __future__ import annotations

from dataclasses import dataclass, field


# --------------------------------------------------------------------------- IR
@dataclass(frozen=True)
class G:
    """A gate. name in {h,s,sdg,x,y,z,cx,cz,barrier}; qubits is its operands."""
    name: str
    qubits: tuple[int, ...] = ()


@dataclass(frozen=True)
class M:
    """Measure `qubit` in `basis` (z|x|y) into register[`position`]."""
    qubit: int
    basis: str
    register: str
    position: int


@dataclass
class CircuitSpec:
    n_qubits: int
    ops: list                      # list[G | M]
    registers: dict[str, int]      # name -> width
    name: str = "circ"

    def measured_layout(self) -> dict[str, tuple[int, ...]]:
        """register name -> tuple of qubit indices, in position order."""
        layout = {r: [None] * w for r, w in self.registers.items()}
        for op in self.ops:
            if isinstance(op, M):
                layout[op.register][op.position] = op.qubit
        out = {}
        for r, qs in layout.items():
            if any(q is None for q in qs):
                raise ValueError(f"register {r!r} has unmeasured positions: {qs}")
            out[r] = tuple(qs)
        return out

    def count_2q(self) -> int:
        return sum(1 for op in self.ops if isinstance(op, G) and len(op.qubits) == 2)

    # -- emitters -----------------------------------------------------------
    def to_qiskit(self):
        from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
        qr = QuantumRegister(self.n_qubits, "q")
        cregs = {n: ClassicalRegister(w, n) for n, w in self.registers.items()}
        qc = QuantumCircuit(qr, *cregs.values(), name=self.name)
        for op in self.ops:
            if isinstance(op, G):
                if op.name == "barrier":
                    qc.barrier()
                elif op.name in ("h", "s", "sdg", "x", "y", "z"):
                    getattr(qc, op.name)(qr[op.qubits[0]])
                elif op.name in ("cx", "cz"):
                    getattr(qc, op.name)(qr[op.qubits[0]], qr[op.qubits[1]])
                else:
                    raise ValueError(f"unknown gate {op.name!r}")
            else:  # M
                if op.basis == "x":
                    qc.h(qr[op.qubit])
                elif op.basis == "y":
                    qc.sdg(qr[op.qubit]); qc.h(qr[op.qubit])
                qc.measure(qr[op.qubit], cregs[op.register][op.position])
        return qc

    def to_braket(self):
        from braket.circuits import Circuit
        c = Circuit()
        _BK = {"h": "h", "s": "s", "sdg": "si", "x": "x", "y": "y", "z": "z"}
        for op in self.ops:
            if isinstance(op, G):
                if op.name == "barrier":
                    continue
                elif op.name in _BK:
                    getattr(c, _BK[op.name])(op.qubits[0])
                elif op.name == "cx":
                    c.cnot(op.qubits[0], op.qubits[1])
                elif op.name == "cz":
                    c.cz(op.qubits[0], op.qubits[1])
                else:
                    raise ValueError(f"unknown gate {op.name!r}")
            else:  # M — Braket measures all used qubits implicitly at shots>0;
                   # only the basis-change rotation is emitted here.
                if op.basis == "x":
                    c.h(op.qubit)
                elif op.basis == "y":
                    c.si(op.qubit); c.h(op.qubit)
        return c


# ------------------------------------------------------------------- Steane
# Standard [[7,1,3]] encoder + stabilizer supports (mirror of
# nwt_substrate.heron.steane_pair_synthesis, kept as the single IR-level copy).
STEANE_PREP_H = (0, 1, 3)
STEANE_PREP_CX = [(0, 2), (0, 4), (0, 6), (1, 2), (1, 5), (1, 6), (3, 4), (3, 5), (3, 6)]
STEANE_STAB = [[0, 2, 4, 6], [1, 2, 5, 6], [3, 4, 5, 6]]
STEANE_LZ = [0, 1, 2]


def steane_base_ops(x_bits, z_bits) -> list:
    """|0_L> prep (3 H + 9 CX) then the walk's Pauli word on the 7 data qubits."""
    ops: list = [G("h", (q,)) for q in STEANE_PREP_H]
    ops += [G("cx", (c, t)) for c, t in STEANE_PREP_CX]
    for q in range(7):
        if x_bits[q] and z_bits[q]:
            ops.append(G("y", (q,)))
        elif x_bits[q]:
            ops.append(G("x", (q,)))
        elif z_bits[q]:
            ops.append(G("z", (q,)))
    return ops


def control_base_ops(kind: str) -> list:
    ops: list = [G("h", (q,)) for q in STEANE_PREP_H]
    ops += [G("cx", (c, t)) for c, t in STEANE_PREP_CX]
    if kind == "x7":
        ops += [G("x", (q,)) for q in range(7)]
    elif kind == "z7":
        ops += [G("z", (q,)) for q in range(7)]
    elif kind == "h7":
        ops += [G("h", (q,)) for q in range(7)]
    elif kind != "identity":
        raise ValueError(f"unknown control {kind!r}")
    return ops


# ----------------------------------------------- measurement-scheme strategies
def destructive_css(base_ops, name: str) -> tuple[CircuitSpec, CircuitSpec]:
    """Two 7-qubit circuits: Z-readout (Z-stabs + logical-Z) and X-readout
    (X-stabs). The portable, low-depth scheme (15 CZ on Garnet)."""
    z = CircuitSpec(7, list(base_ops) + [G("barrier")] + [M(q, "z", "c", q) for q in range(7)],
                    {"c": 7}, f"{name}_z")
    x = CircuitSpec(7, list(base_ops) + [G("barrier")] + [M(q, "x", "c", q) for q in range(7)],
                    {"c": 7}, f"{name}_x")
    return z, x


def ancilla_syndrome(base_ops, name: str) -> CircuitSpec:
    """Single 13-qubit circuit: 6 ancilla stabilizer extraction (c_syn) +
    destructive logical-Z readout of the data (c_lz). The original, depth-heavy
    scheme — fine on heavy-hex+DD, prohibitive on sparse grids."""
    ops = list(base_ops) + [G("barrier")]
    for k in range(3):                      # X-type stabilizers -> ancillas 7,8,9
        anc = 7 + k
        ops.append(G("h", (anc,)))
        ops += [G("cx", (anc, q)) for q in STEANE_STAB[k]]
        ops.append(G("h", (anc,)))
        ops.append(M(anc, "z", "c_syn", k))
    for k in range(3):                      # Z-type stabilizers -> ancillas 10,11,12
        anc = 10 + k
        ops += [G("cx", (q, anc)) for q in STEANE_STAB[k]]
        ops.append(M(anc, "z", "c_syn", 3 + k))
    ops += [M(q, "z", "c_lz", q) for q in range(7)]  # logical-Z destructive readout
    return CircuitSpec(13, ops, {"c_syn": 6, "c_lz": 7}, name)
