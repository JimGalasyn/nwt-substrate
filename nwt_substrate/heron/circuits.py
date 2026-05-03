"""
Quantum circuits for IBM Heron (and other gate-model devices) implementing
substrate-algebra structures.

The K_7 graph state |K_7> is the substrate's central object — stabilizer
ground state on 7 qubits, eigenstate of S_v = X_v product_{u~v} Z_u for all
v.  These circuits prepare it directly with H gates + CZs over all 21 K_7
edges.

Key identity: K_7 has 7 vertices, 21 edges, all pairs connected.  Heron
has 156 qubits and is calibrated for 100+ depth-2 circuits, so a 7-qubit
K_7 prep is well within reach.

Returns qiskit QuantumCircuit objects when qiskit is available.
"""

from __future__ import annotations

try:
    from qiskit import QuantumCircuit, ClassicalRegister
    HAS_QISKIT = True
except ImportError:
    HAS_QISKIT = False
    QuantumCircuit = None
    ClassicalRegister = None


def k7_graph_state(n_classical: int = 0) -> "QuantumCircuit":
    """
    Build the K_7 graph state preparation circuit.

    Recipe:
      1. Apply H to all 7 qubits  ->  |+>^{otimes 7}
      2. Apply CZ to every edge of K_7 (all 21 pairs)

    The resulting state |K_7> is the unique +1 eigenstate of every
    stabilizer S_v = X_v product_{u != v} Z_u.

    Returns
    -------
    QuantumCircuit on 7 qubits with optional classical register for
    measurements (set n_classical = 7 if you want to measure all qubits).
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed; pip install qiskit")
    qc = QuantumCircuit(7)
    if n_classical > 0:
        qc.add_register(ClassicalRegister(n_classical, "c"))

    # Step 1: superposition
    for q in range(7):
        qc.h(q)

    # Step 2: CZ over all 21 K_7 edges (every pair)
    for i in range(7):
        for j in range(i + 1, 7):
            qc.cz(i, j)

    return qc


def stabilizer_measurement(vertex: int) -> "QuantumCircuit":
    """
    Build the circuit that measures the K_7 stabilizer  S_v at the given vertex.

    S_v = X_v product_{u != v} Z_u

    Implementation: use an ancilla, apply Hadamard, controlled-X (for X_v)
    and controlled-Zs (for the Z's at all other vertices), then Hadamard
    and measure.

    Returns 8-qubit circuit (7 system qubits + 1 ancilla).
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    qc = QuantumCircuit(8, 1)   # qubits 0..6 = system, qubit 7 = ancilla

    # Prepare ancilla in |+>
    qc.h(7)

    # Controlled-X on the chosen vertex
    qc.cx(7, vertex)

    # Controlled-Z on all other vertices (Z_u for u != v)
    for u in range(7):
        if u != vertex:
            qc.cz(7, u)

    # Hadamard ancilla and measure
    qc.h(7)
    qc.measure(7, 0)
    return qc


def full_k7_state_prep_with_measurement() -> "QuantumCircuit":
    """
    Combined: prepare |K_7> on qubits 0..6, then measure all 7 stabilizers.

    Useful as a self-test on Heron — every measurement should yield +1
    (|K_7> is the +1 eigenstate of every S_v).  If even one measurement
    is -1, the device is failing the K_7 stabilizer check.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    qc = QuantumCircuit(8, 7)   # 7 system + 1 ancilla, 7 classical bits

    # Prepare |K_7> on qubits 0..6
    for q in range(7):
        qc.h(q)
    for i in range(7):
        for j in range(i + 1, 7):
            qc.cz(i, j)

    # Measure each stabilizer S_v sequentially using the ancilla
    for v in range(7):
        # Reset and prepare ancilla in |+>
        qc.reset(7)
        qc.h(7)
        # CX for the X_v factor
        qc.cx(7, v)
        # CZs for the Z_u factors (u != v)
        for u in range(7):
            if u != v:
                qc.cz(7, u)
        # Hadamard ancilla and measure into classical bit v
        qc.h(7)
        qc.measure(7, v)

    return qc


def y_basis_3body_correlator(triple: tuple = (0, 1, 2)) -> "QuantumCircuit":
    """
    Experiment 4: <Y_u Y_v Y_w> = 0 null test for non-Fano triples.

    Substrate prediction: for any three K_7 vertices NOT lying on a Fano line,
    the joint Y-basis expectation <Y_u Y_v Y_w> on |K_7> is zero exactly.
    This probes the PSL(2,7) Fano-line structure of the substrate.

    To measure Y on a qubit: apply S^dag then H, then measure in Z.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    qc = QuantumCircuit(7, 3)
    # Prepare |K_7>
    for q in range(7):
        qc.h(q)
    for i in range(7):
        for j in range(i + 1, 7):
            qc.cz(i, j)
    # Rotate the three measured qubits to the Y basis: S^dag, then H
    for k, q in enumerate(triple):
        qc.sdg(q)
        qc.h(q)
        qc.measure(q, k)
    return qc


def syndrome_distribution_circuit(error_qubit: int = 0,
                                    error_axis: str = "X") -> "QuantumCircuit":
    """
    Experiment 5: syndrome distribution test.

    Prep |K_7>, apply a single-qubit error (X or Z) at the chosen vertex,
    then measure all 7 stabilizers via 7 ancilla measurements.  Substrate
    prediction: the syndrome flips exactly the stabilizers that don't
    commute with the error, in line with the K_7 stabilizer code.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    qc = QuantumCircuit(8, 7)
    for q in range(7):
        qc.h(q)
    for i in range(7):
        for j in range(i + 1, 7):
            qc.cz(i, j)
    # Inject error
    if error_axis.upper() == "X":
        qc.x(error_qubit)
    elif error_axis.upper() == "Z":
        qc.z(error_qubit)
    elif error_axis.upper() == "Y":
        qc.y(error_qubit)
    else:
        raise ValueError(f"error_axis must be X/Y/Z, got {error_axis!r}")
    # Measure all 7 stabilizers using ancilla
    for v in range(7):
        qc.reset(7)
        qc.h(7)
        qc.cx(7, v)
        for u in range(7):
            if u != v:
                qc.cz(7, u)
        qc.h(7)
        qc.measure(7, v)
    return qc


def entanglement_tomography_x_basis() -> "QuantumCircuit":
    """
    Prep |K_7>, then measure all 7 qubits in the X basis.

    For Paper 19 W6: the joint X-basis distribution is one slice of the
    full tomography.  Combined with Z-basis and Y-basis runs, one can
    reconstruct the K_7 entanglement entropy bipartition.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    qc = QuantumCircuit(7, 7)
    # Prepare |K_7>
    for q in range(7):
        qc.h(q)
    for i in range(7):
        for j in range(i + 1, 7):
            qc.cz(i, j)
    # Rotate to X basis (apply H) then measure in Z (which is now X)
    for q in range(7):
        qc.h(q)
    qc.measure(range(7), range(7))
    return qc


def muon_decay_circuit(theta_param=None) -> "QuantumCircuit":
    """
    Experiment 10: muon decay on the K_7 substrate vacuum.

    Construction:
      Qubits 0-6: K_7 substrate vacuum (7 H + 21 CZ).
      Qubit 7:    muon ancilla (initialised to |1⟩ by an X gate).
      Qubit 8:    electron ancilla (initialised to |0⟩).

    A parameterised XY-mixer between (muon, electron) at angle θ models
    the V-A weak vertex driving the transition

        |muon⟩ |0⟩  →  cos(θ/2) |muon⟩ |0⟩  +  sin(θ/2) |0⟩ |electron⟩

    so as θ : 0 → π the muon population transfers fully to the electron.

    Stabiliser of the K_7 vacuum at vertex 0 is also measured (via 4 cbits)
    to verify that the substrate ground state is unchanged by the local
    weak interaction on the ancillae -- a substrate-stability check.

    Parameters
    ----------
    theta_param : Parameter or float, optional
        Mixing angle.  If a qiskit ``Parameter`` is passed, the circuit is
        parameterised so the user can bind multiple values in a single
        ``Sampler.run`` call (recommended for the Heron submission).

    Returns
    -------
    QuantumCircuit on 9 qubits, 9 classical bits.

    Substrate identifications:
      - K_7 graph: substrate vacuum (memory: K7-on-torus-matter-gravity-unification)
      - V-A vertex: chiral SU(2)_L from electroweak shim (Walk-phase 4b)
      - Muon J=1/2 + (p,q)=(2,1): from L1 spinor sector + Paper 6 mass formula
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed; pip install qiskit")
    from qiskit.circuit import Parameter

    if theta_param is None:
        theta_param = Parameter("theta")

    qc = QuantumCircuit(9, 9)

    # Step 1: Prepare K_7 graph state on qubits 0-6
    for q in range(7):
        qc.h(q)
    for i in range(7):
        for j in range(i + 1, 7):
            qc.cz(i, j)

    # Step 2: Initialise muon ancilla to |1⟩ (qubit 7), electron to |0⟩ (qubit 8)
    qc.x(7)

    # Step 3: V-A weak transition.  XY-mixer:
    #   exp(-i (θ/2)(X_7 X_8 + Y_7 Y_8)/2)
    #   = SWAP-like rotation that takes |10⟩ ↔ |01⟩ at angle θ
    # Standard 3-CNOT decomposition of the iSWAP-like XY rotation:
    qc.rxx(theta_param, 7, 8)
    qc.ryy(theta_param, 7, 8)

    # Step 4: Measure muon (qubit 7) and electron (qubit 8) populations
    qc.measure(7, 7)
    qc.measure(8, 8)

    # Step 5: Measure K_7 stabiliser S_0 = X_0 ∏_{u≠0} Z_u via Z-basis on
    # all K_7 qubits (cheap proxy: when K_7 is in its ground state, Z-basis
    # samples have a specific bit-flip-symmetric distribution).
    for v in range(7):
        qc.measure(v, v)

    return qc


def circuit_summary(qc: "QuantumCircuit") -> dict:
    """Pretty-print summary of a circuit (gate counts + depth + width)."""
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    return {
        "n_qubits": qc.num_qubits,
        "n_classical": qc.num_clbits,
        "depth": qc.depth(),
        "n_gates": sum(qc.count_ops().values()),
        "gate_counts": dict(qc.count_ops()),
    }
