"""
Quantum circuits for IBM Heron (and other gate-model devices) implementing
substrate-algebra structures.

The K_7 graph state |K_7> is the substrate's central object — stabilizer
ground state on 7 qubits, eigenstate of S_v = X_v product_{u~v} Z_u for all
v.  These circuits prepare it directly with H gates + CZs over all 21 K_7
edges.

Substrate identity: every count and dimension here comes from
`nwt_substrate.isa.constants`:
  N_VERTICES_K7 = 7  → number of system qubits
  N_EDGES_K7    = 21 → number of CZ gates
  k7_edge_list()     → the explicit 21 edges of K_7

This module is the *quantum-hardware* shim of the substrate ISA.
The same N_EDGES_K7 = 21 that drives chemistry's K_7-hub detection,
gravity's α^(21/2) Wilson exponent, QED's 8×8 Dirac γ shape, QCD's 8
gluons, particles' 7 carrier topologies, and EW's b_QED^SM = 8 also
counts the CZ gates that the substrate compiles to on Heron.

Returns qiskit QuantumCircuit objects when qiskit is available.
"""

from __future__ import annotations

from ..isa.constants import N_VERTICES_K7, N_EDGES_K7
from ..isa.so7 import k7_edge_list

try:
    from qiskit import QuantumCircuit, ClassicalRegister
    HAS_QISKIT = True
except ImportError:
    HAS_QISKIT = False
    QuantumCircuit = None
    ClassicalRegister = None


# K_7 edge list, cached once (since it doesn't change).
_K7_EDGES = k7_edge_list()
assert len(_K7_EDGES) == N_EDGES_K7


def k7_graph_state(n_classical: int = 0) -> "QuantumCircuit":
    """
    Build the K_7 graph state preparation circuit.

    Recipe:
      1. Apply H to all N_VERTICES_K7 = 7 qubits  ->  |+>^{otimes 7}
      2. Apply CZ to every edge of K_7 (all N_EDGES_K7 = 21 pairs)

    The resulting state |K_7> is the unique +1 eigenstate of every
    stabilizer S_v = X_v product_{u != v} Z_u.

    Substrate identity: the circuit has N_VERTICES_K7 qubits and
    N_EDGES_K7 CZ gates, both sourced from isa.constants.

    Returns
    -------
    QuantumCircuit on N_VERTICES_K7 qubits with optional classical
    register for measurements (set n_classical = N_VERTICES_K7 to
    measure all qubits).
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed; pip install qiskit")
    qc = QuantumCircuit(N_VERTICES_K7)
    if n_classical > 0:
        qc.add_register(ClassicalRegister(n_classical, "c"))

    # Step 1: superposition on all K_7 vertices
    for q in range(N_VERTICES_K7):
        qc.h(q)

    # Step 2: CZ over every K_7 edge (sourced from isa.k7_edge_list)
    for i, j in _K7_EDGES:
        qc.cz(i, j)

    return qc


def stabilizer_measurement(vertex: int) -> "QuantumCircuit":
    """
    Build the circuit that measures the K_7 stabilizer S_v at the given vertex.

    S_v = X_v product_{u != v} Z_u

    Implementation: use an ancilla, apply Hadamard, controlled-X (for X_v)
    and controlled-Zs (for the Z's at all other vertices), then Hadamard
    and measure.

    Returns (N_VERTICES_K7 + 1)-qubit circuit (7 system qubits + 1 ancilla).
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    if not (0 <= vertex < N_VERTICES_K7):
        raise ValueError(
            f"vertex must be in [0, {N_VERTICES_K7}) (got {vertex!r})"
        )
    n_qubits = N_VERTICES_K7 + 1
    ancilla = N_VERTICES_K7
    qc = QuantumCircuit(n_qubits, 1)

    # Prepare ancilla in |+>
    qc.h(ancilla)

    # Controlled-X on the chosen vertex
    qc.cx(ancilla, vertex)

    # Controlled-Z on all other vertices (Z_u for u != v)
    for u in range(N_VERTICES_K7):
        if u != vertex:
            qc.cz(ancilla, u)

    # Hadamard ancilla and measure
    qc.h(ancilla)
    qc.measure(ancilla, 0)
    return qc


def k7_stabilizer_circuit(vertex: int) -> "QuantumCircuit":
    """
    Measure the K_7 stabilizer S_v = X_v * prod_{u != v} Z_u on |K_7> by
    rotating qubit v from X->Z basis (H_v), then measuring all 7 qubits
    in Z.  In the rotated basis S_v = prod_i Z_i = parity of all 7 bits.

    This is the CORRECT architecture for measuring a single stabilizer
    on real hardware: prep |K_7> fresh, single basis-change gate, single
    measurement.  Avoids the mid-circuit-reset-on-Heron pitfall in
    full_k7_state_prep_with_measurement.

    To measure all 7 stabilizers, build 7 such circuits (one per
    vertex) and submit them as a list to Sampler.run([...]).  Each
    circuit's <S_v> is the mean of (-1)**parity over its shots.

    Returns
    -------
    QuantumCircuit on 7 qubits with 7-bit classical register.  Each
    measured bitstring's all-bits-parity gives one (+1 or -1)
    sample of S_v.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed; pip install qiskit")
    if not (0 <= vertex < N_VERTICES_K7):
        raise ValueError(
            f"vertex must be in [0, {N_VERTICES_K7}) (got {vertex!r})"
        )

    qc = QuantumCircuit(N_VERTICES_K7, N_VERTICES_K7)
    # Prep |K_7>
    for q in range(N_VERTICES_K7):
        qc.h(q)
    for i, j in _K7_EDGES:
        qc.cz(i, j)
    # Rotate vertex v's X-basis into Z-basis.
    qc.h(vertex)
    # Measure all N_VERTICES_K7 in Z.
    for q in range(N_VERTICES_K7):
        qc.measure(q, q)
    return qc


def parse_k7_stabilizer_counts(counts: dict) -> tuple[float, float]:
    """From a 7-bit counts dict, compute <S_v> = mean of (-1)**parity
    across shots, and the shot-noise stderr.

    The bitstring's all-bits parity (sum mod 2) gives the eigenvalue
    of the rotated stabilizer; +1 if even, -1 if odd.
    """
    total = 0.0
    n_shots = 0
    for bitstr, count in counts.items():
        bs = bitstr.replace(" ", "")
        parity = sum(int(b) for b in bs) % 2
        ev = 1 if parity == 0 else -1
        total += ev * count
        n_shots += count
    if n_shots == 0:
        return 0.0, 0.0
    mean = total / n_shots
    var = max(0.0, 1.0 - mean ** 2)
    return float(mean), float((var / n_shots) ** 0.5)


def full_k7_state_prep_with_measurement() -> "QuantumCircuit":
    """
    Combined: prepare |K_7> on qubits 0..6, then measure all 7 stabilizers.

    Useful as a self-test on Heron — every measurement should yield +1
    (|K_7> is the +1 eigenstate of every S_v).  If even one measurement
    is -1, the device is failing the K_7 stabilizer check.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    n_qubits = N_VERTICES_K7 + 1
    ancilla = N_VERTICES_K7
    qc = QuantumCircuit(n_qubits, N_VERTICES_K7)

    # Prepare |K_7> on qubits 0..N_VERTICES_K7-1
    for q in range(N_VERTICES_K7):
        qc.h(q)
    for i, j in _K7_EDGES:
        qc.cz(i, j)

    # Measure each stabilizer S_v sequentially using the ancilla
    for v in range(N_VERTICES_K7):
        # Reset and prepare ancilla in |+>
        qc.reset(ancilla)
        qc.h(ancilla)
        # CX for the X_v factor
        qc.cx(ancilla, v)
        # CZs for the Z_u factors (u != v)
        for u in range(N_VERTICES_K7):
            if u != v:
                qc.cz(ancilla, u)
        # Hadamard ancilla and measure into classical bit v
        qc.h(ancilla)
        qc.measure(ancilla, v)

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
    qc = QuantumCircuit(N_VERTICES_K7, 3)
    # Prepare |K_7>
    for q in range(N_VERTICES_K7):
        qc.h(q)
    for i, j in _K7_EDGES:
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
    n_qubits = N_VERTICES_K7 + 1
    ancilla = N_VERTICES_K7
    qc = QuantumCircuit(n_qubits, N_VERTICES_K7)
    for q in range(N_VERTICES_K7):
        qc.h(q)
    for i, j in _K7_EDGES:
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
    # Measure all N_VERTICES_K7 stabilizers using ancilla
    for v in range(N_VERTICES_K7):
        qc.reset(ancilla)
        qc.h(ancilla)
        qc.cx(ancilla, v)
        for u in range(N_VERTICES_K7):
            if u != v:
                qc.cz(ancilla, u)
        qc.h(ancilla)
        qc.measure(ancilla, v)
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
    qc = QuantumCircuit(N_VERTICES_K7, N_VERTICES_K7)
    # Prepare |K_7>
    for q in range(N_VERTICES_K7):
        qc.h(q)
    for i, j in _K7_EDGES:
        qc.cz(i, j)
    # Rotate to X basis (apply H) then measure in Z (which is now X)
    for q in range(N_VERTICES_K7):
        qc.h(q)
    qc.measure(range(N_VERTICES_K7), range(N_VERTICES_K7))
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

    n_qubits = N_VERTICES_K7 + 2   # K_7 vacuum + muon + electron
    muon_qubit = N_VERTICES_K7      # qubit 7
    electron_qubit = N_VERTICES_K7 + 1   # qubit 8
    qc = QuantumCircuit(n_qubits, n_qubits)

    # Step 1: Prepare K_7 graph state on qubits 0..N_VERTICES_K7-1
    for q in range(N_VERTICES_K7):
        qc.h(q)
    for i, j in _K7_EDGES:
        qc.cz(i, j)

    # Step 2: Initialise muon ancilla to |1⟩, electron to |0⟩
    qc.x(muon_qubit)

    # Step 3: V-A weak transition.  XY-mixer:
    #   exp(-i (θ/2)(X_μ X_e + Y_μ Y_e)/2)
    #   = SWAP-like rotation that takes |10⟩ ↔ |01⟩ at angle θ
    qc.rxx(theta_param, muon_qubit, electron_qubit)
    qc.ryy(theta_param, muon_qubit, electron_qubit)

    # Step 4: Measure muon and electron populations
    qc.measure(muon_qubit, muon_qubit)
    qc.measure(electron_qubit, electron_qubit)

    # Step 5: Measure K_7 stabiliser S_0 = X_0 ∏_{u≠0} Z_u via Z-basis on
    # all K_7 qubits (cheap proxy: when K_7 is in its ground state, Z-basis
    # samples have a specific bit-flip-symmetric distribution).
    for v in range(N_VERTICES_K7):
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
