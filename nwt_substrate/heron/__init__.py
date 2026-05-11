"""
nwt.heron
=========

Quantum-circuit interface for substrate experiments on IBM Heron and other
gate-model quantum devices.

The substrate's central object is the K_7 graph state |K_7> — the unique
+1 eigenstate of stabilizers S_v = X_v product_{u != v} Z_u.  Heron is
the natural hardware home for substrate experiments: 156 qubits, fast
two-qubit gates, calibrated for the depths needed (K_7 prep is depth ~7,
plus stabilizer measurements bring full sequences to depth ~50).

Quick start::

    import nwt_substrate.heron as heron

    # Build the K_7 preparation circuit
    qc = heron.k7_graph_state()
    print(qc.draw())

    # Build a single stabilizer measurement
    qc2 = heron.stabilizer_measurement(vertex=0)

    # Full self-test (prep + measure all 7 stabilizers)
    qc_test = heron.full_k7_state_prep_with_measurement()

    # X-basis tomography
    qc_x = heron.entanglement_tomography_x_basis()

    # Inspect the experiment registry
    heron.list_experiments()             # all 9
    heron.list_experiments("run")        # 5 already on hardware
    heron.list_experiments("proposed")   # 4 in the queue
    e9 = heron.experiment(9)             # the neutron decay teaser
"""

from .circuits import (
    k7_graph_state,
    stabilizer_measurement,
    k7_stabilizer_circuit,
    parse_k7_stabilizer_counts,
    full_k7_state_prep_with_measurement,
    entanglement_tomography_x_basis,
    y_basis_3body_correlator,
    syndrome_distribution_circuit,
    muon_decay_circuit,
    circuit_summary,
    HAS_QISKIT,
)

from .experiments import (
    HeronExperiment,
    EXPERIMENTS,
    list_experiments,
    experiment,
    summary,
)

from .exporters import (
    export_experiment_script,
    export_all_experiments,
)

from .sidereal_geometry import (
    Observatory,
    YORKTOWN,
    EHNINGEN,
    lab_to_icrs_matrix,
    lab_to_icrs,
    default_k7_lab_positions,
    predicted_sigma_pattern,
    directional_match_score,
    lst_hours,
    next_lst_match_unix,
    schedule_triplet_at_lst,
)

from .exp11_controls import (
    readout_control_circuit,
    t1_idle_control_circuit,
    parse_zbasis_expectations,
    control_drift_corrected,
)

# ---- Substrate-ISA identities ----
from ..isa import (
    N_VERTICES_K7,
    N_EDGES_K7,
    DEGREE_K7,
    DIM_OCTONION,
    k7_edge_list,
)


def substrate_breakdown() -> str:
    """Pretty-print the heron-circuit substrate identities.

    The K_7 graph state circuit on Heron is the most direct hardware
    realization of the substrate algebra: N_VERTICES_K7 = 7 qubits,
    N_EDGES_K7 = 21 CZ gates, exactly one per K_7 edge.

    This is the substrate ISA → quantum hardware closure: the same
    21 that drives gravity's α^(21/2) Wilson exponent and chemistry's
    K_7-hub indicator literally compiles to 21 CZ gates on Heron.
    """
    lines = ["IBM Heron K_7 circuit from substrate ISA:"]
    lines.append("")
    lines.append(f"    K_7 graph state preparation:")
    lines.append(f"      N_VERTICES_K7 = {N_VERTICES_K7} system qubits")
    lines.append(f"      {N_VERTICES_K7} Hadamard gates → |+>^{{⊗{N_VERTICES_K7}}}")
    lines.append(f"      N_EDGES_K7 = {N_EDGES_K7} CZ gates over k7_edge_list()")
    lines.append(f"      (every K_7 vertex pair receives one CZ)")
    lines.append("")
    lines.append(f"    Stabilizer measurements:")
    lines.append(f"      {N_VERTICES_K7} stabilizers S_v = X_v Π_{{u≠v}} Z_u")
    lines.append(f"      Each S_v has 1 X + (N_VERTICES_K7 - 1) = "
                 f"{1 + (N_VERTICES_K7 - 1)} Pauli factors")
    lines.append(f"      DEGREE_K7 = {DEGREE_K7} (= N_VERTICES_K7 - 1, each")
    lines.append(f"        vertex's K_7 valence)")
    lines.append("")
    lines.append(f"    Circuit footprint:")
    lines.append(f"      k7_graph_state():  {N_VERTICES_K7} qubits, "
                 f"{N_VERTICES_K7} H + {N_EDGES_K7} CZ gates")
    lines.append(f"      stabilizer_measurement(v): {N_VERTICES_K7 + 1} qubits "
                 f"(system + 1 ancilla)")
    lines.append(f"      muon_decay_circuit():  {N_VERTICES_K7 + 2} qubits "
                 f"(K_7 + muon + electron)")
    lines.append("")
    lines.append(f"    Cross-shim ID — same {N_EDGES_K7} edges of K_7 drive:")
    lines.append(f"      chemistry's K_7-hub Tr(M²) detection (coronene)")
    lines.append(f"      gravity's α^({N_EDGES_K7}/2) Wilson exponent")
    lines.append(f"      qed's so(7) adjoint that holds the γ-matrix algebra")
    lines.append(f"      qcd's so(7) adjoint containing SU(3) ⊂ G_2 ⊂ Spin(7)")
    lines.append(f"      particles' {N_VERTICES_K7} carrier-knot types")
    lines.append(f"      ew's so(7) generator decomposition 21 = 6 + 3 + 12")
    lines.append(f"      heron's CZ-gate count on real quantum hardware")
    lines.append("")
    lines.append(f"    The K_7 substrate is now active-encoded all the way")
    lines.append(f"    from substrate algebra to qiskit quantum circuits.")
    return "\n".join(lines)


def verify_k7_circuit_substrate(qc) -> dict:
    """Verify that a K_7 circuit has the substrate-correct counts.

    Counts CZ gates and H gates and returns dict with:
      n_qubits, n_h, n_cz, n_edges_match, n_vertices_match

    Used by tests to assert that k7_graph_state() has exactly
    N_EDGES_K7 = 21 CZs and N_VERTICES_K7 = 7 H's.
    """
    if not HAS_QISKIT:
        raise RuntimeError("qiskit not installed")
    n_h = sum(1 for instr in qc.data if instr.operation.name == "h")
    n_cz = sum(1 for instr in qc.data if instr.operation.name == "cz")
    return {
        "n_qubits": qc.num_qubits,
        "n_h": n_h,
        "n_cz": n_cz,
        "n_edges_match": (n_cz == N_EDGES_K7),
        "n_vertices_match": (qc.num_qubits == N_VERTICES_K7),
    }


class _SubstrateNamespace:
    """Substrate identities visible from the Heron shim."""
    N_VERTICES_K7 = N_VERTICES_K7
    N_EDGES_K7 = N_EDGES_K7
    DEGREE_K7 = DEGREE_K7
    DIM_OCTONION = DIM_OCTONION
    k7_edge_list = staticmethod(k7_edge_list)
    breakdown = staticmethod(substrate_breakdown)
    verify_k7_circuit = staticmethod(verify_k7_circuit_substrate)


substrate = _SubstrateNamespace()


__all__ = [
    # circuits
    "k7_graph_state",
    "stabilizer_measurement",
    "k7_stabilizer_circuit",
    "parse_k7_stabilizer_counts",
    "full_k7_state_prep_with_measurement",
    "entanglement_tomography_x_basis",
    "y_basis_3body_correlator",
    "syndrome_distribution_circuit",
    "muon_decay_circuit",
    "circuit_summary",
    "HAS_QISKIT",
    # experiments registry
    "HeronExperiment",
    "EXPERIMENTS",
    "list_experiments",
    "experiment",
    "summary",
    # exporters
    "export_experiment_script",
    "export_all_experiments",
    # sidereal geometry (Exp 11 directional layer)
    "Observatory",
    "YORKTOWN",
    "EHNINGEN",
    "lab_to_icrs_matrix",
    "lab_to_icrs",
    "default_k7_lab_positions",
    "predicted_sigma_pattern",
    "directional_match_score",
    "lst_hours",
    "next_lst_match_unix",
    "schedule_triplet_at_lst",
    # Exp 11 control channels (drift attribution)
    "readout_control_circuit",
    "t1_idle_control_circuit",
    "parse_zbasis_expectations",
    "control_drift_corrected",
    # substrate identities
    "substrate", "substrate_breakdown", "verify_k7_circuit_substrate",
]
