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
]
