# nwt_substrate.heron

> Quantum-circuit interface for substrate experiments on IBM Heron and other gate-model devices. The central object is the K_7 graph state `|K_7⟩` — the +1 eigenstate of the stabilizers `S_v = X_v ∏_{u≠v} Z_u` — which compiles to **7 qubits + 21 CZ gates** (one CZ per K_7 edge). Provides circuit builders, an **11-experiment registry** (6 run on hardware + 5 proposed), OPENQASM/script exporters, sidereal-geometry scheduling for the Experiment 11 anisotropy program, and Experiment 11 control circuits. Circuit construction requires the `[heron]` extra (`qiskit` / `qiskit-ibm-runtime`); the sidereal layer requires `astropy`.

[← Back to index](../index.md) · Source: [`nwt_substrate/heron/`](../../nwt_substrate/heron/) · Papers: [17](https://zenodo.org/records/15445103) (§13 K_7 stabilizer tests), [19](https://zenodo.org/communities/nwt) (W1-W7 quantum-substrate experiments + substrate-monism anisotropy probe), [20](https://zenodo.org/communities/nwt) (muon/neutron decay on the substrate vacuum), [22](https://zenodo.org/communities/nwt) (BH-cosmogenesis joint test)

## Common questions

### What is the K_7 graph state and how does it map to Heron?

`|K_7⟩` is the substrate's central object: the unique +1 eigenstate of every stabilizer `S_v = X_v ∏_{u≠v} Z_u`. As a graph state for the complete graph `K_7` it compiles directly to hardware — apply a Hadamard to each of the **7 = `N_VERTICES_K7`** qubits to make `|+⟩^{⊗7}`, then a CZ on every K_7 edge, **21 = `N_EDGES_K7`** of them. The substrate ISA → hardware closure is exact:

```python
import nwt_substrate.heron as heron
print(heron.substrate_breakdown())
# IBM Heron K_7 circuit from substrate ISA:
#
#     K_7 graph state preparation:
#       N_VERTICES_K7 = 7 system qubits
#       7 Hadamard gates → |+>^{⊗7}
#       N_EDGES_K7 = 21 CZ gates over k7_edge_list()
#       (every K_7 vertex pair receives one CZ)
#
#     Stabilizer measurements:
#       7 stabilizers S_v = X_v Π_{u≠v} Z_u
#       Each S_v has 1 X + (N_VERTICES_K7 - 1) = 7 Pauli factors
#       DEGREE_K7 = 6 (= N_VERTICES_K7 - 1, each vertex's K_7 valence)
#     ...
```

The same `21` that drives gravity's `α^(21/2)` Wilson exponent and chemistry's K_7-hub indicator literally compiles to 21 CZ gates on Heron. `k7_graph_state()` is depth 12; a full stabilizer self-test is depth ~82.

### What experiments are in the registry?

There are **11** registered experiments — **6 run** on real hardware (IDs 1–5, 10) and **5 proposed** (IDs 6–9, 11). `list_experiments()` returns the full list, `list_experiments("run")` / `list_experiments("proposed")` filter by status, and `summary()` gives the counts:

```python
heron.summary()
# {'total': 11,
#  'by_status': {'run': 6, 'proposed': 5},
#  'run_numbers': [1, 2, 3, 4, 5, 10],
#  'proposed_numbers': [6, 7, 8, 9, 11]}
```

The five run experiments (all on `ibm_marrakesh`) cover graph-state prep + stabilizer check, Z-basis entropy, X-basis tomography, the `⟨Y_u Y_v Y_w⟩ = 0` Fano null test, and a syndrome-distribution test; Experiment 10 is the parameterised muon-decay run from 2026-05-01. The proposed set runs from full tomography (6) and metrology (7) up to the headliners: two-universe slicing (8), "watch a neutron decay" (9), and the Experiment 11 sidereal anisotropy probe.

### Do I need qiskit?

Yes, for circuit construction. The circuit builders (`k7_graph_state`, `stabilizer_measurement`, `muon_decay_circuit`, …) and the control circuits raise `RuntimeError("qiskit not installed")` unless `qiskit` is importable; check `heron.HAS_QISKIT`. Install the `[heron]` extra (`pip install nwt-substrate[heron]`, which pulls `qiskit` + `qiskit-ibm-runtime`). The experiment **registry** (`HeronExperiment`, `EXPERIMENTS`, `list_experiments`, `experiment`, `summary`), the **exporters**, and the pure-numpy `directional_match_score` work without qiskit. The **sidereal geometry** layer additionally needs `astropy` (check `heron.HAS_ASTROPY`) for any function touching Earth orientation or LST.

### What is the sidereal / Experiment 11 anisotropy program?

Experiment 11 runs three identical K_7 stabilizer-measurement triplets at sidereal times `t`, `t + 12h` sidereal, and `t + 24h` sidereal (slots A / B / C). Slots A and C share lab-to-inertial orientation; B is rotated ~180° in the inertial frame. The drift-corrected signal is `σ_v = (B − A) − (C − A)/2` per stabilizer. The **substrate-monism prediction** is that if the `SO(8) → Spin(7) → 3+1` reduction picks an inertial-frame axis (the parent-BH spin axis = CMB "axis of evil", Paper 22), then B differs from A,C by a coherent per-vertex pattern. The `sidereal_geometry` module forward-models that pattern (`predicted_sigma_pattern`), schedules triplets at a fixed local sidereal time (`schedule_triplet_at_lst`), and scores observed-vs-predicted (`directional_match_score`). The `exp11_controls` module adds readout-floor and T1-floor control channels so any signal can be attributed to the substrate axis rather than thermal/readout drift.

### Are these simulations or direct observations?

The framing throughout the registry is that a Heron experiment is a substrate-algebra computation realized in qubits — substrate predictions about K_7 stabilizer expectations, entanglement bipartitions, syndrome distributions, or decay patterns become specific quantum-circuit runs whose outcomes are direct observations of substrate physics rather than analog simulations. Experiment 9 ("watch a neutron decay") and Experiment 10 (muon decay) carry this framing explicitly.

## Prediction / experiment table

The 11 registered experiments, from the real `EXPERIMENTS` registry (`status` ∈ {run, proposed}; run experiments use `ibm_marrakesh`):

| ID | Name | Status | Paper §  |
|---|---|---|---|
| 1 | K_7 graph-state preparation + stabilizer check | **run** | P17 §13 / P19 W1+W2 |
| 2 | Z-basis joint distribution + entanglement entropy | **run** | P19 W1+W2 |
| 3 | X-basis tomography | **run** | P19 W1+W2 |
| 4 | 3-body `⟨Y_u Y_v Y_w⟩ = 0` null test | **run** | P17 §13.3 |
| 5 | Syndrome distribution test | **run** | P17 §13.4 / P19 W3 |
| 6 | K_7 entanglement tomography | proposed | P19 W6 |
| 7 | Substrate metrology (α from K_7 Wilson amplitude) | proposed | P18 §X (TBD) |
| 8 | Two-universe slicing of `\|G_aug⟩` | proposed | P19 W7 / P22 |
| 9 | Watch a neutron decay | proposed | P20 (proposed) |
| 10 | Muon decay on K_7 vacuum (run 2026-05-01) | **run** | P20 / EW shim demo |
| 11 | Sidereal A/B/C substrate-anisotropy probe | proposed | P19 (anisotropy) / P22 |

Circuit-footprint facts (real `circuit_summary(...)` output, `qiskit` present):

| Builder | Qubits | Classical | Depth | Gates | Gate counts |
|---|---|---|---|---|---|
| `k7_graph_state()` | 7 | 0 | 12 | 28 | 7 H + 21 CZ |
| `stabilizer_measurement(0)` | 8 (+1 ancilla) | 1 | 10 | 10 | 2 H + 1 CX + 6 CZ + 1 meas |
| `k7_stabilizer_circuit(0)` | 7 | 7 | 13 | 36 | 8 H + 21 CZ + 7 meas |
| `full_k7_state_prep_with_measurement()` | 8 | 7 | 82 | 105 | 21 H + 63 CZ + 7 reset + 7 CX + 7 meas |
| `entanglement_tomography_x_basis()` | 7 | 7 | 14 | 42 | 14 H + 21 CZ + 7 meas |
| `y_basis_3body_correlator()` | 7 | 3 | 12 | 37 | 10 H + 21 CZ + 3 Sdg + 3 meas |
| `syndrome_distribution_circuit()` | 8 | 7 | 83 | 106 | 21 H + 63 CZ + 7 reset + 7 CX + 1 X + 7 meas |
| `muon_decay_circuit()` | 9 | 9 | 13 | 40 | 7 H + 21 CZ + 1 X + 1 RXX + 1 RYY + 9 meas |
| `readout_control_circuit()` | 7 | 7 | 1 | 7 | 7 meas |
| `t1_idle_control_circuit()` | 7 | 7 | 3 | 23 | 7 X + 7 delay + 7 meas + 2 barrier |

Every count above is a `pytest` assertion in [`nwt_substrate/tests/test_heron.py`](../../nwt_substrate/tests/test_heron.py); the sidereal layer is covered by [`test_heron_sidereal.py`](../../nwt_substrate/tests/test_heron_sidereal.py) and [`test_heron_sidereal_geometry.py`](../../nwt_substrate/tests/test_heron_sidereal_geometry.py) (44 tests total, all passing; astropy-gated tests skip cleanly without `astropy`).

## Quick start

```python
import nwt_substrate.heron as heron

heron.HAS_QISKIT      # → True  (False without the [heron] extra)
heron.HAS_ASTROPY     # → True  (False without astropy)

# Build the K_7 graph state preparation circuit (requires qiskit)
qc = heron.k7_graph_state()
heron.circuit_summary(qc)
# → {'n_qubits': 7, 'n_classical': 0, 'depth': 12, 'n_gates': 28,
#    'gate_counts': {'cz': 21, 'h': 7}}

# Verify substrate-correct counts (used by tests)
heron.verify_k7_circuit_substrate(qc)
# → {'n_qubits': 7, 'n_h': 7, 'n_cz': 21,
#    'n_edges_match': True, 'n_vertices_match': True}

# A single stabilizer-measurement circuit (CORRECT single-S_v architecture)
qc_sv = heron.k7_stabilizer_circuit(vertex=0)   # 7 qubits, 7 cbits, depth 13

# Inspect the experiment registry (works without qiskit)
heron.summary()
# → {'total': 11, 'by_status': {'run': 6, 'proposed': 5},
#    'run_numbers': [1, 2, 3, 4, 5, 10], 'proposed_numbers': [6, 7, 8, 9, 11]}

# Export a standalone, dependency-free runner for hardware
heron.export_experiment_script(1, "/tmp/exp01.py")   # → Path('/tmp/exp01.py')
```

## API by topic

### K_7 circuits (`circuits.py`)

| Function | Returns | Notes |
|---|---|---|
| `k7_graph_state(n_classical=0)` | `QuantumCircuit` | 7 H + 21 CZ; optional `n_classical` register |
| `stabilizer_measurement(vertex)` | `QuantumCircuit` | Ancilla-based `S_v` measurement (8 qubits) |
| `k7_stabilizer_circuit(vertex)` | `QuantumCircuit` | Prep + single basis-change + Z-measure all 7 (no mid-circuit reset) |
| `parse_k7_stabilizer_counts(counts)` | `(⟨S_v⟩, stderr)` | Parity-based stabilizer expectation from a counts dict |
| `full_k7_state_prep_with_measurement()` | `QuantumCircuit` | Prep + all 7 stabilizers via ancilla; self-test |
| `entanglement_tomography_x_basis()` | `QuantumCircuit` | Prep + X-basis readout of all 7 |
| `y_basis_3body_correlator(triple=(0,1,2))` | `QuantumCircuit` | Exp 4 `⟨Y_u Y_v Y_w⟩` Fano null test |
| `syndrome_distribution_circuit(error_qubit=0, error_axis="X")` | `QuantumCircuit` | Exp 5 syndrome after a single-qubit X/Y/Z error |
| `muon_decay_circuit(theta_param=None)` | `QuantumCircuit` | Exp 10 V-A XY-mixer between muon + electron ancillae |
| `circuit_summary(qc)` | `dict` | `n_qubits`, `n_classical`, `depth`, `n_gates`, `gate_counts` |
| `HAS_QISKIT` | `bool` | Whether qiskit is importable |

### Experiment registry (`experiments.py`)

| Symbol | Returns | Notes |
|---|---|---|
| `HeronExperiment` | dataclass | `number`, `name`, `status`, `description`, `backend`, `paper_section`, `notes` |
| `EXPERIMENTS` | `list[HeronExperiment]` | All 11 registered experiments |
| `list_experiments(status=None)` | `list[HeronExperiment]` | Optionally filter by `"run"` / `"proposed"` |
| `experiment(n)` | `HeronExperiment` | Look up by number (raises `KeyError` if absent) |
| `summary()` | `dict` | `total`, `by_status`, `run_numbers`, `proposed_numbers` |

### Exporters (`exporters.py`)

| Function | Returns | Notes |
|---|---|---|
| `export_experiment_script(n, output_path, backend="ibm_marrakesh", shots=10000)` | `Path` | Self-contained, dependency-free runner with the circuit builder inlined |
| `export_all_experiments(output_dir, backend="ibm_marrakesh", shots=10000)` | `list[Path]` | One script per registered experiment |

The exported scripts have **no** `nwt_substrate` runtime dependency: `pip install qiskit qiskit-ibm-runtime`, set `QISKIT_IBM_TOKEN`, and run. Without a token the script just prints the circuit structure.

### Sidereal geometry (`sidereal_geometry.py`, astropy-gated)

| Symbol | Returns | Notes |
|---|---|---|
| `Observatory` | frozen dataclass | `name`, `lat_deg`, `lon_deg`, `height_m`, `notes` |
| `YORKTOWN`, `EHNINGEN` | `Observatory` | IBM Yorktown Heights (hosts `ibm_marrakesh`) / IBM Ehningen (`ibm_aachen`) |
| `lab_to_icrs_matrix(unix_time, observatory)` | `(3,3)` ndarray | Lab Cartesian → ICRS rotation (needs astropy) |
| `lab_to_icrs(lab_vec, unix_time, observatory)` | `(3,)` ndarray | Single unit-vector transform |
| `default_k7_lab_positions(radius_m=0.005)` | `(7,3)` ndarray | Hexagon-plus-center K_7 lab embedding (pure math) |
| `predicted_sigma_pattern(t_A, t_B, t_C, observatory, …)` | `dict` | Forward-modelled `σ_v`; coupling `C_v = ((R p_v)·ŝ)²` |
| `predicted_sigma_pattern_asymmetric(…, epsilon=0.0)` | `dict` | Dipolar-jet refinement; `ε=0` recovers the symmetric model |
| `directional_match_score(observed, predicted, observed_std=None)` | `dict` | `pearson_r`, `cosine_similarity`, `best_amplitude`, `best_amplitude_chi2`, `dof` (pure numpy) |
| `lst_hours(unix_time, observatory)` | `float` | Local apparent sidereal time in hours [0, 24) |
| `next_lst_match_unix(target_lst_hours, observatory, …)` | `float` | Next unix time hitting a target LST |
| `schedule_triplet_at_lst(target_lst_hours, observatory, …)` | `dict` | A/B/C unix timestamps + their LSTs |
| `HAS_ASTROPY` | `bool` | Whether astropy is importable |

### Experiment 11 controls (`exp11_controls.py`)

| Function | Returns | Notes |
|---|---|---|
| `readout_control_circuit()` | `QuantumCircuit` | Prep `\|0⟩^7`, Z-measure all 7 (readout-floor control) |
| `t1_idle_control_circuit(idle_duration_ns=8000, use_delay_instruction=True)` | `QuantumCircuit` | Prep `\|1⟩^7`, idle, Z-measure all 7 (T1-floor control) |
| `parse_zbasis_expectations(counts)` | `list[float]` | Per-qubit `⟨Z_q⟩` from a 7-bit counts dict |
| `control_drift_corrected(A, B, C)` | `list[float]` | Same `(B−A) − (C−A)/2` correction as the K_7 channel |

### Substrate identities

| Symbol | Value / Returns | Notes |
|---|---|---|
| `substrate` | namespace | `N_VERTICES_K7`, `N_EDGES_K7`, `DEGREE_K7`, `DIM_OCTONION`, `k7_edge_list`, `breakdown`, `verify_k7_circuit` |
| `substrate_breakdown()` | `str` | Pretty-printed ISA → circuit closure (7 qubits, 21 CZ) |
| `verify_k7_circuit_substrate(qc)` | `dict` | Asserts a circuit has 7 H + 21 CZ on 7 qubits |
| `N_VERTICES_K7`, `N_EDGES_K7`, `DEGREE_K7`, `DIM_OCTONION` | `7`, `21`, `6`, `8` | Re-exported from `nwt_substrate.isa` |

## Worked examples

### Experiment registry lookup + the neutron-decay teaser

```python
import nwt_substrate.heron as heron

for e in heron.list_experiments():
    print(f"{e.number:2d} | {e.status:8s} | {e.name}")
#  1 | run      | K_7 graph-state preparation + stabilizer check
#  2 | run      | Z-basis joint distribution + entanglement entropy
#  3 | run      | X-basis tomography
#  4 | run      | 3-body <Y_u Y_v Y_w> = 0 null test
#  5 | run      | Syndrome distribution test
#  6 | proposed | K_7 entanglement tomography (proposed)
#  7 | proposed | Substrate metrology (proposed)
#  8 | proposed | Two-universe slicing of |G_aug> (proposed)
#  9 | proposed | Watch a neutron decay (proposed)
# 10 | run      | Muon decay on K_7 vacuum (run 2026-05-01)
# 11 | proposed | Sidereal A/B/C substrate-anisotropy probe (proposed)

e9 = heron.experiment(9)
print(e9.status, "|", e9.paper_section)   # proposed | Paper 20 (proposed)
```

### Substrate ISA → circuit closure

```python
from nwt_substrate.heron import substrate, k7_graph_state, verify_k7_circuit_substrate

substrate.N_VERTICES_K7, substrate.N_EDGES_K7, substrate.DEGREE_K7
# → (7, 21, 6)

qc = k7_graph_state()
verify_k7_circuit_substrate(qc)
# → {'n_qubits': 7, 'n_h': 7, 'n_cz': 21,
#    'n_edges_match': True, 'n_vertices_match': True}
```

The 21 CZ gates are exactly `len(substrate.k7_edge_list())` — one per K_7 edge.

### Scheduling and forward-modelling an Experiment 11 triplet (astropy)

```python
import numpy as np
import nwt_substrate.heron as heron

# Plan an A/B/C triplet anchored at local sidereal time 6.0 h at Yorktown
sched = heron.schedule_triplet_at_lst(6.0, heron.YORKTOWN, after_unix=1.764547200e9)
round(sched["A_lst"], 3), round(sched["B_lst"], 3), round(sched["C_lst"], 3)
# → (6.0, 18.0, 6.0)   # A and C share LST; B is +12 h sidereal

# Forward-model the per-vertex signal under the Axis-of-Evil hypothesis
patt = heron.predicted_sigma_pattern(
    sched["A_unix"], sched["B_unix"], sched["C_unix"], heron.YORKTOWN,
)
patt["predicted_axis_name"]        # → 'Axis of Evil (consensus)'
patt["predicted_axis_lb"]          # → (245.0, 61.5)   galactic (l, b)
np.round(patt["sigma_v_pred"], 4)
# → array([-0.0646, -0.0668, -0.0658, -0.0636, -0.0624, -0.0635, -0.0657])

# Score a (hypothetical) observed pattern against the prediction
score = heron.directional_match_score(
    observed_sigma_v=[0.10, 0.12, 0.11, 0.09, 0.08, 0.09, 0.11],
    predicted_sigma_v=[0.05, 0.06, 0.055, 0.045, 0.04, 0.045, 0.055],
)
score["pearson_r"], score["cosine_similarity"]   # → (1.0, 1.0)  (perfect match)
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `7` | K_7 vertex count = number of system qubits | `isa.N_VERTICES_K7` |
| `21` | K_7 edge count = number of CZ gates (one per edge) | `isa.N_EDGES_K7` |
| `6` | K_7 valence = `N_VERTICES_K7 − 1` = Pauli-`Z` factors per `S_v` | `isa.DEGREE_K7` |
| `8` | octonion / Spin(7) spinor dim (cross-shim ID) | `isa.DIM_OCTONION` |
| `S_v = X_v ∏_{u≠v} Z_u` | K_7 stabilizer generators (`|K_7⟩` is the +1 eigenstate) | `circuits.py`, `substrate_breakdown()` |
| `k7_edge_list()` | the explicit 21 edges driving the CZ layer | `isa.so7.k7_edge_list` |

`verify_k7_circuit_substrate` and the `test_heron.py` assertions refuse to pass if these counts drift — preventing silent divergence between `nwt_substrate.isa.constants` and the compiled circuits.

## Papers

- **Paper 17** §13 — K_7 stabilizer existence tests + the `⟨Y_u Y_v Y_w⟩ = 0` PSL(2,7) Fano-line null test (Experiments 1, 4, 5); also the `n = 21 → m_e` Wilson-amplitude metrology target (Experiment 7).
- **Paper 19** W1–W7 — quantum-substrate experiment series (graph-state prep, Z/X tomography, entanglement entropy, syndrome distribution, two-universe slicing) and the substrate-monism anisotropy probe (Experiment 11).
- **Paper 20** — muon decay (Experiment 10, run 2026-05-01 on `ibm_marrakesh`) and the proposed neutron-decay experiment (Experiment 9), via the substrate-algebraic `|M|²`.
- **Paper 22** — BH-cosmogenesis: the parent-BH spin axis sets the inertial-frame substrate axis tested by Experiment 11; cross-checked against the CMB "axis of evil" (galactic `l=245°, b=61.5°`).

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`cosmology`](cosmology.md) — anisotropy axes (`AXIS_OF_EVIL_CONSENSUS`, `HPA_CONSENSUS`) used by the Experiment 11 forward model
- [`qpu`](qpu.md) — device-agnostic QPU runner / spec / adapters layer
- [`topology`](topology.md) — K_7 graph, colored Jones, torus knots (substrate topology backing the carriers)
- [`isa`](../../nwt_substrate/isa/README.md) — `N_VERTICES_K7 = 7`, `N_EDGES_K7 = 21`, `DEGREE_K7 = 6`, `DIM_OCTONION = 8`
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — substrate-prediction benchmark suite
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
