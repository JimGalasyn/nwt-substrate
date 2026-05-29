# nwt_substrate.qpu

> A vendor-neutral QPU interface that separates substrate physics (`spec`, `decode`) from vendor plumbing (adapters for IBM, AWS Braket, and a local simulator). The canonical `Counts` contract normalizes measurement results across vendors — most importantly bit-ordering — so that the *same* substrate experiment decodes identically regardless of backend. This is an **infrastructure** shim: it carries the K_7 / SU(2)_5 experiments onto real hardware rather than predicting a physics observable. M1 milestone: `spec` + `decode` + adapters with the canonical counts contract.

[← Back to index](../index.md) · Source: [`nwt_substrate/qpu/`](../../nwt_substrate/qpu/) · Design: `analysis/qpu_interface_design.md` (null-worldtube-private) · [NWT community](https://zenodo.org/communities/nwt)

## Common questions

### What problem does the qpu shim solve?

Substrate experiments (Steane `[[7,1,3]]` walk-encodings, K_7 stabilizer reads) had ~26 hand-rolled decoders scattered across `analysis/`, each re-interpreting raw counts and each re-tripping the same vendor bit-ordering footgun. The `qpu` shim collapses that into one path: a tiny vendor-neutral circuit IR (`spec`), one decoder (`decode`), and per-vendor adapters whose *only* job is to compile, submit, and normalize raw counts into a canonical `Counts`. After normalization, all physics is vendor-agnostic.

### What is the canonical `Counts` contract?

`Counts` is a frozen dataclass with two fields: `register_order: tuple[str, ...]` (classical registers in a fixed order) and `data: dict[tuple[str, ...], int]` (a tuple of per-register bitstrings → shot count). The single load-bearing convention is:

> Within each register's bitstring, **character index `i` is qubit/cbit `i` — index 0 is the LEFTMOST character.**

Every adapter is responsible for normalizing its SDK's raw output to this convention, so the order is interpreted *exactly once*:

- **Braket** keys (`measurement_counts`) are already in this order (leftmost = qubit 0), so `canonicalize_braket` is a pure slice into named registers using the circuit's measured layout.
- **Qiskit** `SamplerV2` returns little-endian per-register bitstrings (rightmost char = cbit 0), so `canonicalize_ibm` reverses each register's string (`[::-1]`) and zips them per shot.

Because bit order is reversed in *one place* (`canonicalize_ibm`), `decode` cannot be vendor-sensitive. The contract is asserted directly in `test_qpu.py::test_canonical_counts_vendor_equivalence`: a Braket-style raw key and the matching reversed Qiskit-style string normalize to the *same* `Counts` and decode to the same Fano point.

### Which backends are supported?

Three, all behind the same `Backend`/`Job` protocol:

| Backend | SDK | Notes |
|---|---|---|
| `IBMBackend` | `qiskit-ibm-runtime` | `SamplerV2`; one batched job-id per submission |
| `BraketBackend` | `amazon-braket-sdk` | IQM, Rigetti, AQT devices via ARN; one ARN per circuit |
| `SimulatorBackend` | `braket` LocalSimulator | free dry-runs; decodes through the *same* `canonicalize_braket` path as real Braket |

### How does it relate to the `heron` shim?

`heron` *builds* the K_7 graph-state and Steane circuits for IBM Heron (it owns the substrate circuit constructions). `qpu` is the lower, vendor-neutral run/decode layer: it takes a `CircuitSpec`, emits Qiskit *or* Braket from it, submits to any of the three backends, and decodes the normalized counts. `qpu` is where IBM stops being special — the same spec runs on Braket/IQM and the local simulator without touching the physics.

### Do I need vendor SDKs installed?

Only for the adapter you actually use. `import nwt_substrate.qpu` pulls in `spec` and `decode` with **zero** SDK dependencies — `to_qiskit`/`to_braket` and every adapter import their SDK lazily. `canonicalize_braket` and `canonicalize_ibm` are pure functions (SDK-free, unit-tested against fixtures), so you can normalize and decode saved counts without any SDK. Constructing a `BraketBackend` needs `amazon-braket-sdk`; `IBMBackend` needs `qiskit-ibm-runtime`; `SimulatorBackend` needs Braket's `LocalSimulator`.

## Architecture

Three layers, with the canonical `Counts` as the seam between vendor plumbing and substrate physics:

| Layer | Module | Responsibility |
|---|---|---|
| **spec** (vendor-neutral) | `qpu.spec` | Tiny circuit IR (`G`, `M`, `CircuitSpec`) + Steane `[[7,1,3]]` builders + measurement-scheme strategies. Emits Qiskit or Braket on demand. |
| **adapters** (vendor plumbing) | `qpu.adapters.{ibm,braket,simulator}` | Per-vendor compile + submit + **canonicalize raw counts → `Counts`**. |
| **decode** (vendor-neutral) | `qpu.decode` | The single Steane decoder: Fano-point syndrome + logical-Z + pass/fail verdicts, consuming only canonical `Counts`. |

Two more modules orchestrate around these (later milestones): `qpu.capabilities` (`preflight()` routing/cost/window guardrail) and `qpu.runner` (`run()` lifecycle: build → preflight/stage → submit → poll → decode → JSON, with on-disk resume handles). `qpu.experiments` defines the exp12 Phase-2 walk set against the interface.

### Dataclass / protocol contract (`adapters/base.py`)

| Type | Kind | Fields / methods |
|---|---|---|
| `Counts` | `@dataclass(frozen=True)` | `register_order`, `data`; `total()`, `register(name)`, `single()`, `from_canonical(...)` |
| `Capabilities` | `@dataclass(frozen=True)` | `name`, `sdk`, `n_qubits`, `native_2q="cz"`, `coupling_map`, `per_shot_usd`, `task_fee_usd`, `median_2q_fidelity`, `execution_windows` |
| `Compiled` | `@dataclass` | `sdk_circuit` (Any), `register_layout: dict[str, tuple[int,...]]`, `name`, `n_2q`, `depth` |
| `Backend` | `Protocol` (runtime-checkable) | attr `capabilities`; `compile(spec, opt_level=3)`, `submit(compiled, shots)`, `resume(handles, compiled)` |
| `Job` | `Protocol` (runtime-checkable) | `result() -> list[Counts]`; property `handles -> list[str]` |

The contract is exercised by [`nwt_substrate/tests/test_qpu.py`](../../nwt_substrate/tests/test_qpu.py) — IR round-trip emission, the canonical-counts equivalence test (the structural fix for the Braket-vs-Qiskit endianness footgun), and an end-to-end decode of the electron walk on the local simulator.

## Quick start

```python
from nwt_substrate.qpu import spec, decode
from nwt_substrate.qpu.adapters import SimulatorBackend

# Electron (2,1) walk: its X/Z Pauli word (from the compendium lookup)
ELECTRON_X = [0, 1, 0, 0, 0, 0, 0]
ELECTRON_Z = [0, 0, 1, 1, 1, 0, 1]

# Build the two destructive-readout circuits (Z-stabs+logical-Z, X-stabs)
base = spec.steane_base_ops(ELECTRON_X, ELECTRON_Z)
z, x = spec.destructive_css(base, "electron")
z.n_qubits, z.count_2q()           # → 7, 9   (9 CZ/CX pre-routing; ~15 routed on Garnet)

# Run on the free local simulator (Braket LocalSimulator under the hood)
backend = SimulatorBackend()
backend.capabilities.sdk           # → 'simulator'
zc = backend.submit([backend.compile(z)], shots=400).result()[0]   # canonical Counts
xc = backend.submit([backend.compile(x)], shots=400).result()[0]

# Decode (vendor-agnostic) → verdict
x_dist, z_dist = decode.destructive_dists(zc.single(), xc.single())
v = decode.verdict_destructive("electron", x_dist, z_dist, "F2", "E2", -1)
v.passed                            # → True
v.x_half.modal                      # → 'F2'        (prob 1.00, noiseless)
v.z_half.modal                      # → ('E2', -1)  (prob 1.00, noiseless)
```

A `BraketBackend("iqm/Garnet")` or `IBMBackend("ibm_fez")` is a drop-in replacement for the `SimulatorBackend` here — same `compile`/`submit`/`result()` surface, same `decode` call.

## API by topic

### Core contracts (`qpu.adapters.base`)

| Symbol | Returns / fields |
|---|---|
| `Counts(register_order, data)` | Frozen canonical counts; `total()`, `register(name) → {bitstring: n}`, `single()`, `from_canonical(order, data)` |
| `Capabilities(...)` | Backend description: `name`, `sdk`, `n_qubits`, `native_2q`, `coupling_map`, `per_shot_usd`, `task_fee_usd`, `median_2q_fidelity`, `execution_windows` |
| `Compiled(sdk_circuit, register_layout, ...)` | SDK circuit + `register_layout` (name → measured qubit indices) + `name`, `n_2q`, `depth` |
| `Backend` (Protocol) | `capabilities`; `compile(spec, opt_level=3) → Compiled`; `submit(compiled, shots) → Job`; `resume(handles, compiled) → Job` |
| `Job` (Protocol) | `result() → list[Counts]`; `handles → list[str]` (per-circuit resume handles) |

### Adapters (`qpu.adapters`)

| Symbol | Notes |
|---|---|
| `BraketBackend(backend="iqm/Garnet", region=None)` | AWS Braket QPU; accepts device name or full ARN. `submit` returns one ARN per circuit. |
| `IBMBackend(backend_name, service=None, token_path=None)` | Qiskit Runtime `SamplerV2`. `handles` is the single batch job-id. |
| `SimulatorBackend(n_qubits=32)` | Braket `LocalSimulator`; immediate results, `handles == []`, `resume` raises. |

### Canonicalization (`qpu.adapters`)

| Function | Signature → behavior |
|---|---|
| `canonicalize_braket(raw_counts, register_layout)` | `{leftmost=qubit0 bitstring: n}` + `{reg: (qubit indices)}` → `Counts`. Pure slice; Braket is already canonical order. SDK-free. |
| `canonicalize_ibm(per_register_bitstrings, register_order)` | `{reg: [little-endian str per shot]}` → `Counts`. **Reverses** each register (`[::-1]`) and zips per shot. SDK-free. |

### Subpackages

| Module | Public surface |
|---|---|
| `spec` | IR `G`, `M`, `CircuitSpec` (`measured_layout()`, `count_2q()`, `to_qiskit()`, `to_braket()`); `steane_base_ops`, `control_base_ops`, `destructive_css`, `ancilla_syndrome`; `STEANE_PREP_H/PREP_CX/STAB/LZ` |
| `decode` | `fano`, `logical_z`, `destructive_dists`, `ancilla_dist`, `verdict_destructive`, `verdict_ancilla`; `HalfResult`, `Verdict`; `QUBIT_TO_FANO` |
| `capabilities` | `preflight`, `Preflight`, `window_status`, `fetch_capabilities`, `REGISTRY` (M2 routing/cost/window guardrail) |
| `runner` | `run`, `build_specs`, `choose_scheme`, `Item` (M3/M4 experiment lifecycle) |
| `experiments` | `exp12_phase2_items`, `PHASE2_TARGETS` (M3 acceptance consumer) |

## Worked examples

### Cross-vendor counts agree (the canonical contract)

A Braket raw key and the matching little-endian Qiskit per-shot string normalize to the *same* `Counts` and decode to the same Fano point:

```python
from nwt_substrate.qpu import canonicalize_braket, canonicalize_ibm, decode

qubit_bits = [1, 0, 1, 1, 0, 0, 1]                          # qubit q -> bit value
braket_key = "".join(str(b) for b in qubit_bits)            # "1011001" (leftmost=qubit0)
qiskit_le  = "".join(str(b) for b in reversed(qubit_bits))  # "1001101" (rightmost=cbit0)

cb = canonicalize_braket({braket_key: 5}, {"c": tuple(range(7))})
ci = canonicalize_ibm({"c": [qiskit_le] * 5}, ("c",))

cb.single()                 # → {'1011001': 5}
ci.single()                 # → {'1011001': 5}   ← same, after IBM reversal
cb.single() == ci.single()  # → True
decode.fano(decode._bits(cb.single().popitem()[0]))   # → 'E3'
decode.fano(decode._bits(ci.single().popitem()[0]))   # → 'E3'  ← identical decode
```

### Multi-register IBM normalization (ancilla scheme)

Two registers, each independently reversed (verified against `test_canonical_counts_multi_register_ibm`):

```python
from nwt_substrate.qpu import canonicalize_ibm

syn_le = "011010"   # cbit5..cbit0  (little-endian)
lz_le  = "1000000"  # cbit6..cbit0
c = canonicalize_ibm({"c_syn": [syn_le], "c_lz": [lz_le]}, ("c_syn", "c_lz"))

c.register_order        # → ('c_syn', 'c_lz')
c.register("c_syn")     # → {'010110': 1}    ← reversed to index-0-left
c.register("c_lz")      # → {'0000001': 1}
c.total()               # → 1
```

### Local-simulator round-trip + decode primitives

```python
from nwt_substrate.qpu import spec, decode
from nwt_substrate.qpu.adapters import SimulatorBackend

backend = SimulatorBackend()
base = spec.steane_base_ops([0,1,0,0,0,0,0], [0,0,1,1,1,0,1])  # electron walk
z, x = spec.destructive_css(base, "electron")
zc = backend.submit([backend.compile(z)], shots=400).result()[0]
xc = backend.submit([backend.compile(x)], shots=400).result()[0]
xd, zd = decode.destructive_dists(zc.single(), xc.single())
decode.verdict_destructive("electron", xd, zd, "F2", "E2", -1).passed   # → True

# Decode primitives in isolation: a valid |0_L> codeword reads trivially
decode.fano([0] * 7)        # → 'I'   (no syndrome)
decode.logical_z([0] * 7)   # → +1
```

## Substrate connection

This layer is what carries the K_7 / SU(2)_5 substrate experiments onto hardware. The substrate's central object is the K_7 graph state `|K_7⟩`, and the per-particle "walk" encodings are run as Steane `[[7,1,3]]` circuits — `7` data qubits, one per K_7 vertex. The numbers that recur here are substrate primitives, not arbitrary:

| Constant | Appears as | Substrate identity |
|---|---|---|
| `7` | data-qubit count; `range(7)` in every Steane spec | K_7 vertex count (`isa.N_VERTICES_K7`) |
| `[[7,1,3]]` | Steane code (`STEANE_STAB`, 3 stabilizer supports) | one logical qubit on the 7-vertex substrate |
| Fano points | `QUBIT_TO_FANO = {0:E3, 1:E2, 2:F1, 3:E1, 4:F2, 5:F3, 6:P}` | Fano-plane labels of the K_7 / PSL(2,7) structure |
| 7-vertex syndrome → particle | `decode.fano` Fano point + `logical_z` | the substrate walk's `(x_fano, z_fano, logical_z)` prediction |

The Hopf-link colored-Jones and full SU(2)_5 modular-S-matrix chip runs (see project memory: WRT D21/D44) are the topology-side experiments this run/decode layer was built to make vendor-portable. `qpu` does not add new physics — it removes the vendor-specific footguns that previously made each chip run a one-off script.

## Papers / references

This is an infrastructure shim, not a physics-prediction paper. The authoritative spec is:

- **Design doc** — `analysis/qpu_interface_design.md` (null-worldtube-private): the full cross-architecture proposal, the canonical `Counts` rationale, the destructive-vs-ancilla scheme decision (2026-05-25 IQM Garnet diagnosis), and the M1–M4 milestone plan.
- **Physics carried** — the substrate experiments these circuits encode are documented in the NWT paper series ([Zenodo NWT community](https://zenodo.org/communities/nwt)) and the [`topology`](topology.md) shim (K_7 Heffter embedding, colored Jones, Hopf links).

## See also

- [`heron`](heron.md) — builds the K_7 graph-state + Steane circuits for IBM Heron (the construction layer above `qpu`)
- [`topology`](topology.md) — K_7 Heffter embedding, colored Jones, Hopf links, SU(2)_5 (the physics these circuits encode)
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`N_VERTICES_K7 = 7`, Fano / PSL(2,7) structure)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — substrate-vs-reference benchmark harness
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
