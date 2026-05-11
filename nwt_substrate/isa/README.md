# `nwt_substrate.isa` — substrate-algebraic Instruction Set Architecture

The substrate ISA layer is the universal compute primitive for NWT.
The 21 generators of so(7) correspond one-to-one with the 21 edges of
K_7, and every substrate-derived observable factors through trace
contractions on antisymmetric 7×7 matrices.

This module exposes:

1. The substrate **alphabet** — so(7) basis + graph embedding
2. The substrate **kernel** — batched trace contractions Tr(M^k)
3. **Polynomial assembly** — observables as polynomials of trace invariants
4. **Backend dispatch** — numpy, torch_cpu, torch_cuda

## DNA/RNA analog

The substrate-ISA architecture is the codebase version of the
Level 0 → Level 1 abiogenesis transition. The substrate (so(7))
is passive; the ISA is the active-encoding layer that turns
substrate elements into observable predictions.

```
┌─────────────────────────────────────────────────────────────┐
│  SUBSTRATE (passive)                                         │
│  21 so(7) generators ↔ K_7 edges ↔ "amino acid alphabet"     │
├─────────────────────────────────────────────────────────────┤
│  ISA (active encoding)        ← THIS MODULE                  │
│  Batched einsum on shape-(N, 7, 7), trace polynomials        │
├─────────────────────────────────────────────────────────────┤
│  SHIMS (translation)                                         │
│  chemistry/qed/qcd/electroweak/gravity/qft/string/heron      │
│  Each shim turns its domain vocabulary into so(7) input.     │
└─────────────────────────────────────────────────────────────┘
```

## Quick start

```python
import numpy as np
import nwt_substrate.isa as isa

# Build a batch of 3 ring-adjacency matrices (7×7, 0/1)
A = np.zeros((3, 7, 7), dtype=int)

# Coronene-like: W_6 wheel (1 hub + 6-cycle around it)
for i in range(1, 7):
    A[0, 0, i] = A[0, i, 0] = 1
for i in range(1, 6):
    A[0, i, i+1] = A[0, i+1, i] = 1
A[0, 1, 6] = A[0, 6, 1] = 1

# Compute trace invariants in a single batched einsum
inv = isa.graphs_to_invariants(A)
print(inv["Tr_M2"][0])   # -24.0 (= -2 × 12 edges of W_6)

# Aromatic RE prediction with chemistry-supplied π pairs
pi = np.array([12, 5, 3])  # coronene, pyrene-shape, benzene-shape
print(isa.aromaticity_score(A, pi))   # [200., 60., 36.]
```

## API summary

### so(7) primitives (`isa.so7`)

```python
isa.SO7_BASIS                        # (21, 7, 7) generators
isa.so7_basis()                      # rebuild basis
isa.edge_index(a, b)                 # K_7 edge → so(7) basis index
isa.adjacency_to_so7(adj)            # symmetric A → antisymmetric M
isa.edge_list_to_so7(edges)          # build M from edge list
isa.trace_invariants(M)              # Tr(M²), Tr(M⁴), Tr(M⁶)
isa.canonical_invariants()           # reference signatures
```

### Batched kernel (`isa.batched`)

```python
isa.adjacency_to_so7_batched(adj_batch, backend="numpy")
isa.trace_invariants_batched(M_batch, orders=(2, 4, 6), backend="numpy")
isa.graphs_to_invariants(adj_batch, backend="numpy")    # one-shot
```

### Observables (`isa.observables`)

```python
isa.aromaticity_score(adj, pi_pairs, backend="numpy")
isa.hopf_pair_count(adj, pi_electrons, backend="numpy")
isa.k7_indicator(adj, n_vertices_used=None, backend="numpy")
isa.classify_signature(adj, backend="numpy")
```

### Backends (`isa.backends`)

```python
isa.available_backends()              # ['numpy', 'torch_cpu', 'torch_cuda']
isa.get_backend("torch_cuda")         # Backend descriptor
```

## Adding a new observable

Any observable that factors through trace invariants of an
antisymmetric 7×7 matrix can be added in `observables.py` in three
steps:

1. Define the polynomial:

   ```python
   def my_observable(adj_batch, ..., backend="numpy"):
       inv = isa.graphs_to_invariants(adj_batch, orders=(2, 4, 6),
                                       backend=backend)
       # polynomial in inv["Tr_M2"], inv["Tr_M4"], ...
       return my_polynomial(inv)
   ```

2. Add a verification test in `tests/test_isa.py` that compares
   ISA output to a reference implementation (e.g. the legacy
   chemistry-shim version, or a hand-computed value on the
   reference set).

3. Register in `isa/__init__.py` if it's part of the public API.

## What the substrate ISA can / can't do

✓ **What it does well:**
- Compute trace contractions Tr(M^k) on batches of (N, 7, 7) tensors
  in nanoseconds-per-molecule on CPU or GPU
- Verify cross-shim consistency: the same trace invariants are
  computed identically across chemistry, qed, qcd, etc.
- Amortize parse + extraction cost over multiple observables
  computed on the same parsed batch

✗ **What it doesn't do:**
- Accelerate domain-specific parsing (SMILES, scattering kinematics)
- Help with single-molecule calls — the kernel is faster than the
  function-call overhead at N=1
- Replace richer linear-algebra operations (full diagonalization,
  determinants) where they're genuinely needed by an observable

## Performance characteristics

Measured on the chemistry-aromaticity benchmark (analysis/isa_scale_benchmark.py):

| Scale  | numpy kernel | torch_cpu kernel | torch_cuda kernel |
|--------|--------------|------------------|-------------------|
| 1k     | 2 µs/mol     | 3 µs/mol         | 2 µs/mol          |
| 10k    | 2 µs/mol     | 433 ns/mol       | 539 ns/mol        |
| 100k   | 2 µs/mol     | 506 ns/mol       | 315 ns/mol        |
| 1M     | 2 µs/mol     | 395 ns/mol       | 346 ns/mol        |

Kernel-only speedup over the legacy networkx path is 500–5000× depending
on scale and backend. End-to-end speedup is bottlenecked by SMILES
parsing.

For 7×7 matrices, GPU acceleration is modest — the kernel is too small
to saturate the device. GPU wins materially when many observables run
on the same parsed batch, or when so(n) is extended to n ≥ 16 for
larger graph projections.

## Related

- `analysis/so7_aromaticity_probe.py` — original probe that motivated this module
- `analysis/isa_scale_benchmark.py` — scale benchmark code
- `nwt_substrate.chemistry.smiles.smiles_resonance_energy_batch` — first shim integration
- `memory/so7-substrate-isa-probe.md` — design memo
