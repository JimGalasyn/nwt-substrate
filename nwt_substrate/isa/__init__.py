"""
nwt_substrate.isa — substrate-algebraic Instruction Set Architecture.

The 21 generators of so(7), corresponding 1-to-1 with the 21 edges of
K_7, form a universal "active-encoding alphabet" for substrate-derived
observables. This module exposes:

    1. The substrate alphabet: so(7) basis + edge-graph embedding
    2. The kernel: batched trace contractions Tr(M^k) on (N, 7, 7)
       tensors via numpy / torch (CPU/CUDA)
    3. Polynomial assembly: observables as polynomials of trace invariants

Three-layer view (DNA/RNA analog):

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

Quick start:

    >>> import nwt_substrate.isa as isa
    >>> import numpy as np
    >>>
    >>> # Build a batch of 3 ring-adjacency matrices (7×7, 0/1)
    >>> A = np.zeros((3, 7, 7), dtype=int)
    >>> # Coronene-like: W_6 wheel (1 hub + 6-cycle)
    >>> for i in range(1, 7):
    ...     A[0, 0, i] = A[0, i, 0] = 1
    >>> for i in range(1, 6):
    ...     A[0, i, i+1] = A[0, i+1, i] = 1
    >>> A[0, 1, 6] = A[0, 6, 1] = 1
    >>>
    >>> # Trace invariants in a single batched einsum
    >>> inv = isa.graphs_to_invariants(A)
    >>> inv["Tr_M2"][0]
    -24.0
    >>>
    >>> # Aromatic RE prediction with chemistry-supplied π pairs
    >>> pi = np.array([12, 5, 3])  # coronene, pyrene, benzene
    >>> isa.aromaticity_score(A, pi)
    array([200.,  60.,  36.])

Backends:

    >>> isa.available_backends()
    ['numpy', 'torch_cpu', 'torch_cuda']
    >>> isa.aromaticity_score(A, pi, backend="torch_cuda")
    array([200.,  60.,  36.])
"""
from .backends import (
    Backend,
    BackendName,
    available_backends,
    get_backend,
)
from .constants import (
    # Spin(7) rep dimensions
    DIM_V_SPIN7,
    DIM_S_SPIN7,
    DIM_ADJ_SPIN7,
    RANK_SO7,
    H_V_SO7,
    H_COXETER_SO7,
    N_POS_ROOTS_SO7,
    # Substrate α and Macken aspect ratio
    ALPHA_SUBSTRATE,
    ALPHA_NWT,
    M_PL_GEV,
    KAPPA_MACKEN,
    # K_7 combinatorics
    N_VERTICES_K7,
    N_EDGES_K7,
    DEGREE_K7,
    # Derived ratios
    SPINOR_VECTOR_RATIO,
    WILSON_EXPONENT_K7,
    NLO_VERTEX_COEFFICIENT,
    NNLO_BRACKET_COEFFICIENT,
    # Octonion / Cl(0,7) semantics
    DIM_OCTONION,
    N_CL07_IMAGINARY,
    N_LORENTZ_FROM_CL07,
    DIM_INTERNAL_SU2,
    # QCD / SU(3)
    N_C_SU3,
    N_GLUONS,
    C_F_SU3,
    C_A_SU3,
    T_R_SU3,
    # Particle carrier-knot identities
    N_CARRIER_TYPES,
    MAX_CROSSING_NUMBER,
    CARRIER_NAMES,
    # Electroweak / SM fermion content
    N_GENERATIONS,
    N_FERMION_TYPES_PER_GENERATION,
    N_SM_FERMION_FLAVORS,
    N_LORENTZ_GENERATORS,
    B_QED_SM,
    # Casimir signatures
    TR_M2_K7,
    TR_M2_W6,
)
from .so7 import (
    SO7_BASIS,
    so7_basis,
    edge_index,
    k7_edge_list,
    adjacency_to_so7,
    edge_list_to_so7,
    trace_invariants,
    canonical_invariants,
    CANONICAL_GRAPHS,
)
from .batched import (
    adjacency_to_so7_batched,
    trace_invariants_batched,
    graphs_to_invariants,
)
from .observables import (
    aromaticity_score,
    hopf_pair_count,
    k7_indicator,
    k7_wilson_amplitude,
    k7_wilson_breakdown,
    classify_signature,
    CALIBRATION_KCAL_PER_PI_PAIR,
    K7_STABILIZATION_KCAL,
    K7_TR_M2_THRESHOLD,
)

__all__ = [
    # Backends
    "Backend",
    "BackendName",
    "available_backends",
    "get_backend",
    # Substrate structural constants
    "DIM_V_SPIN7",
    "DIM_S_SPIN7",
    "DIM_ADJ_SPIN7",
    "RANK_SO7",
    "H_V_SO7",
    "H_COXETER_SO7",
    "N_POS_ROOTS_SO7",
    "ALPHA_SUBSTRATE",
    "ALPHA_NWT",
    "M_PL_GEV",
    "KAPPA_MACKEN",
    "N_VERTICES_K7",
    "N_EDGES_K7",
    "DEGREE_K7",
    "SPINOR_VECTOR_RATIO",
    "WILSON_EXPONENT_K7",
    "NLO_VERTEX_COEFFICIENT",
    "NNLO_BRACKET_COEFFICIENT",
    "DIM_OCTONION",
    "N_CL07_IMAGINARY",
    "N_LORENTZ_FROM_CL07",
    "DIM_INTERNAL_SU2",
    "N_C_SU3",
    "N_GLUONS",
    "C_F_SU3",
    "C_A_SU3",
    "T_R_SU3",
    "N_CARRIER_TYPES",
    "MAX_CROSSING_NUMBER",
    "CARRIER_NAMES",
    "N_GENERATIONS",
    "N_FERMION_TYPES_PER_GENERATION",
    "N_SM_FERMION_FLAVORS",
    "N_LORENTZ_GENERATORS",
    "B_QED_SM",
    "TR_M2_K7",
    "TR_M2_W6",
    # so(7) primitives
    "SO7_BASIS",
    "so7_basis",
    "edge_index",
    "k7_edge_list",
    "adjacency_to_so7",
    "edge_list_to_so7",
    "trace_invariants",
    "canonical_invariants",
    "CANONICAL_GRAPHS",
    # Batched kernels
    "adjacency_to_so7_batched",
    "trace_invariants_batched",
    "graphs_to_invariants",
    # Observables
    "aromaticity_score",
    "hopf_pair_count",
    "k7_indicator",
    "k7_wilson_amplitude",
    "k7_wilson_breakdown",
    "classify_signature",
    "CALIBRATION_KCAL_PER_PI_PAIR",
    "K7_STABILIZATION_KCAL",
    "K7_TR_M2_THRESHOLD",
]
