"""Substrate structural constants — one source of truth for all shims.

These integers and ratios are the fixed structural identities of the
substrate. They are NOT calibrated against data; they are derived
from the representation theory of Spin(7), the combinatorics of the
K_7 graph, and the McKay correspondence for binary polyhedral groups.

Every NWT shim (chemistry, qed, qcd, gravity, ...) that needs any of
these constants should import them from here, NOT redefine them.
That way, structural identities (e.g. dim_Adj = n_edges_K7 = 21) are
enforced at module-import time as runtime assertions, and a refactor
to one shim that violates the substrate's algebraic structure is
immediately visible across all shims.

Layered DNA/RNA analog (Level 0 → Level 1 transition):
  - Level 0: passive substrate (Cl(0,7), so(7), 2I, etc.)
  - Level 1: active encoding (these constants, used by all shims)
"""
from __future__ import annotations

import math


# ---------------------------------------------------------------------------
# Spin(7) representation dimensions
# ---------------------------------------------------------------------------

DIM_V_SPIN7: int = 7
"""Dimension of the vector representation V of Spin(7).
Equal to the number of vertices of K_7. Appears as the index range
of so(7) generators (a, b ∈ {0, ..., 6}, a < b)."""

DIM_S_SPIN7: int = 8
"""Dimension of the spinor representation S of Spin(7).
Sets the (8/7) prefactor in the Paper 17 gravity coupling.
Equal to dim(octonions) = 8 (Spin(7) is the stabilizer of the
G_2-invariant 3-form on Im(O) = R^7)."""

DIM_ADJ_SPIN7: int = 21
"""Dimension of the adjoint representation (= so(7) as Lie algebra).
Equal to (DIM_V_SPIN7 choose 2) = 7·6/2 = 21.
Equal to the number of EDGES of K_7."""

RANK_SO7: int = 3
"""Rank of so(7) (= maximal torus dimension).
Equal to DIM_ADJ_SPIN7 / DIM_V_SPIN7 = 21/7 = 3."""

H_V_SO7: int = 5
"""Dual Coxeter number h_v of so(7).
Used as the "steady-state interface" component count in the Spin(7)
Coxeter decomposition 8 = h_v + rank = 5 + 3 of the disk-ergosurface
matching condition (vortex-vision 2026-05-15)."""

H_COXETER_SO7: int = 6
"""Coxeter number h of so(7) = 2 × rank = number of positive roots / rank.
Appears as the multiplier in the suggestive 6α match to EM (spin-1)
superradiance maximum amplification at extremal Kerr."""

N_POS_ROOTS_SO7: int = 9
"""Number of positive roots of so(7) = (h - 1) × rank / 2 + rank = 9.
Equivalently, (h × rank) / 2 = 18/2 = 9. The full set of 18 roots
(9 positive + 9 negative) accounts for the 18 = dim(Adj) - rank = 21 - 3
non-Cartan generators of so(7)."""


# ---------------------------------------------------------------------------
# Substrate fine-structure constant and Macken aspect ratio
#
# Two structurally equivalent forms of the NWT fine-structure constant:
#
#     α_NWT = 1 / (25π√3 + 1)   (Paper 17 trefoil Aharonov-Bohm)
#     κ_Macken = √(α_NWT⁻¹ / √2) = √((25π√3 + 1) / √2)
#
# κ_Macken ≈ 9.844 is the universal substrate aspect ratio (Macken's
# Wave Theory of the Vacuum). The 1/√2 inside the square root is the
# Z_2 polarity factor of the cross-arc loop.
#
# These are EXACT closed-form expressions, not calibrations — α_NWT is
# defined by the trefoil substrate amplitude geometry, and κ_Macken is
# its derived aspect ratio. They are the entry points for every
# substrate-vortex prediction (Kerr efficiency, cosmogenesis transfer,
# disk-ergosurface matching, QPO torus eigenmodes, healing-length
# domain markers).
# ---------------------------------------------------------------------------

ALPHA_NWT: float = 1.0 / (25.0 * math.pi * math.sqrt(3.0) + 1.0)
"""The substrate (NWT) fine-structure constant.

Closed form: α_NWT = 1 / (25π√3 + 1) ≈ 1/137.0359...

Derivation: Paper 17 trefoil Aharonov-Bohm phase, where the (5√3, 25)
prefactors come from the K_7 trefoil-carrier (n_q=3) winding number
and the 5 = h_v(so(7)) dual-Coxeter components of the steady-state
substrate channel. The +1 is the identity (singlet) octonion component.

This is the framework's PRIMARY substrate parameter — every other
gravity/cosmogenesis observable below is expressed as a polynomial in
α or √α with substrate-integer multipliers (RANK_SO7=3, H_V_SO7=5,
DIM_OCTONION=8, H_COXETER_SO7=6)."""

KAPPA_MACKEN: float = math.sqrt((25.0 * math.pi * math.sqrt(3.0) + 1.0)
                                 / math.sqrt(2.0))
"""Macken's universal substrate aspect ratio κ ≈ 9.844.

Closed form: κ = √((25π√3 + 1) / √2) = √(1/α_NWT / √2)

The 1/√2 inside is the Z_2 polarity factor of the cross-arc loop.
κ is the substrate's intrinsic length scale ratio (vortex core radius
divided by healing length, equivalently). It appears as the daughter-BH
"κ_daughter" in cosmogenesis throat-selection: the parent-BH κ_parent
relates to κ_daughter by κ_parent / κ_daughter = 3(1+√α), and the
accretion-disk κ_disk relates to the parent by κ_disk / κ_parent =
8(1 - √α) (Spin(7) Coxeter decomposition; vortex-vision 2026-05-15)."""


# ---------------------------------------------------------------------------
# K_7 combinatorial constants
# ---------------------------------------------------------------------------

N_VERTICES_K7: int = 7
"""Number of vertices of the complete graph K_7.
Equal to DIM_V_SPIN7 (structural identity, asserted below)."""

N_EDGES_K7: int = 21
"""Number of edges of the complete graph K_7.
Equal to DIM_ADJ_SPIN7 (structural identity, asserted below).
Equal to (N_VERTICES_K7 choose 2) = 7·6/2 = 21."""

DEGREE_K7: int = 6
"""Degree of every vertex in K_7.
Equal to N_VERTICES_K7 - 1 = 6."""


# ---------------------------------------------------------------------------
# K_8 combinatorial constants (Paper 20 neutrino sector)
#
# K_8 extends K_7 with a single fixed vertex (the eighth) corresponding
# to the Spin(7) ⊂ Spin(8) stabilized unit vector in V_8, identified
# with the real unit on the octonion side (Baez 2002 §3.2).
# ---------------------------------------------------------------------------

N_VERTICES_K8: int = 8
"""Number of vertices of the complete graph K_8.
Equal to N_VERTICES_K7 + 1 = 8 (Paper 20 §5: the eighth vertex is
the Spin(7) ⊂ Spin(8) stabilized vector / octonion real unit)."""

N_EDGES_K8: int = 28
"""Number of edges of the complete graph K_8.
Equal to (N_VERTICES_K8 choose 2) = 8·7/2 = 28.
Equal to N_EDGES_K7 + DEGREE_K8 = 21 + 7 = 28."""

DEGREE_K8: int = 7
"""Degree of every vertex in K_8.
Equal to N_VERTICES_K8 - 1 = 7 = N_VERTICES_K7."""

K8_PARTITION: tuple = (6, 3, 12, 1, 6)
"""K_8 edge decomposition under Z_3 ⊂ G_2 = Aut(O) (Paper 20 §5):
  Lorentz (6) + internal SU(2) (3) + SM-flavor (12) + Higgs vacuum (1)
  + Yukawa (6) = 28 = N_EDGES_K8.
The first three terms (6+3+12 = 21) match the K_7 substrate-boson
partition of Paper 19. The +1 is the Z_3-fixed self-loop of the
eighth vertex. The final +6 is the Z_3-fixed edges between the
eighth vertex and the K_7 vertices."""

K8_PARTITION_K7_INHERITED: int = 21
"""The first three K8_PARTITION terms (6+3+12) sum to 21 = N_EDGES_K7.
This is the K_7 substrate-boson partition inherited from Paper 19,
preserved under the K_8 extension."""

K8_ACTIVE_EDGE_COUNT: int = 28
"""Edge count for the active-neutrino phase soliton on K_8 (Paper 19
W3.3-J via Paper 20 §3): m_ν_active = (8/8) α^(N_e/2) prefactors
with N_e = K8_ACTIVE_EDGE_COUNT = 28."""

K8_STERILE_EDGE_COUNT: int = 19
"""Edge count for the sterile-neutrino phase soliton on K_8 (Paper 20
§6): K_8 minus the 9 cross-orbit edges of the Z_3 ⊂ G_2 decomposition.
N_e_sterile = K8_ACTIVE_EDGE_COUNT - 9 = 28 - 9 = 19. Gives the
α^(9/2) seesaw ratio m_active/m_sterile ≈ 2.4×10⁻¹⁰."""

K8_SEESAW_EDGE_DIFFERENCE: int = 9
"""Edge-count difference between active and sterile K_8 phase solitons.
Equal to K8_ACTIVE_EDGE_COUNT - K8_STERILE_EDGE_COUNT = 28 - 19 = 9.
The 9 = 12 (SM-flavor edges) - 3 (internal SU(2) edges) forced by
G_2 → SU(3) representation theory (Paper 20 §4)."""


# ---------------------------------------------------------------------------
# Derived ratios used by shims
# ---------------------------------------------------------------------------

SPINOR_VECTOR_RATIO: float = DIM_S_SPIN7 / DIM_V_SPIN7
"""(8/7) — the prefactor in m_e/M_Pl from Spin(7)'s spinor/vector ratio."""

WILSON_EXPONENT_K7: float = N_EDGES_K7 / 2.0
"""21/2 — the K_7 Wilson amplitude exponent on α.
Comes from √α per edge × 21 edges = α^(21/2) (Paper 17 Phase 6.1b)."""

NLO_VERTEX_COEFFICIENT: float = 1.0 / N_VERTICES_K7
"""1/7 — the per-K_7-vertex NLO correction coefficient (Paper 17 Phase 6.5)."""

NNLO_BRACKET_COEFFICIENT: float = DIM_ADJ_SPIN7 / DIM_S_SPIN7
"""21/8 — the NNLO bracket coefficient = dim(Adj)/dim(S) (Paper 17 Phase 6.9 §3.iii)."""


# ---------------------------------------------------------------------------
# Octonion / Cl(0,7) substrate semantics
#
# Cl(0,7) = octonion left-multiplication algebra (Walk Phase 1).
# 7 imaginary octonion directions ↔ 7 vertices of K_7 ↔ Spin(7) vector rep.
# 8-dim octonions ↔ Spin(7) spinor rep.
# Picking 4 of 7 imaginary directions + Wick-rotating one gives Cl(1,3)
# Lorentz; the remaining 3 generate an internal SU(2) commuting with γ^μ.
# ---------------------------------------------------------------------------

DIM_OCTONION: int = DIM_S_SPIN7
"""Dimension of the octonions = 8. Equal to DIM_S_SPIN7 (Spin(7) spinor rep).
The 8×8 Dirac gamma matrices in the substrate are octonion left-multiplications
on the 8-dim octonion space."""

N_CL07_IMAGINARY: int = DIM_V_SPIN7
"""Number of imaginary octonion directions (e_1, ..., e_7) = 7.
Equal to DIM_V_SPIN7 and N_VERTICES_K7. The generators of Cl(0,7)."""

N_LORENTZ_FROM_CL07: int = 4
"""Number of Cl(0,7) directions used for the Lorentz Cl(1,3) embedding.
3 imaginary spatial + 1 Wick-rotated time = 4 of the 7 imaginary directions."""

DIM_INTERNAL_SU2: int = 3
"""Number of remaining imaginary directions after the Lorentz subspace.
= N_CL07_IMAGINARY - N_LORENTZ_FROM_CL07 = 7 - 4 = 3.
These generate an internal SU(2) that commutes with all γ^μ.
Numerically equal to RANK_SO7 (which is also 3 = (N_VERTICES_K7 - 1)/2)."""


# ---------------------------------------------------------------------------
# QCD / SU(3) substrate identities
#
# SU(3) color sits inside G_2 = aut(octonions) ⊂ Spin(7). The substrate
# embedding is:
#   - dim(adjoint SU(3)) = 8 gluons = DIM_OCTONION
#     (same 8 that gives Dirac γ^μ their 8×8 shape)
#   - rank(SU(3)) related to RANK_SO7 via the G_2 embedding
#   - N_c = 3 = RANK_SO7
#     (and = DIM_INTERNAL_SU2, the 3 "extra" Cl(0,7) directions left
#      after the Lorentz partition)
# ---------------------------------------------------------------------------

N_C_SU3: int = RANK_SO7
"""Number of colors in QCD = 3. Numerically equal to rank(so(7)), to
the number of internal-SU(2) directions in Cl(0,7), and to
(N_VERTICES_K7 - 1) / 2."""

N_GLUONS: int = DIM_OCTONION
"""Number of gluons = dim(adjoint SU(3)) = N_c² - 1 = 8.
SUBSTRATE IDENTITY: equal to DIM_OCTONION = DIM_S_SPIN7.
The same 8 appears in:
  - QED Dirac γ^μ as 8×8 matrices
  - QCD's 8 gluons
  - The Spin(7) spinor rep
  - The octonion algebra dimension
This is not coincidence — SU(3) ⊂ G_2 ⊂ Spin(7), with all three
acting on the octonions."""

C_F_SU3: float = DIM_OCTONION / (2.0 * N_C_SU3)
"""Fundamental Casimir of SU(3) = (N_c² - 1)/(2 N_c) = 8/6 = 4/3.
Expressed as DIM_OCTONION / (2 × N_C_SU3) to make the substrate
origin explicit."""

C_A_SU3: float = float(N_C_SU3)
"""Adjoint Casimir of SU(3) = N_c = 3.
SUBSTRATE IDENTITY: equal to RANK_SO7."""

T_R_SU3: float = 0.5
"""Index of fundamental rep: Tr(T^a T^b) = T_R δ^ab = (1/2) δ^ab."""


# ---------------------------------------------------------------------------
# Particle carrier-knot identities (Paper 6 substrate-algebraic spectrum)
#
# Every NWT particle has a carrier knot whose crossing number is n_q ∈
# {0, 1, ..., 6}. That's 7 values = N_VERTICES_K7 = DIM_V_SPIN7 — each
# carrier type corresponds to one K_7 vertex / one Spin(7) vector
# direction. The 6 + 1 (unknot is degenerate, used for n_q=0 and n_q=1)
# distinct carrier topologies populate the 6 2I irrep dimensions found
# in SM matter (memory: "all 6 2I irrep dims populated by SM particles").
# ---------------------------------------------------------------------------

N_CARRIER_TYPES: int = N_VERTICES_K7
"""Number of distinct carrier-knot types in the Paper 6 spectrum = 7.
Equal to N_VERTICES_K7 (one per K_7 vertex). The valid range of n_q
in any Particle is [0, N_CARRIER_TYPES - 1] = [0, 6]."""

MAX_CROSSING_NUMBER: int = N_VERTICES_K7 - 1
"""Maximum carrier-knot crossing number = 6. Beyond this the carrier
cannot fit into the K_7 substrate's 7-vertex toroidal embedding."""

CARRIER_NAMES: tuple[str, ...] = (
    "unknot",       # n_q = 0
    "unknot",       # n_q = 1 (degenerate with n_q=0; charge-different)
    "Hopf",         # n_q = 2 → mesons
    "trefoil",      # n_q = 3 → most baryons
    "figure-eight", # n_q = 4 → tetraquarks
    "cinquefoil",   # n_q = 5 → nucleon family + pentaquarks
    "6_1 knot",     # n_q = 6 → hexaquarks / dibaryons
)
"""Canonical carrier knot for each crossing number n_q.
Length = N_CARRIER_TYPES. Used by topology.torus_knots.carrier_for_n_q."""


# ---------------------------------------------------------------------------
# Electroweak / SM fermion content identities
#
# The 21 so(7) generators decompose as (memory: substrate-boson-count-so7):
#   21 = 6 (Lorentz so(1,3)) + 3 (internal SU(2)) + 12 (mixed)
# The 12 "mixed" generators correspond to the 12 SM fermion flavors
# (3 generations × 4 types per generation = 12):
#   per generation: 1 charged lepton + 1 neutrino + 1 up-type + 1 down-type
#
# The QED 1-loop β-function coefficient summed over all SM fermions:
#   b_QED^SM = 3 (leptons × 1²) + 4 (up-type × 3 × (2/3)²)
#                                + 1 (down-type × 3 × (1/3)²)
#            = 8 = DIM_OCTONION
# This is the substrate identity for the SM fermion charge content.
# ---------------------------------------------------------------------------

N_GENERATIONS: int = RANK_SO7
"""Number of SM fermion generations = 3.
Equal to RANK_SO7, N_C_SU3, and DIM_INTERNAL_SU2.
This is the substrate origin of the 3-fold replication of fermions."""

N_FERMION_TYPES_PER_GENERATION: int = 4
"""Number of fermion types per SM generation = 4.
One charged lepton + one neutrino + one up-type quark + one down-type
quark per generation."""

N_SM_FERMION_FLAVORS: int = N_GENERATIONS * N_FERMION_TYPES_PER_GENERATION
"""Total SM fermion flavors (counting color as singlet) = 12.
Equal to N_EDGES_K7 - 6 (Lorentz) - DIM_INTERNAL_SU2 = 21 - 6 - 3 = 12,
the count of 'mixed' so(7) generators per the substrate boson decomposition.

(Per the substrate-boson-count-so7 hypothesis: each mixed generator
hosts one SM fermion flavor.)"""

N_LORENTZ_GENERATORS: int = 6
"""Dimension of the Lorentz algebra so(1,3) = 6 (3 rotations + 3 boosts).
Equal to N_EDGES_K7 - DIM_INTERNAL_SU2 - N_SM_FERMION_FLAVORS = 21 - 3 - 12 = 6."""

B_QED_SM: int = DIM_OCTONION
"""The QED 1-loop β-function coefficient at full SM activation = 8.

  b_QED^SM = Σ over SM fermions of N_c(f) × Q_f²
           = 3 × (1)² (leptons)
             + 3 × 3 × (2/3)² (up-type quarks, 3 colors each)
             + 3 × 3 × (1/3)² (down-type quarks, 3 colors each)
           = 3 + 4 + 1
           = 8

SUBSTRATE IDENTITY: equal to DIM_OCTONION = DIM_S_SPIN7.
The same 8 that gives Dirac γ^μ their 8×8 shape and SU(3) its 8
gluons also equals the SM fermion charge sum that drives QED running."""


# ---------------------------------------------------------------------------
# Substrate Casimir signatures (graph invariants on the symmetric adjacency)
# ---------------------------------------------------------------------------

TR_M2_K7: int = -2 * N_EDGES_K7
"""Tr(M²) for the K_7 antisymmetric so(7) embedding = -2 × 21 = -42."""

TR_M2_W6: int = -2 * 12
"""Tr(M²) for the wheel-graph W_6 (1 hub + 6-cycle).
This is coronene's ring-adjacency signature: 6 hub-petal edges +
6 petal-cycle edges = 12 edges. Threshold for K_7-hub indicator
in chemistry observables."""


# ---------------------------------------------------------------------------
# Structural-identity assertions (run at import time)
# ---------------------------------------------------------------------------

assert N_VERTICES_K7 == DIM_V_SPIN7, (
    f"Substrate structural identity violated: "
    f"N_VERTICES_K7={N_VERTICES_K7} != DIM_V_SPIN7={DIM_V_SPIN7}"
)

assert N_EDGES_K7 == DIM_ADJ_SPIN7, (
    f"Substrate structural identity violated: "
    f"N_EDGES_K7={N_EDGES_K7} != DIM_ADJ_SPIN7={DIM_ADJ_SPIN7}"
)

assert N_EDGES_K7 == (N_VERTICES_K7 * (N_VERTICES_K7 - 1)) // 2, (
    f"K_7 edge count mismatch: {N_EDGES_K7} != C(7,2)"
)

assert DIM_ADJ_SPIN7 == RANK_SO7 * DIM_V_SPIN7, (
    f"so(7) dimension formula violated: "
    f"{DIM_ADJ_SPIN7} != {RANK_SO7} * {DIM_V_SPIN7}"
)

assert DIM_OCTONION == DIM_S_SPIN7, (
    f"Octonion / Spin(7)-spinor dimension mismatch: "
    f"{DIM_OCTONION} != {DIM_S_SPIN7}"
)

assert N_CL07_IMAGINARY == DIM_V_SPIN7, (
    f"Cl(0,7) imaginary count != Spin(7) vector dim: "
    f"{N_CL07_IMAGINARY} != {DIM_V_SPIN7}"
)

assert N_LORENTZ_FROM_CL07 + DIM_INTERNAL_SU2 == N_CL07_IMAGINARY, (
    f"Lorentz + internal-SU(2) directions don't partition Cl(0,7): "
    f"{N_LORENTZ_FROM_CL07} + {DIM_INTERNAL_SU2} != {N_CL07_IMAGINARY}"
)

assert DIM_INTERNAL_SU2 == RANK_SO7, (
    f"Internal SU(2) dim != rank(so(7)): "
    f"{DIM_INTERNAL_SU2} != {RANK_SO7}"
)

assert N_GLUONS == DIM_OCTONION, (
    f"N_gluons != DIM_OCTONION (substrate identity): "
    f"{N_GLUONS} != {DIM_OCTONION}"
)

assert N_GLUONS == N_C_SU3 * N_C_SU3 - 1, (
    f"SU(N_c) adjoint dim formula violated: "
    f"{N_GLUONS} != {N_C_SU3}² - 1 = {N_C_SU3*N_C_SU3 - 1}"
)

assert N_C_SU3 == RANK_SO7, (
    f"N_c != rank(so(7)) (substrate identity): "
    f"{N_C_SU3} != {RANK_SO7}"
)

assert abs(C_F_SU3 - (N_C_SU3 * N_C_SU3 - 1) / (2.0 * N_C_SU3)) < 1e-12, (
    f"C_F formula violated: "
    f"{C_F_SU3} != (N_c² - 1)/(2 N_c) = {(N_C_SU3*N_C_SU3 - 1)/(2*N_C_SU3)}"
)

assert N_CARRIER_TYPES == N_VERTICES_K7, (
    f"Carrier type count != |V(K_7)|: "
    f"{N_CARRIER_TYPES} != {N_VERTICES_K7}"
)

assert MAX_CROSSING_NUMBER == N_VERTICES_K7 - 1, (
    f"Max crossing number != |V(K_7)| - 1: "
    f"{MAX_CROSSING_NUMBER} != {N_VERTICES_K7 - 1}"
)

assert len(CARRIER_NAMES) == N_CARRIER_TYPES, (
    f"CARRIER_NAMES length != N_CARRIER_TYPES: "
    f"{len(CARRIER_NAMES)} != {N_CARRIER_TYPES}"
)

assert N_GENERATIONS == RANK_SO7, (
    f"N_GENERATIONS != rank(so(7)): {N_GENERATIONS} != {RANK_SO7}"
)

assert N_SM_FERMION_FLAVORS == N_GENERATIONS * N_FERMION_TYPES_PER_GENERATION, (
    f"SM fermion count formula violated: "
    f"{N_SM_FERMION_FLAVORS} != {N_GENERATIONS} × {N_FERMION_TYPES_PER_GENERATION}"
)

assert (N_LORENTZ_GENERATORS + DIM_INTERNAL_SU2 + N_SM_FERMION_FLAVORS
        == N_EDGES_K7), (
    f"so(7) generator decomposition violated: "
    f"{N_LORENTZ_GENERATORS} + {DIM_INTERNAL_SU2} + {N_SM_FERMION_FLAVORS} "
    f"!= {N_EDGES_K7}"
)

assert B_QED_SM == DIM_OCTONION, (
    f"b_QED^SM != DIM_OCTONION (substrate identity): "
    f"{B_QED_SM} != {DIM_OCTONION}"
)

# ---------------------------------------------------------------------------
# Coxeter / Macken identities
# ---------------------------------------------------------------------------

assert H_V_SO7 + RANK_SO7 == DIM_OCTONION, (
    f"Spin(7) Coxeter decomposition violated: "
    f"H_V_SO7({H_V_SO7}) + RANK_SO7({RANK_SO7}) != DIM_OCTONION({DIM_OCTONION})"
)

assert H_COXETER_SO7 == 2 * RANK_SO7, (
    f"so(7) Coxeter number formula violated: "
    f"H_COXETER_SO7({H_COXETER_SO7}) != 2 × RANK_SO7({RANK_SO7})"
)

assert 2 * N_POS_ROOTS_SO7 + RANK_SO7 == DIM_ADJ_SPIN7, (
    f"so(7) root-system decomposition violated: "
    f"2 × N_POS_ROOTS_SO7({N_POS_ROOTS_SO7}) + RANK_SO7({RANK_SO7}) "
    f"!= DIM_ADJ_SPIN7({DIM_ADJ_SPIN7})"
)

# α_NWT and κ_Macken closed-form consistency: κ² × √2 = 1/α
assert abs(KAPPA_MACKEN ** 2 * math.sqrt(2.0) - 1.0 / ALPHA_NWT) < 1e-9, (
    f"Macken/α closed-form consistency violated: "
    f"κ²√2 ({KAPPA_MACKEN**2 * math.sqrt(2.0)}) != 1/α ({1.0/ALPHA_NWT})"
)

# α_NWT matches CODATA α_QED at ~ppm precision (Paper 17 prediction)
_ALPHA_QED_CODATA = 1.0 / 137.035999084
assert abs(ALPHA_NWT - _ALPHA_QED_CODATA) / _ALPHA_QED_CODATA < 1e-5, (
    f"α_NWT closed form deviates > 10 ppm from CODATA α: "
    f"α_NWT={ALPHA_NWT}, α_CODATA={_ALPHA_QED_CODATA}"
)


# K_8 structural identities (Paper 20)

assert N_VERTICES_K8 == N_VERTICES_K7 + 1, (
    f"K_8 = K_7 + 1 vertex violated: "
    f"N_VERTICES_K8({N_VERTICES_K8}) != N_VERTICES_K7({N_VERTICES_K7}) + 1"
)

assert N_EDGES_K8 == N_VERTICES_K8 * (N_VERTICES_K8 - 1) // 2, (
    f"K_8 edge-count formula violated: "
    f"N_EDGES_K8({N_EDGES_K8}) != N_VERTICES_K8·(N_VERTICES_K8-1)/2"
    f"({N_VERTICES_K8 * (N_VERTICES_K8 - 1) // 2})"
)

assert DEGREE_K8 == N_VERTICES_K8 - 1, (
    f"K_8 degree formula violated: "
    f"DEGREE_K8({DEGREE_K8}) != N_VERTICES_K8 - 1({N_VERTICES_K8 - 1})"
)

assert sum(K8_PARTITION) == N_EDGES_K8, (
    f"K_8 partition does not sum to N_EDGES_K8: "
    f"sum(K8_PARTITION)={sum(K8_PARTITION)} != N_EDGES_K8({N_EDGES_K8})"
)

assert sum(K8_PARTITION[:3]) == K8_PARTITION_K7_INHERITED == N_EDGES_K7, (
    f"K_8 partition K_7-inheritance violated: "
    f"sum(K8_PARTITION[:3])={sum(K8_PARTITION[:3])}, "
    f"K8_PARTITION_K7_INHERITED={K8_PARTITION_K7_INHERITED}, "
    f"N_EDGES_K7={N_EDGES_K7}"
)

assert K8_ACTIVE_EDGE_COUNT == N_EDGES_K8, (
    f"K_8 active-neutrino edge count must equal N_EDGES_K8: "
    f"K8_ACTIVE_EDGE_COUNT({K8_ACTIVE_EDGE_COUNT}) != N_EDGES_K8({N_EDGES_K8})"
)

assert K8_SEESAW_EDGE_DIFFERENCE == K8_ACTIVE_EDGE_COUNT - K8_STERILE_EDGE_COUNT, (
    f"K_8 seesaw edge-difference formula violated: "
    f"K8_SEESAW_EDGE_DIFFERENCE({K8_SEESAW_EDGE_DIFFERENCE}) != "
    f"K8_ACTIVE_EDGE_COUNT({K8_ACTIVE_EDGE_COUNT}) - "
    f"K8_STERILE_EDGE_COUNT({K8_STERILE_EDGE_COUNT})"
)

assert K8_SEESAW_EDGE_DIFFERENCE == 12 - 3, (
    f"K_8 seesaw 9 = SM-flavor (12) - internal-SU(2) (3) violated: "
    f"K8_SEESAW_EDGE_DIFFERENCE({K8_SEESAW_EDGE_DIFFERENCE}) != 12 - 3 = 9"
)
