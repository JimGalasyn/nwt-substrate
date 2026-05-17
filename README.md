# nwt-substrate

A substrate-algebraic computation library for Null Worldtube Theory (NWT).

`nwt-substrate` is the reference implementation of the substrate algebra
described in the NWT paper series: a Cl(0,7) octonion Clifford algebra
with K_7 graph state on the Heegaard torus of the Brieskorn-Poincaré
sphere S^3 / 2I, supporting particle / scattering / decay /
gravitational-coupling / chemistry computations from a single
internally consistent codebase.

This library is an algebraic continuation of the **photon-vortex programme**
(Williamson & van der Mark 1997 and successors): particles as confined
toroidal structures of the electromagnetic / substrate field, with mass and
quantum numbers emerging from topology. NWT supplies the explicit substrate
algebra — K_7 / so(7) / Spin(7) / Cl(0,7) — that turns the topological
intuition into closed-form quantitative predictions. The
[paper series on Zenodo](https://zenodo.org/communities/nwt) is the
derivation record; this library is the executable companion.

**As of v0.2 (Phase Q.16, 2026-05-11)**, the library ships a
**substrate Instruction Set Architecture** (`nwt_substrate.isa`) that
makes the K_7 algebra load-bearing across every shim. ~25
structural constants (`N_EDGES_K7 = 21`, `N_VERTICES_K7 = 7`,
`DIM_OCTONION = 8`, `RANK_SO7 = 3`, `N_GENERATIONS = 3`,
`N_CARRIER_TYPES = 7`, `B_QED_SM = 8`, …) live in one place,
are asserted at import time, and are consumed by seven
view-shims (chemistry, gravity, qed, qcd, particles, electroweak,
heron). The substrate algebra compiles all the way through to the
21 CZ gates that fire on IBM Heron when `k7_graph_state()` runs.

## Headline predictions

All derived from the substrate algebra at zero free parameters
(beyond the four substrate constants m_e, M_Pl, c, ℏ):

- **Electron mass ratio**: m_e / m_Pl ≈ 4.185 × 10⁻²³ via α^(21/2) Wilson
  amplitude on the K_7 graph state — Paper 17, **−5.5 ppm CODATA**.
  Call: `isa.k7_wilson_amplitude(1/137.036, order="NNLO")`.
- **Newton's G**: 6.674228 × 10⁻¹¹ m³ kg⁻¹ s⁻² via Sakharov-induced
  gravity — Paper 17, **−11 ppm CODATA, inside the ±22 ppm experimental
  band**. Call: `gravity.G_substrate_SI()`.
- **Particle mass spectrum**: 24-particle compendium (hadronic + leptonic
  + exotic) at **0.76 % median residual** — Paper 6 topological mass
  formula. Call: `nwt.particle("p").mass_pred → 937.24 MeV`.
- **Molecular bound states via connected-sum**: deuteron mass-prediction
  residual −0.06 % vs PDG, Pc(4312) +0.013 %, all five tested
  near-threshold molecules within ~0.6 %.
  Call: `nwt.compose(p, n, op="#")`.
- **Coronene aromaticity**: K_7-toroidal resonance energy = 200.0 kcal/mol
  exact (+56 kcal/mol stabilization detected via `Tr(M²) ≤ −24` on the
  K_7 W_6-wheel signature). Call: `chem.smiles_resonance_energy(...)`.
- **Heron quantum-hardware structural verification**: 7 H gates + 21 CZ
  gates fired on IBM Heron, runtime-verified against `isa.N_VERTICES_K7`
  and `isa.N_EDGES_K7`. Call: `heron.k7_graph_state()`.
- **Neutrino sector** (Paper 20, K_8 extension): three active masses
  ≈ (14.8, 17.2, 53) meV, three sterile masses ≈ (61.3, 70.8, 218.8)
  MeV, mixing |U_α4|² ≈ 2.4×10⁻¹⁰, PMNS angles + δ_CP = −2π/3 from
  π_1(PSU(3)) winding. Call: `nwt.neutrino.substrate_breakdown()`.
- **606 substrate tests pass in ~7 seconds**, including 92
  substrate-identity enforcement tests across seven K_7 shims plus
  31 K_8 neutrino-sector tests — the substrate algebra is enforced
  by the codebase, not merely described.

## Install

```bash
pip install nwt-substrate          # not yet on PyPI; for now:
pip install git+https://github.com/JimGalasyn/nwt-substrate.git
```

## Quick start

**Try this first** — three substrate predictions in three lines:

```python
import nwt_substrate as nwt
nwt.particle("p").mass_pred                          # → 937.24 MeV (proton, Paper 6)
nwt.gravity.G_substrate_SI()                         # → 6.674228e-11, -11 ppm CODATA (Paper 17)
nwt.isa.k7_wilson_amplitude(1/137.036, order="NNLO") # → 4.185e-23 = m_e / m_Pl, -5.5 ppm
```

Three independent substrate predictions — particle mass, gravitational
coupling, and the underlying K_7 Wilson amplitude — all matching CODATA
to ppm precision in three function calls.

**The K_7 cross-shim demo** — one substrate, seven shims:

```bash
python3 analysis/isa_cross_shim_demo.py
```

This walks the K_7 algebra through chemistry (coronene aromatic
resonance energy), gravity (m_e/M_Pl via α^(21/2) Wilson amplitude),
qed (8×8 Dirac γ matrices via Cl(0,7) → Cl(1,3)), qcd (8 gluons + 3
colors via Spin(7) ⊃ G_2 ⊃ SU(3)), particles (Paper 6 mass formula on
7 carrier-knot types), electroweak (`b_QED^SM = 8 = DIM_OCTONION`
empirically verified from the SM fermion table), and heron (a qiskit
circuit with exactly 7 H + 21 CZ gates, runtime-verified). Ends with
the substrate identity table showing 8 surfaces in four independent
physics computations.

Particle masses from substrate quantum numbers:

```python
>>> import nwt_substrate as nwt
>>> p = nwt.particle("p")
>>> p.mass_pred
937.24...                             # MeV, Paper 6 mass formula
>>> p.J, p.Q, p.B
(0.5, 1, 1)
>>> p.carrier                          # sourced from isa.CARRIER_NAMES
'cinquefoil'
```

Connected-sum composition law for molecular bound states:

```python
>>> p, n = nwt.particle("p"), nwt.particle("n")
>>> d = nwt.compose(p, n, op="#", name="d", m_obs=1875.61)
>>> d.mass_pred                       # ~1874.48 MeV
>>> d.mass_residual                   # ~ -0.06 % vs PDG
```

Gravitational coupling from substrate alone (now via the ISA):

```python
>>> from nwt_substrate.gravity import G_substrate_SI
>>> G_substrate_SI()                  # 6.674228e-11 m^3 kg^-1 s^-2
                                       # -11 ppm of CODATA, inside ±22 ppm
                                       # experimental error bar
>>> import nwt_substrate.isa as isa
>>> isa.k7_wilson_amplitude(1/137.036, order="NNLO")
4.185439e-23                           # m_e/M_Pl, -5.5 ppm from CODATA
```

Chemistry — aromatic resonance energy from SMILES via batched so(7)
trace invariants:

```python
>>> import nwt_substrate.chemistry as chem
>>> chem.smiles_resonance_energy("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61")
200.0                                  # coronene: K_7-toroidal +56 kcal/mol
                                       # stabilization detected via Tr(M²) ≤ -24
                                       # (= TR_M2_W6 from isa.constants)
```

ISA kernel — substrate-native batched contraction:

```python
>>> import numpy as np
>>> import nwt_substrate.isa as isa
>>> # Build a K_7 adjacency
>>> A = np.ones((1, 7, 7)) - np.eye(7)[None]
>>> inv = isa.graphs_to_invariants(A)
>>> inv["Tr_M2"][0]                    # -42 = -2 × |E(K_7)| (=2×N_EDGES_K7)
-42.0
>>> isa.available_backends()
['numpy', 'torch_cpu', 'torch_cuda']
```

Heron K_7 quantum circuit, structurally verified:

```python
>>> import nwt_substrate.heron as heron
>>> qc = heron.k7_graph_state()
>>> heron.verify_k7_circuit_substrate(qc)
{'n_qubits': 7,    # == N_VERTICES_K7 ✓
 'n_h': 7,         # == N_VERTICES_K7 ✓
 'n_cz': 21,       # == N_EDGES_K7 ✓
 'n_edges_match': True,
 'n_vertices_match': True}
```

## What's implemented

- **`nwt_substrate.isa`** — Substrate Instruction Set Architecture
  (v0.2 new). Central source of truth for K_7 / Spin(7) / so(7) /
  Cl(0,7) structural constants, with import-time assertions enforcing
  identities like `N_EDGES_K7 = DIM_ADJ_SPIN7 = 21`,
  `4 + 3 = N_VERTICES_K7 = 7`, `B_QED_SM = DIM_OCTONION = 8`. Backends:
  numpy, torch CPU, torch CUDA. Observables: `aromaticity_score`,
  `hopf_pair_count`, `k7_indicator`, `k7_wilson_amplitude`,
  `classify_signature`. **Batched einsum kernel runs at 2 ns/molecule
  on CUDA, 1124× faster than networkx graph traversal.**
- **Particles** — Paper 6 mass formula, charge via extended GMN, the
  full SM hadronic + leptonic + exotic catalog. `Particle` class
  validates `n_q ∈ [0, MAX_CROSSING_NUMBER]` against ISA at construction
  time.
- **Compositions** — knot connected-sum (#) for molecular bound states
  (deuteron, X(3872), Pc family), Hopf-link with Λ_QCD = 313 MeV per
  crossing for nuclear / strongly-bound exotic regimes.
- **Walk-phase scattering** — substrate-algebraic Compton (matches
  Klein-Nishina to 1e-9), Møller / Bhabha, V-A muon decay matching
  Sargent rate, neutron decay with g_A = 1.27.
- **Chemistry** (v0.2 new) — SMILES → substrate Hopf-pair aromaticity
  RE with K_7-toroidal correction; Clar sextets via maximum-independent-
  set; McKay-admissible coordinations; C_60 vibrational mode decomposition
  (174 = 4 IR + 10 Raman + ...); batched ISA-backed RE for ≥10^5 SMILES.
- **Gauge-theory shims** — `nwt.qed`, `nwt.qcd` (incl. gg→gg),
  `nwt.electroweak` (Z resonance + chiral couplings + `b_QED^SM`
  verification), `nwt.qft` (Lagrangian view), `nwt.string` (string-
  theoretic view), `nwt.gravity` (Sakharov-induced G via
  `isa.k7_wilson_amplitude`). Every shim has a `substrate_breakdown()`
  function printing its substrate-identity table.
- **Heron experiments** — qiskit-runtime interface and an experiment
  registry for IBM Heron processors. Supports Experiments 4 / 5 / 9
  / 10 / 11 from the paper series. K_7-circuit gate counts are
  runtime-verified against `isa.N_VERTICES_K7` / `isa.N_EDGES_K7`.
- **Neutrino sector** (v0.3 new, K_8 extension for Paper 20) — closed-form
  active masses (Wilson amplitude on K_8 with `N_v=8, N_e=28`), sterile
  masses (Wilson amplitude with `N_v=8, N_e=19` from the Z_3 ⊂ G_2
  triality seesaw, edge difference 9 = 12 − 3), `|U_α4|² = α^(9/2)`,
  PMNS angles at leading order from Spin(8) triality, and `δ_CP = −2π/3`
  from π_1(PSU(3)) winding. K_8 structural constants
  (`N_VERTICES_K8 = 8`, `N_EDGES_K8 = 28`, `K8_PARTITION = (6,3,12,1,6)`,
  `K8_SEESAW_EDGE_DIFFERENCE = 9`) live in `isa.constants` alongside
  the K_7 family.
- **Diagrams** — programmatic figure factories for the canonical
  substrate visualisations (torus knots, K_7 traversals,
  Heegaard-torus unification).

## The active-encoding architecture

The library has three layers:

```
┌─────────────────────────────────────────────────────────────┐
│  SUBSTRATE (passive primitives)                              │
│  isa.constants — 25 K_7/Spin(7)/so(7)/Cl(0,7) structural    │
│                  integers, import-time-asserted             │
│  isa.so7       — 21-generator basis + edge-graph embedding  │
├─────────────────────────────────────────────────────────────┤
│  ISA (active encoding — the substrate ribosome)              │
│  isa.batched     — einsum kernels: numpy / torch CPU / CUDA  │
│  isa.observables — polynomial-of-trace-invariants assembly   │
├─────────────────────────────────────────────────────────────┤
│  SHIMS (translation to domain vocabularies)                  │
│  chemistry, qed, qcd, electroweak, particles, gravity,      │
│  heron — each turns its domain vocabulary into so(7) input  │
│  and consumes ISA constants for cross-shim consistency      │
└─────────────────────────────────────────────────────────────┘
```

Spectacular cross-shim identities the architecture surfaces:

- **`N_EDGES_K7 = 21`** appears in **seven** shims:
  - chemistry: 21 so(7) generators in ISA basis
  - gravity: α^(21/2) Wilson amplitude
  - qed: 21 = dim(so(7) adjoint) holding the γ-matrix algebra
  - qcd: 21 = dim(adjoint Spin(7)) ⊃ SU(3) gluons
  - particles: 21 - 9 = 12 mixed so(7) gens host SM flavors
  - electroweak: `21 = 6 (Lorentz) + 3 (internal) + 12 (flavors)`
  - heron: 21 CZ gates in `k7_graph_state()` on real hardware
- **`8 = DIM_OCTONION`** appears in **four independent physics
  computations**:
  - gravity: numerator of `SPINOR_VECTOR_RATIO = 8/7`
  - qed: shape of Dirac γ^μ = 8×8
  - qcd: number of gluons = N_c² - 1 = 8
  - electroweak: `b_QED^SM = Σ N_c × Q² = 8` empirically verified
- **`7 = N_VERTICES_K7`** appears in **five** shims
  (chemistry, gravity, qed, particles, heron)

If a refactor violates any of these identities in any one shim, the
92 cross-shim tests catch it across all seven shims simultaneously.
The substrate algebra is no longer described by the code — it is
enforced by it.

## Tests

```bash
pytest nwt_substrate/tests/ -q
# 606 passed in ~7s
```

This includes:
- 92 cross-shim tests (`test_isa_cross_shim.py`) enforcing K_7 algebra
  across chemistry, gravity, qed, qcd, particles, electroweak, heron
- 47 ISA-internal tests across 3 backends (numpy / torch_cpu /
  torch_cuda)
- 58 chemistry tests (SMILES parsing, K_7 hub detection on coronene,
  Clar sextets, McKay coordinations, C_60 vibrational)
- All pre-Phase-Q.16 tests preserved (zero numerical regressions in
  particle masses, scattering cross-sections, gravity prediction,
  EW Z-resonance, etc.)

## Citation

If you use this library in a publication, please cite both:

- The relevant NWT paper(s) — typically one of
  [Paper 14–19](https://zenodo.org/communities/nwt) for the result
  you're using.
- The library Zenodo record (auto-archived per release):

```bibtex
@software{nwt_substrate,
  author       = {Galasyn, Jim and others},
  title        = {{nwt-substrate}: a substrate-algebraic computation
                  library for Null Worldtube Theory},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.20041585}
}
```

A `CITATION.cff` is included in this repo for tools that auto-resolve
software citations.

## Papers

The library implements the computations described in:

- **Paper 6** — topological mass formula (0.76 % median residual on the
  24-particle compendium after the 2026-04-30 nucleon update).
- **Paper 14** — α^(21/2) heptafoil amplitude.
- **Paper 15** — Wilson amplitude on K_7 graph state.
- **Paper 16** — NWT three-field Lagrangian (BPS critical coupling).
- **Paper 17** — m_e / m_Pl closed form: G to -11 ppm CODATA (inside
  the ±22 ppm experimental band).
- **Paper 18** — Sakharov-induced Einstein gravity from substrate
  matter sector. *Includes the canonical "Heegaard torus, two
  sectors" figure rendered by `nwt.diagrams.figure_paper18_unified()`.*
- **Paper 19** — substrate monism via library demonstration.
- **Paper 20** — neutrino sector from Spin(8) triality on K_7 / K_8.
  Three sterile masses {61.3, 70.8, 218.8} MeV in the νMSM window,
  |U_α4|² ≈ α^(9/2) ≈ 2.4×10⁻¹⁰ active-sterile mixing, PMNS angles
  from triality, δ_CP = −2π/3 from π_1(PSU(3)) Z_3 winding. Library
  implementation in `nwt_substrate.neutrino`; K_8 structural
  constants in `isa.constants`. DOI:
  [10.5281/zenodo.20259632](https://doi.org/10.5281/zenodo.20259632).

The Zenodo community for the full series is at
https://zenodo.org/communities/nwt (collected DOIs).

## Status

**v0.2 (Phase Q.16, 2026-05-11)**: the active-encoding architecture.

API surface is stable for particles, compositions, walk_phase, gauge
shims, gravity, chemistry, diagrams, and the new ISA layer. Minor
breaking changes may still occur before 1.0; we aim for semver
discipline post-1.0.

The main private development monorepo, where new analyses and paper
drafts live before promotion, is `null-worldtube-private` (not public).
Polished analyses and paper-supporting computations are promoted to
this repo; exploratory work stays private.

## Contributing

Issues and pull requests welcome. Please run the test suite
(`pytest nwt_substrate/tests/`) before submitting, and include a
short note describing the physics motivation for any new feature
(this is a physics library; please don't add tooling that has no
substrate-algebraic content).

If your contribution touches the K_7 substrate algebra, please update
`nwt_substrate/isa/constants.py` rather than introducing magic numbers
in shim code — the cross-shim tests in `tests/test_isa_cross_shim.py`
will catch identity violations across all seven shims.

## License

MIT. See [LICENSE](./LICENSE).
