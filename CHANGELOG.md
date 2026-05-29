# Changelog

All notable changes to `nwt-substrate` are recorded here.

The format follows [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/), and the project follows [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html). Narrative release notes for tagged versions live under [`docs/releases/`](docs/releases/).

## [Unreleased]

### Benchmarks — correction layers (follow-up to P. Kaboth's v0.3.0 report)

Pasquale Kaboth flagged the two largest benchmark deviations (Z-boson width, muon lifetime) as places where "additional coupling or correction layers may be required." Both are now resolved — and the diagnosis differs between them:

- **Z-boson width** was pure LO (2.93% low). Hadronic (N_c = 3) partial widths now carry the standard massless-quark QCD radiative correction `R_QCD = 1 + a + 1.409 a² − 12.77 a³` (`a = α_s/π`), via new `electroweak.decays.qcd_correction_factor()`. Γ_Z: **2.4223 → 2.4874 GeV (2.93% → 0.31% of PDG)**; BR(Z→ℓℓ) shifts from 3.443% toward PDG 3.366%. Leptonic and invisible widths are unchanged.
- **Muon lifetime** (was 9.94%) is the **m⁵ amplification of the Paper-6 m_μ residual** (m_μ is −1.97%; τ ∝ m_μ⁵ → ≈ 5× → ~9.8%), *not* a missing weak layer. `benchmark_muon_decay_rate` — the weak-sector closure, isolated from the mass formula with PDG m_μ — gains the standard phase-space `f(x)` + 1-loop QED radiative correction layer: **0.45% (tree) → 0.007%**. `benchmark_muon_lifetime` is re-labeled the *compound* probe and its deviation correctly attributed to the mass formula, cross-referencing the weak-closure benchmark. Both corrections are parameter-free (α is substrate-predicted).
- Pasquale's second point — sensitivity of predictions to small changes in the structural integers — is answered by the v0.3.1 ISA sensitivity sweep (`python -m nwt_substrate.sensitivity`).

### Sensitivity — structural-criticality layer (M. Wende's d12rg insight)

Marcel Wende observed that *sensitivity* (how strongly an observable responds) is not the only interesting quantity — *structural criticality* (which integers carry the most load, and which observables move together) reveals the architecture. Added to `SensitivityReport`:

- `structural_load()` — integers ranked by the number of observables each moves (the load-bearing vs localized distinction).
- `comovement(min_shared)` — benchmark pairs that move *together* under the same integers, ranked by shared count: correlated structural dependencies, even where the individual shifts are modest.
- `criticality_summary()` and a `python -m nwt_substrate.sensitivity --criticality` flag that prints both layers on top of the raw counts.

### Self-consistency audit (library-internal, value-preserving)

A full audit of the library against the canonical v0.3.1 `isa/constants.py`:

- **Fixed**: `dark_sector.wimp_98gev` docstring σ_SI tier table (now matches the module's own formula and tests — only the α⁴ tier survives LZ-2024); `string.is_trefoil` no longer returns True for the unknot T(2,1)/(1,2); `gravity.healing_length` ξ_cosmo docstring (≈100 kpc, was an inconsistent 30); `particles.stability_ratio` docstring (22 catalogued particles; removed a ρ-meson example absent from the catalog); pruned 26 stale ghost exports from `electroweak.__all__` (meson decay constants live in `particles.decay_constants` since v0.3.1); `amplitudes.vertices.G_WEAK` single-sourced from `electroweak.G_W` (was 0.6536 vs canonical 0.6495); a dead `M_Pl` assignment in `gravity.coupling`.
- **Changed**: the neutrino-sector NNLO Wilson bracket unified on the canonical *additive* `(1 + α/7 + (21/8)α²)` form (was multiplicative, α²-coefficient 3); two `gravity.einstein` docstrings likewise.
- **Added**: `isa.N_TRIANGLES_K7 = C(7,3) = 35` (import-time asserted), closing an orphan structural integer the chemistry f-block rule depended on. Routed ~a dozen bare structural literals through `isa` (SU(3) Casimirs, QED-loop trace dims, `clifford` K8 edges, `colored_jones` level, cosmology Λ-exponent, gravity G prefactor, qft colour/gluon counts, CKM 7³/√7). De-duplicated `K7_CORRECTION_KCAL_PER_RING_SET`, `STEANE_LZ`, and `G_F`.
- **Fixed**: a latent `amplitudes ↔ electroweak ↔ qcd` circular import (a module-level cross-package import made function-local, matching the existing `algebra.su3`/`dirac` idiom).

## [0.3.1] - 2026-05-28

[Full narrative release notes](docs/releases/v0.3.1.md). Version DOI [10.5281/zenodo.20438240](https://doi.org/10.5281/zenodo.20438240), concept DOI [10.5281/zenodo.20012027](https://doi.org/10.5281/zenodo.20012027).

Prompted by the first independent reproduction of the v0.3.0 benchmark suite (P. Kaboth, d12rg). This release adds a sensitivity-analysis tool and — in building it — found and fixed two predictions that were reporting measured inputs or overstated accuracies. A total-transparency patch.

### Added
- `nwt_substrate.sensitivity`: **ISA structural-integer sensitivity sweep**. `integer_sweep()` perturbs each leaf `NAME: int = V` constant in `isa/constants.py` by ±1 in a fresh `python -O` subprocess (`-O` strips the import-time structural asserts), recomputes all 38 benchmarks, and reports which predictions move — the robustness / look-elsewhere map. `python -m nwt_substrate.sensitivity` prints the table; `SensitivityReport` exposes `.movers()`, `.inert_integers`, `.never_moved`, `.coupling()`. Transiently source-patches `isa/constants.py` and always restores it (writable / editable install required).
- Per-shim reference docs for the final 10 shims — `qed`, `qcd`, `qft`, `neutrino`, `atomic`, `dark_sector`, `chemistry`, `topology`, `heron`, `qpu` — completing the v0.3.0 "remaining documentation" item. Every shim now has an API reference page; the `docs/index.md` status table is all "written".

### Changed
- `electroweak/substrate_gf.py`, `electroweak/substrate_ckm.py`: prediction formulas now **import their structural integers from `isa`** instead of re-encoding them as bare literals (`25 → H_V_SO7**2`, `625 → H_V_SO7**4`; `N_VERTICES_K7` and `DIM_ADJ_SPIN7` imported). Value-preserving. Predictions coupled to the ISA in the sensitivity sweep rose **17/38 → 21/38**; `benchmark_fermi_constant`, `benchmark_higgs_vev`, and `benchmark_full_ckm` are now load-bearing under integer perturbation.

### Fixed
- **`benchmark_sin2_theta_W` now reports the honest substrate prediction** `sin²θ_W = (2 + α)/(DIM_OCTONION + 1) = (2 + α)/9 = 0.22303` at its true **~3.5 %** residual vs PDG 0.23122. Earlier versions reported the hardcoded **measured** effective angle (`SIN2_THETA_W = 0.23121`, ~43 ppm) as if it were the prediction. The measured angle is retained — now explicitly labeled — as the Z-pole coupling **input** (Z observables unchanged). New `SIN2_THETA_W_SUBSTRATE`. The overstated accuracy claim is corrected across `README.md`, `llms.txt`, `llms-full.txt`, `docs/index.md`, `docs/FAQ.md`, `docs/shims/electroweak.md`, and `benchmarks/README.md`.
- `docs/shims/electroweak.md` Z-width row corrected to `Γ_Z = 2.4223 GeV` (2.9 %), not the previously-claimed `2.4979 GeV` (0.1 %) — surfaced by the same reproduction.
- `qft.Lagrangian.beta_0()` now recognises a generic `SU(N)` Yang-Mills gauge tag (previously returned `None` for `yang_mills(N)`); `yang_mills(3).beta_0(n_f_dirac=6) = 7.0`.
- `neutrino` docstring example masses updated to the current NuFIT/PDG Δm²₃₁ anchor (m₃ 53.00 → 52.25 meV, m_N3 218.8 → 215.7 MeV).

## [0.3.0] - 2026-05-28

[Full narrative release notes](docs/releases/v0.3.0.md). Version DOI [10.5281/zenodo.20435950](https://doi.org/10.5281/zenodo.20435950), concept DOI [10.5281/zenodo.20012027](https://doi.org/10.5281/zenodo.20012027).

### Added
- `nwt_substrate.benchmarks` subpackage: **38 forward-prediction benchmarks** spanning particle physics, atomic physics, QED/QCD, electroweak precision, cosmology, gravity, black-hole thermodynamics, and chemistry. Full suite runs in ~100 ms via `python -m nwt_substrate.benchmarks` and serves as the empirical anti-numerology argument made concrete. Each benchmark returns a `BenchmarkResult` dataclass with substrate timing, substrate value, accuracy vs reference, and the traditional method that would otherwise be needed.
- `nwt_substrate.benchmarks.run_all()` and a CLI entry point (`python -m nwt_substrate.benchmarks --summary | --json | --max-time-ms`).
- `nwt_substrate.dark_sector.wimp_98gev`: L_NWT Higgs-portal calculation for the 98 GeV WIMP — `WIMP_98GeV` dataclass, `sigma_si_higgs_portal`, `lhc_production_cross_section`, `predict_all`. Substrate prediction `σ_SI ~ 10⁻⁴⁶ cm²`, testable at LZ-G3.
- `nwt_substrate.particles.decay_constants`: unified canonical home for **light + heavy + vector + B_c** meson decay constants (P7b §2-3, §7.5, §7.6). Consolidates the previously-duplicated three-module `electroweak.{light,heavy,vector}_meson_decay_constants` into one tested module with dataclass specs (`LightPseudoscalarSpec`, `HeavyMesonSpec`, `VectorMesonSpec`), name-based lookups (`light_meson_fX_for`, `heavy_meson_fX_for`, `vector_meson_fX_for`), ratio diagnostics (`fX_ratio_strange_over_nonstrange`, `c_ratio_precision_chain`), per-sector precision chains, and unified verification. 18 mesons total covering π, K, η, D, D_s, B, B_s, ρ, ω, K*, φ, J/ψ, Υ, D*, D_s*, B*, B_s*, B_c.
- `nwt_substrate.heron.sidereal_geometry`: Exp 11 directional layer with Observatory dataclass, lab → ICRS rotation, K_7 lab-frame embedding, predicted-sigma-pattern forward models (symmetric + asymmetric-jet variants), directional match score, LST anchoring, and multi-day folding helpers.
- `nwt_substrate.gravity.nhek`: Near-Horizon Extremal Kerr geometry. Metric, inverse, determinant, signature, Killing structure, Christoffels (numerical + symbolic via sympy), vacuum-Einstein verification, substrate-vortex centerline + bifurcation-sphere helpers. Background for the cosmogenic bridge in Paper 22.
- `nwt_substrate.cosmology` shim: `eta_B` (`(3/14)α⁴`, 0.38 % Planck), `omega_b_c` (`25α + 75α²`, 67 ppm Planck — 240× tighter than the measurement), `lambda_cc` (`(m_e/M_Pl)⁴ · α¹⁶ · h_Cox`, 0.74 % — closes the 123-orders-of-magnitude problem). Plus `anisotropy_axes` submodule with canonical CMB anomaly directions, citations, great-circle fits, and the 12-pair AoE × HPA separation check used by the Paper 22 cosmogenesis program.
- `nwt_substrate.qpu`: vendor-neutral QPU adapter (M1 + M2 + M3 + M4). `Spec`, preflight guardrail, experiment runner, IBM / AWS Braket / simulator adapters, destructive CSS readout as cross-vendor default.
- `nwt_substrate.condensate`: abelian-Higgs + Bogoliubov phases A through K-2 (vortex-line tension, σ-orbit Wilson dynamics, NR superfluid ξ = λ̄_C, walk-to-Pauli lookup, K_7 walks at L ≤ 25).
- `nwt_substrate.algebra.codes.k8`: K_8 Steane-code framework used by `cosmology.eta_B` derivation.
- `nwt_substrate.qcd.exotic_states`: substrate universal mass formula `m² = (4 m_π⁰)² · N` for exotic hadronic states.
- `nwt_substrate.atomic.hydrogen`: substrate hydrogen spectroscopy from Paper 8 Coulomb (Bohr radius, Lyman α, 21 cm, Lamb shift, Rydberg constant).
- `nwt_substrate.algebra.clifford`: substrate Pauli theorem `L²(8) = 7 + 21 = 28 = |E(K_8)|`.
- `nwt_substrate.topology.colored_jones`: colored-Jones / Rosso-Jones cabled-knot states for torus knots.
- `nwt_substrate.particles.stability_ratio`: substrate pattern-stability ratio `ρ = m/Γ` classification (passive / BPS / active), with `all_k7_walks_are_passive_or_BPS` consistency check.
- `nwt_substrate.particles.charge`: Gell-Mann-Nishijima `Q = T₃ + Y/2` integer-pair consistency.
- `nwt_substrate.particles`: topological isospin ceiling `I ≤ n_q / 2` (WRT D17, non-breaking).
- CI: GitHub Actions `tests` workflow on Python 3.10, 3.11, 3.12 with pytest-cov; `benchmarks` workflow runs `run_all()` and uploads JSON artifacts.
- Codecov integration with `codecov.yml` patch target 80 %, hardware-SDK paths in ignore list. Coverage gated explicitly via `CODECOV_TOKEN`.
- AEO/AIO-friendly documentation: `llms.txt`, `llms-full.txt`, `AGENTS.md`, `docs/FAQ.md`, `docs/index.md`, `docs/shims/{gravity,particles,electroweak,cosmology}.md`.
- `CONTRIBUTING.md` (this file's companion) for human contributors.
- `astropy>=6.0` added to `[test]` extras (unlocks heron.sidereal_geometry tests in CI).
- Production-code test coverage pass: `test_gravity_nhek.py` (23 tests), `test_heron_sidereal_geometry.py` (23 tests), `test_qed_qcd_diagram.py` (37 tests), `test_particles_decay_constants.py` (30 tests). Total test count: **1352** (was 1239).

### Changed
- Codecov line coverage: **73.67 % → 80.03 %** since v0.2.0 tag, via dedup of decay-constants modules + new tests + honest ignore of hardware-SDK paths.
- `nwt_substrate.electroweak.__init__`: meson decay constants are no longer re-exported here; canonical home is `nwt_substrate.particles.decay_constants`. The `_SubstrateNamespace` meson-decay convenience bindings removed (would create a circular import with `electroweak.substrate_gf`).
- `nwt_substrate/benchmarks/compute_speed.py`: `benchmark_decay_constants` and `benchmark_vector_meson_decay` now import from `nwt_substrate.particles.decay_constants` (was `electroweak.{light,heavy,vector}_meson_decay_constants`).

### Removed
- `nwt_substrate.electroweak.heavy_meson_decay_constants` — consolidated into `particles.decay_constants`.
- `nwt_substrate.electroweak.light_meson_decay_constants` — consolidated into `particles.decay_constants`.
- `nwt_substrate.electroweak.vector_meson_decay_constants` — consolidated into `particles.decay_constants`.

These were dead code at the import-graph level (only the now-deleted test files imported them); no external consumer is affected.

## [0.2.0] - 2026-05-26

[Full narrative release notes](docs/releases/v0.2.0.md)

### Added
- **`nwt_substrate.isa`**: substrate Instruction Set Architecture — the central source of truth for K_7 / Spin(7) / so(7) / Cl(0,7) structural constants. About 25 integers (`N_VERTICES_K7 = 7`, `N_EDGES_K7 = 21`, `DIM_OCTONION = 8`, `RANK_SO7 = 3`, `N_GENERATIONS = 3`, `N_CARRIER_TYPES = 7`, `B_QED_SM = 8`, …) live in one place, are asserted at import time to satisfy structural identities, and are consumed by every shim. This makes cross-shim consistency mechanical rather than coincidental.
- `nwt_substrate.isa.k7_wilson_amplitude(alpha, order="LO"|"NLO"|"NNLO")`: the substrate K_7 Wilson amplitude — drives `m_e/M_Pl`, Newton's `G`, the cosmological constant, and other derivations across multiple shims.
- `nwt_substrate.electroweak.substrate_gf`: substrate Fermi constant `G_F` closure (P7b §7.4). 55 ppm vs PDG.
- `nwt_substrate.electroweak.substrate_ckm`: substrate Wolfenstein CKM (λ, A, ρ̄, η̄, δ_CP, J, V_us, V_cb, V_ub, V_td, full 3×3 matrix) from Paper 6b's `λ² = 7α` Cabibbo-Wilson amplitude.
- `nwt_substrate.electroweak.form_factors`: substrate Dalitz `f_+(0)` form factors (K → π, D → π, D → K, B → π, B → D) via `cos^N(θ_C)` substrate Cabibbo law (P7b §7.7).
- `verify_b_qed_sm` substrate-identity check: `b_QED^SM = Σ N_c × Q²` over all SM fermions equals `isa.B_QED_SM = DIM_OCTONION = 8` exactly.
- 92 substrate-identity enforcement tests across seven K_7 shims plus 31 K_8 neutrino-sector tests.
- Cross-shim consistency checks: same `N_EDGES_K7 = 21` integer appears in qed (running coefficient), gravity (`α^(21/2)` Wilson amplitude), chemistry (K_7 hub stabilization), and substrate ISA.

### Changed
- Released as v0.2.0 — version DOI [10.5281/zenodo.20398451](https://doi.org/10.5281/zenodo.20398451), concept DOI [10.5281/zenodo.20012027](https://doi.org/10.5281/zenodo.20012027).
- `CITATION.cff` and README updated for v0.2.0 release metadata.
- The K_7 substrate algebra now compiles all the way through to the 21 CZ gates that fire on IBM Heron when `k7_graph_state()` runs — gate counts are runtime-verified against `isa.N_VERTICES_K7` and `isa.N_EDGES_K7`.

## [0.1.3] - 2026-05-05

### Added
- Algebra-to-picture color mapping for Feynman diagrams: the `Diagram.render_color_mapped()` method renders each amplitude term in a chosen color and emits a matching colored diagram so the substrate-algebra expression and the Feynman graph are visually paired.
- `nwt_multiview_demo.py` analysis script demonstrating the "one electron, 7 lenses" framing — the same particle viewed through each shim.

### Changed
- README refreshed for the v0.1.3 release.

## [0.1.2] - 2026-05-04

### Fixed
- **NNLO bracket coefficient** in `m_e/M_Pl`: was `(3) α²`, corrected to `(21/8) α²` (from PSL(2,7) Fano-plane analysis). With the correction, substrate Newton's `G` lands at **−11 ppm CODATA**, inside the ±22 ppm experimental band. Previously: a few ppm off, but with the wrong structural coefficient.

This is a substrate-algebra correction, not a fit. The coefficient `21/8 = dim Adj(so(7)) / dim S(Spin(7))` is forced by the substrate; the prior value was a derivation slip.

## [0.1.1] - 2026-05-03

### Added
- Paper 18 G1-G6 reproduction scripts: full Sakharov-induced gravity derivation chain (linearized and full nonlinear Einstein equations from the NWT substrate condensate).

## [0.1.0] - 2026-05-03

### Added
- **Initial public release** of `nwt-substrate`.
- Core substrate-algebra library: K_7 graph state on the Heegaard torus of the Brieskorn-Poincaré sphere `S³ / 2I`, with Cl(0,7) octonion Clifford algebra and `so(7)` gauge structure.
- `nwt_substrate.particles`: Paper 6 carrier-knot mass formula, 24-particle compendium at 0.76 % median residual.
- `nwt_substrate.gravity`: substrate Newton's G via Sakharov-induced derivation, black-hole thermodynamics, Kerr efficiency, FLRW backreaction null tests (Buchert `Q₀`), cosmogenesis (Thorne `a* = 0.998` equilibrium).
- `nwt_substrate.electroweak`: SM coupling table, Z partial / total widths, `σ(e⁺e⁻ → ff̄)` via γ + Z + interference.
- `nwt_substrate.qed`, `nwt_substrate.qcd`: gauge-theory shims (Schwinger 1-loop `a_e`, R-ratio, exotic states).
- `nwt_substrate.chemistry`: aromaticity / NICS / C_60 combinatorics, SMILES → substrate resonance energy.
- `nwt_substrate.heron`: qiskit-runtime IBM Heron interface; `k7_graph_state()` 7 H + 21 CZ circuit, runtime-verified against substrate constants.
- `nwt_substrate.neutrino`: K_8 extension (3 active + 3 sterile masses, PMNS, `δ_CP = −2π/3` from `π₁(PSU(3))` winding).
- `nwt_substrate.topology`: K_7 Heffter genus-1 toroidal embedding, torus knots, colored Jones polynomials.
- `nwt_substrate.compositions`: connected-sum molecule construction; deuteron, P_c(4312) reproduced within sub-percent.
- Test suite covering the core derivations.
- MIT License.

---

## Comparing versions

```bash
git log v0.1.0..v0.2.0     # commits between any two tags
git log v0.2.0..HEAD       # all unreleased commits
```

## How to cite a specific version

See [`CITATION.cff`](CITATION.cff). For pinned-version reproducibility:

```bibtex
@software{nwt_substrate,
  author  = {Galasyn, Jim},
  title   = {nwt-substrate: a substrate-algebraic computation library for Null Worldtube Theory},
  version = {0.3.0},
  doi     = {10.5281/zenodo.20435950},
  url     = {https://github.com/JimGalasyn/nwt-substrate}
}
```

Concept DOI (always-latest, version-agnostic): [`10.5281/zenodo.20012027`](https://doi.org/10.5281/zenodo.20012027).

[Unreleased]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/JimGalasyn/nwt-substrate/releases/tag/v0.1.0
