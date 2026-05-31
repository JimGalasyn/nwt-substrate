# Changelog

All notable changes to `nwt-substrate` are recorded here.

The format follows [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/), and the project follows [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html). Narrative release notes for tagged versions live under [`docs/releases/`](docs/releases/).

## [Unreleased]

### Added

- **Classical EM fields from substrate sources (`nwt_substrate.em`).** FFT–Poisson solvers for the electric and magnetic fields a particle's carrier radiates from its charge density ρ and supercurrent density j: `electric_field` (∇²φ = −ρ, E = −∇φ) and `magnetic_field` (∇²A = −μ₀j, B = ∇×A, Biot–Savart), with `divergence`/`curl` spectral helpers, source deposition on a carrier curve (`deposit_sources`), and field-line tracing (`trace_field_lines`). **`maxwell_eh`** adds the optional **Euler–Heisenberg** nonlinear-vacuum correction (`L_EH = ξ[(E²−B²)² + 7(E·B)²]` → polarization/magnetization, solved by fixed-point iteration); `eh_xi=0` is exact linear Maxwell, and the physical coupling is ~α³ at the Compton scale (negligible — boost `eh_xi` to probe the strong-field regime). **`multipole_moments`** returns charge / electric-dipole / magnetic-dipole — the bridge to electromagnetic form factors and elastic scattering, so the module is reusable beyond visualisation. Natural units (ε₀=μ₀=1). Both **`bc="periodic"`** (FFT, neutralizing background — exact near the source, Gauss-law-clean) and **`bc="open"`** (zero-padded free-space Green's function — no image charges, clean Coulomb field lines at all latitudes) are supported. 11 tests (Coulomb falloff, open-vs-periodic falloff, Gauss's law, dipole axiality, ∇·B=0, EH-reduces-to-linear, multipole recovery).
- **Particle portraits — render a particle from its field configuration (`nwt_substrate.portraits`).** A new module renders each NWT particle as a glowing, phase-coloured filament built from its actual field content: the exact BPS Nielsen-Olesen vortex profile `f(ρ)` (new `condensate.solve_bps_vortex`, the first-order Bogomolny solution at λ = e²/2) is bent along the carrier-knot curve and the abelian-Higgs phase `p·φ + q·θ` is composited as an emission-absorption volume. `portraits.portrait(p, q)` / `portrait(hopf=True)` returns a Figure; `portraits.gallery()` renders the carrier zoo (unknot/lepton, Hopf/meson, trefoil/baryon, cinquefoil/nucleon, …); `portraits.n_q_to_knot()` maps a carrier crossing number to its knot. Reproduce with `python diagrams/particle_portraits.py`.
- **`condensate.solve_bps_vortex` — exact BPS vortex profile.** The radial cross-section `f(ρ), a(ρ)` of the unit (n=1) abelian-Higgs vortex, by shooting on the first-order Bogomolny equations (converges to leading coefficient c1 ≈ 0.6033; the scalar decays at the BPS rate 1/ξ). Complements the existing BPS *scalars* in `condensate.abelian_higgs` (healing length, line tension) with the actual profile; returned as a `BPSVortexProfile` with interpolating `f_at`/`a_at`. Higher windings are documented as future work (shooting is stiff for n≥2).
- **`diagrams.torus_knot_curve` / `diagrams.hopf_link_curves` — knot-curve primitives.** Return the (p, q) torus-knot and Hopf-link curves as point arrays (previously only drawable, not retrievable); `draw_torus_knot` now builds on `torus_knot_curve`.
- **Closure-priority layer — `o10.closure_priority` + `--priority` (M. Wende's prioritization ask).** The directional dual of the diagnostic layers: instead of "what breaks the system?", it answers "where does the next month of effort buy the most structural closure?". Each benchmark scores a closure weight (0 ungrounded / 0.5 single-sector hidden SPOF / 1.0 ≥2 independent sectors, −0.5 if a marked defect edge); the suite currently sits at **23.5 / 38 (62%)**. Open items split into **[A] effort-known** (ranked by `ROI = closure-gain / effort-tier` — verified measured-input leaks are tier-1 mechanical, plain hidden SPOFs and defect edges tier-2) and **[B] a triage queue** of benchmarks the sweep's leaf-integer perturbations can't reach, whose effort is genuinely unknown until a code read separates a cheap leak from a real derivation gap (so the layer does *not* pretend to know — surfacing that distinction is the point). Validated against the v0.4.2 α-fix, which was exactly a tier-1 leak worth +2.0 closure. `python -m nwt_substrate.benchmarks.o10 --suite --sensitivity --priority`.

### Fixed

- **`benchmark_muon_decay_rate` — removed a circular measured-input dependence (closure-priority follow-up).** The closure-priority layer flagged this as the last hidden SPOF; investigating it surfaced a deeper leak than the headline α: the benchmark computed Γ = G_F²m_μ⁵/192π³ using the **measured PDG G_F** (`electroweak.G_F_GEV = 1.1663787e-5`), which is itself *extracted from the muon lifetime* — so it "predicted" τ_μ from a constant measured out of τ_μ (circular, hence the suspicious 0.000%). It now uses the **substrate** Fermi constant (`electroweak.fermi_constant_substrate`, 0.006% from PDG) and the substrate α (`isa.ALPHA_SUBSTRATE`) — zero measured couplings in (m_μ, m_e remain measured inputs by design, isolating the weak vertex from the mass formula). With the standard phase-space + 1-loop QED layer it lands at **0.011%** vs PDG — a genuine prediction replacing the circular 0.000%, finally making the docstring's "substrate weak-sector closure / α is substrate-predicted" claims true.
- **QED cross sections now close on the *derived* α (route-diversity follow-up).** The v0.4.2 diversity readout flagged five leptonic benchmarks (bhabha, Møller, Compton, e⁺e⁻→μ⁺μ⁻, muon-decay-rate) as single-sector "hidden SPOFs" — graph-2-route but collapsing to one effective route. Root cause: a measured-input leak. The QED vertex coupling `amplitudes.vertices.ELECTRIC_CHARGE` (and the textbook closed forms in `amplitudes.cross_sections`, plus the `qed.constants.alpha` shim) hardcoded the **CODATA** `e = √(4πα)` / `α = 1/137.035999084`, so the substrate's *derived* α never reached the matrix elements (`e⁴`) — only the dimensional integers `8`/`4` did, and those co-move (one Cl(0,7) descent → one sector). All three now source α from `isa.ALPHA_SUBSTRATE` (= 1/(25π√3+1), 7.6 ppm from CODATA → ≤15 ppm on σ, value-preserving). **Effect:** the four pure-QED cross sections now move under α → raw 3, **effective 2** routes (α is an independent sector from the dimensional one), so they're genuinely resilient *and* threshold-robust; hidden SPOFs drop **5 → 1**. (`muon_decay_rate` — a weak-sector closure with its own local α literal — is left as a separate follow-up.) This is the same derived-not-measured discipline o10/predict enforce; the diversity metric doubled as a measured-input-contamination detector.

## [0.4.2] - 2026-05-31

[Full narrative release notes](docs/releases/v0.4.2.md). Version DOI [10.5281/zenodo.20476222](https://doi.org/10.5281/zenodo.20476222); concept DOI [10.5281/zenodo.20012027](https://doi.org/10.5281/zenodo.20012027) (resolves to latest).

This release answers M. Wende's follow-up on the d12rg benchmark thread — that route *multiplicity* is not route *independence* — and with it formally ships the derivation-route-redundancy layer first shared as v0.4.1 on the list, plus two scheme-correct benchmark-comparison fixes carried since v0.4.0. Value-preserving throughout; full suite green (1368 tests).

### Added

- **Derivation-route redundancy — independent routes + single points of failure (M. Wende).** `DerivationDAG` gains the leave-one-route-out dual of the load ranking: `independent_routes()` counts internally node-disjoint derivation routes to a target (Menger, via unit-node-capacity max flow — conjunctive premises funnelling through one closed form count as *one* route; a value reached two disjoint ways counts as many), `cut_nodes()` lists the single points of failure that survive removal, `route_redundancy()` reports both per target plus resilience (routes ≥ 2), and `criticality_ranking()` ranks nodes by *criticality* (outputs ungrounded if the node is removed) beside *load* (outputs reached), with `redundancy = load − critical`. Surfaced by `python -m nwt_substrate.benchmarks.o10 [--suite --sensitivity] --redundancy`. Key result: high reach is not high criticality — `DIM_S_SPIN7` carries the largest load (14 benchmarks) but α is the sole route for six pure-α observables.
- **Derivation-*diversity* layer — route independence beyond route count (M. Wende follow-up).** Route multiplicity is not route independence: two graph-disjoint routes can collapse under the same perturbation. Each integer's mover-set (from the sensitivity sweep) is its empirical perturbation-response signature, so `SensitivityReport` now clusters the load-bearing integers into structural **sectors** by response overlap and re-counts routes per sector: `integer_similarity()` (Jaccard of two integers' mover-sets), `integer_sectors(threshold)` (single-linkage sectors), and `route_diversity(threshold)` (per benchmark: `raw` routes vs `effective` = distinct sectors; `effective < raw` is false redundancy, `effective == 1` with `raw ≥ 2` is a hidden single point of failure). `python -m nwt_substrate.benchmarks.o10 --suite --sensitivity --redundancy` appends the diversity readout. **Finding:** five QED/lepton benchmarks (`bhabha`, `moller`, `muon_decay_rate`, `qed_compton`, `qed_eemumu`) that the route count called 2-route "resilient" are single-sector — `DIM_S_SPIN7` and `N_LORENTZ_FROM_CL07` always co-move (Jaccard 0.53) — while `DIM_V_SPIN7 ≡ N_EDGES_K7` and the K8 trio are perturbation-identical (Jaccard 1.0) at every threshold. The five hidden SPOFs are threshold-contingent (the 0.53 spine merge appears at threshold ≤ 0.53); the geometry- and K8-sector collapses are threshold-robust.

### Fixed

- **`benchmark_sin2_theta_W` — scheme-correct comparison.** The substrate prediction `(2+α)/9 ≈ 0.22303` is a leading-order (tree-level) angle, which *is* the **on-shell** `sin²θ_W ≡ 1 − M_W²/M_Z²` by definition — it matches the PDG on-shell value to **<0.1%** (0.009% at M_W=80.379). The benchmark had been comparing this on-shell prediction against the **effective/MS-bar** angle (0.23122), an apples-to-oranges scheme mismatch that read as a 3.5% deviation; the on-shell↔effective separation (+3.68% radiative running) *is* that gap. Now compared like-for-like against the on-shell value (computed from the PDG `M_W, M_Z` already in `constants.py` — not a hand-picked witness), it moves from a marked O10 defect edge to admissible. The prediction is unchanged; only the comparison scheme is corrected, and the notes keep the effective-angle running explicit as an open higher-order item.

### Changed

- **`benchmark_pmns_angles` — clearer defect framing (note-only, no numbers changed).** Documented that the PMNS leading order is tri-bimaximal + a √(3α) reactor angle: `θ_13 = √(3α) = 8.51°` is a first-principles success (0.7% vs NuFIT 8.57°), while `θ_12`/`θ_23` are the LO tri-bimaximal values whose ~5% deviation has a known but *underived* NLO correction (the α-suppressed charged-lepton `U_ℓ` rotation, Paper 20 §7.6). That correction is deliberately **not** applied — its NuFIT-reproducing magnitudes would be a fit, not a derivation — so this stays an honestly-marked O10 defect edge. Prediction and comparison witness unchanged.

## [0.4.0] - 2026-05-29

[Full narrative release notes](docs/releases/v0.4.0.md). Version DOI [10.5281/zenodo.20448443](https://doi.org/10.5281/zenodo.20448443); concept DOI [10.5281/zenodo.20012027](https://doi.org/10.5281/zenodo.20012027) (resolves to latest).

This release responds to the first independent reviews of the v0.3.0 benchmark suite (P. Kaboth, M. Wende, L. Leighton, on d12rg) and runs a full self-consistency audit of the library. It adds physically-motivated correction layers where benchmarks were leading-order, a structural-criticality layer on the sensitivity sweep, and an O10 derivation-separation predictor + DAG cit-readout — while a library-wide audit anchors every structural integer on `isa` and fixes a handful of internal issues. Value-preserving throughout; full suite green.

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

### Benchmarks — derivation-separation predictor (L. Leighton's protocol)

`nwt_substrate.benchmarks.predict` — a standalone emitter of *dimensionless* substrate predictions (1/α, sin²θ_W, Cabibbo λ, η_B, m_e/M_Pl) that takes **no input of any kind**: α is the substrate closed form `1/(25π√3+1)` and every value is a pure function of the K_7/Spin(7) integers. Measured reference values are quarantined in `REFERENCE` and never read by `predictions()`. The two streams are emitted separately and compared *outside* the derivation:

```
python -m nwt_substrate.benchmarks.predict            > predicted.txt
python -m nwt_substrate.benchmarks.predict --reference > measured.txt
diff -u predicted.txt measured.txt
```

This enforces the prediction/measurement boundary structurally rather than by convention — a measured value cannot leak into a prediction, which is the class of bug the v0.3.1 sin²θ_W fix corrected. (The diff also makes plain that sin²θ_W is a ~3.5% leading-order angle, not a sub-1% result.) Implements the standalone-output rung of L. Leighton's O10 "Ladder Derivation Protocol" (DAG cit-readout specialization): the proof-order edges `analytic proof → symbolic parent → evaluator → standalone output → witness` are one-way, and per O10's constants rule the witness layer is CODATA-2018 / PDG (post-SI2019 *defined* constants excluded even as witnesses).

`nwt_substrate.benchmarks.o10` — the whole DAG that `predict.py`'s rung lives in, as a validated structure. Models the constants stack as a directed acyclic graph (stages `STRUCTURAL → SYMBOLIC → EVALUATOR → OUTPUT → WITNESS`) and checks the O10 invariants on it: edges are one-way (`backward_edges()`), the graph is acyclic, witnesses are sinks, and the multi-path *commutative-diagram* identities agree (grounded in isa's import-time asserts, e.g. `21 = C(7,2) = 3·7`). `cit_readout()` is the transversal witness-invariance check — and it *marks* the one defect edge (sin²θ_W, the LO angle) rather than repairing it. `load_ranking()` (which nodes reach the most observables — α and its closed form reach all five) and `horn_frontier()` (the `(s₁ ∧ … ∧ sₙ) → W` premise set) connect O10's multi-parent readout to M. Wende's structural-load idea. `python -m nwt_substrate.benchmarks.o10` prints the acceptance checklist, cit readout, load ranking, and commutativity checks.

`build_suite_dag()` extends the DAG to the **whole 38-benchmark suite** as one O10-validated graph: each benchmark is an OUTPUT node with a WITNESS edge carrying its deviation (parsed from `substrate_accuracy` — ranges take the conservative upper bound; benchmarks with no vs-measurement metric, e.g. the cross-section computations, are flagged *qualitative* rather than fabricated). The cit readout is then the structurally-enforced falsification report Pasquale asked for: 27 admissible, 8 defect edges *marked* (the known >1% benchmarks — sin²θ_W, the muon compound, PMNS, mass spectrum, …), 3 qualitative. The STRUCTURAL→OUTPUT edges are the **computed** sensitivity coupling: pass a `SensitivityReport` (or `--sensitivity`) and `integer → benchmark` edges are drawn for every integer that moves that benchmark, so the DAG's `load_ranking()` *is* the sweep's structural load. `python -m nwt_substrate.benchmarks.o10 --suite [--sensitivity]`.

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

[Unreleased]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.4.2...HEAD
[0.4.2]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.4.0...v0.4.2
[0.4.0]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/JimGalasyn/nwt-substrate/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/JimGalasyn/nwt-substrate/releases/tag/v0.1.0
