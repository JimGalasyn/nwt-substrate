# `nwt_substrate.benchmarks`

[![benchmarks](https://github.com/JimGalasyn/nwt-substrate/actions/workflows/benchmarks.yml/badge.svg)](https://github.com/JimGalasyn/nwt-substrate/actions/workflows/benchmarks.yml)

Substrate algebra vs traditional methods: speed and accuracy comparison
across **38 physical observables** spanning particle physics, atomic
physics, QED/QCD, electroweak precision, cosmology, gravity, black-hole
thermodynamics, and chemistry. Full suite runs in **~110 ms**.

## What this is

The substrate program's anti-numerology argument made concrete and
reproducible. Each benchmark times a substrate-algebra closed-form
prediction, compares it against the PDG/CODATA/Planck observed value,
and quotes the rough computational cost of the traditional method that
would otherwise be required (lattice QCD, CKMfitter, LEP-1 precision
program, etc.).

If you can reproduce this output on your hardware, you have empirical
evidence that the substrate algebra is genuinely predictive — not
fitted, not numerology. **38 substantive physical observables in
~100 ms of CPU time**, zero free parameters tuned to any single one.

## Quick start

```python
from nwt_substrate.benchmarks import run_all
run_all()
```

This prints a comparison table for each of the 38 benchmarks plus a
summary line at the bottom.

For a single benchmark:

```python
from nwt_substrate.benchmarks import benchmark_higgs_vev
result = benchmark_higgs_vev()
print(result.substrate_value)            # v_EW = 246.2128 GeV vs PDG 246.2197 GeV
print(result.substrate_accuracy)         # 27.7 ppm vs PDG
print(result.substrate_time_us)          # ~1 microsecond
```

## What's covered

### Fundamental couplings (5 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Fine structure constant α | 7.6 ppm |
| Newton's G | 11 ppm |
| Fermi constant G_F | 55 ppm |
| Weak mixing angle sin²θ_W (on-shell, (2+α)/9) | <0.1% vs `1−M_W²/M_Z²` (effective angle +3.68% via running) |
| Strong coupling α_s(M_Z), Λ_QCD, Λ_χ | matches PDG |

### Higgs sector (2 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Higgs VEV v_EW | 27.7 ppm |
| Higgs mass m_h via λ_H = 18α | 0.9% (also predicts a second scalar at 98 GeV) |

### Mass spectra (4 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Particle compendium (25 particles) | ~1% median |
| Neutrino sector (3 active + 3 sterile + δ_CP) | 0.04% on ν₁ |
| K_8 dark matter mass tower (11 rungs) | ~0.1% on anchored rungs |
| Vector meson decay constants (11 states) | 1-2% |

### Flavor mixing (3 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Cabibbo angle θ_C | ~0.1% |
| Full CKM matrix (V_us, V_cb, V_ub, V_td, J) | ~1% |
| PMNS leptonic angles | θ_13 = √(3α) at 0.7%; θ_12/θ_23 ~5% (LO tri-bimaximal) |

### Decay constants and rates (3 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Pseudoscalar decay constants (f_π, f_K, f_η, f_D, f_Ds, f_B, f_Bs) | 1-3% |
| Z boson width + lepton universality | 0.31% (with hadronic QCD correction) / **0.9 ppm** |
| Muon lifetime — weak-sector closure (substrate G_F + PDG m_μ + SM correction) | 0.007% (the m⁵-amplified compound benchmark is a separate mass-formula probe) |

### Atomic physics + QED (3 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Hydrogen chain (a₀, R_H, Lyman α, 21 cm, Lamb) | ~7 ppm on a₀ and Lyman α |
| Electron magnetic moment a_e (Schwinger 1-loop) | matches Schwinger formula exactly |
| QED Compton scattering (Thomson limit) | 53 ppm |
| QED e⁺e⁻ → μ⁺μ⁻ at LEP2 (200 GeV) | LO QED at substrate α |

### Cosmology (3 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Cosmological constant Λ | 0.74% (solves 123-OoM problem) |
| Baryon/CDM ratio Ω_b/Ω_c | **0.0067%** (better than Planck systematic!) |
| Baryon asymmetry η_B | 0.38% (no traditional method even attempts it) |

### Gravity and black-hole thermodynamics (2 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Black hole thermodynamics (T_H, r_S, τ_evap) | closed form using substrate G |
| Cosmogenesis (κ_parent, f_J, Thorne a*=0.998) | matches Bardeen-Press-Teukolsky exact |

### Topological / algebraic foundations (2 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| K_7 Heffter genus-1 toroidal embedding | exact (V=7, E=21, F=14) |
| SU(2)_5 modular tensor category data (6 anyons + S/T + c=15/7) | exact closed form |

### Chemistry (3 benchmarks)

| Observable | Substrate accuracy |
|---|---|
| Aromaticity + NICS + C_60 combinatorics | 100% on Hückel/Möbius set, 14/14 NICS, C_60 exact |
| NMR chemical shift sign rule | 14/14 (Hopf-pair parity, O(1) lookup) |
| C_60 vibrational mode decomposition (174 modes) | exact (group theory) |

### O10 derivation-separation & structural criticality (v0.4.0)

Added in response to the d12rg review round (L. Leighton, M. Wende):

- **`python -m nwt_substrate.benchmarks.predict`** — L. Leighton's O10
  derivation-separation rung: emits *only* substrate-derived dimensionless
  predictions (no input of any kind), with measured values quarantined in a
  separate `--reference` stream (CODATA-2018 witness; post-SI2019 defined
  constants excluded). `diff -u` of the two streams is the falsification report —
  a measured value structurally cannot leak into a prediction.
- **`python -m nwt_substrate.benchmarks.o10 --suite`** — the whole 38-benchmark
  suite as one validated O10 DAG: one-way proof-order edges, witness sinks,
  commutative-diagram identities (= isa import-time asserts, e.g. `21=C(7,2)=3·7`),
  and a cit readout with defect edges **marked** (not repaired). With
  `--sensitivity`, the `STRUCTURAL→OUTPUT` edges are the *computed* coupling, so
  the DAG's load ranking equals the sweep's (the sweep also perturbs the derived
  scalar α and κ, so α-anchored benchmarks couple to their root). Add
  **`--redundancy`** for M. Wende's leave-one-route-out layer: how many
  *independent* routes converge on each answer, the single points of failure,
  and node **criticality** (outputs ungrounded if a node is removed) vs **load**
  (outputs reached). It separates the two — α is the master root (load 19, sole
  route for 6 pure-α observables like the electron anomaly), while `DIM_S_SPIN7`
  (=8) is high-*load* (14 benchmarks) but, now that α backs up most of them, the
  sole route for only one (redundancy 13).
- **`python -m nwt_substrate.sensitivity --criticality`** — M. Wende's
  structural-criticality layer: which structural knobs carry the most load. The
  sweep perturbs the ISA integers *and* the derived scalars α and κ; α (the
  master coupling) reaches 19 of 38 benchmarks — more than any integer
  (`DIM_S_SPIN7` moves 14) — and the gravity/cosmology benchmarks co-move.

## Methodology

Each benchmark returns a `BenchmarkResult` dataclass with:

- `name`: short description
- `substrate_time_us`: wall-clock μs of substrate calculation
- `substrate_value`: what was computed (with units)
- `substrate_accuracy`: comparison to PDG/CODATA/Planck
- `traditional_method`: what method would otherwise be used
- `traditional_cost`: order-of-magnitude effort of the traditional method
- `speedup_factor_str`: rough speedup multiplier
- `notes`: physics context, what the substrate prediction means

## Honest framing

The substrate algebra's computational advantage is for **forward
prediction at substrate accuracy** (typically 0.1-7% on observables;
ppm on a handful of dimensionful constants). Traditional methods
retain the precision lead for sub-ppm measurement validation.

**Best use**:
- Substrate algebra → forward prediction, parameter reduction
- Traditional methods → precision validation of substrate predictions
- L_NWT (Paper 16) → collider phenomenology, scattering amplitudes
- H_SU(2)_5/K_7 (D43 chain) → topological-phase verification,
  condensed-matter analogs

The benchmarks demonstrate that the substrate algebra is a **genuine
predictive framework** for forward inference across particle, atomic,
and cosmological physics — and that "numerology" objections fail
empirically: a 2-integer (18, 25) framework correctly predicting 32
independent physical observables in 100 ms of arithmetic is not
pattern matching.

## Adding new benchmarks

Follow the pattern in `compute_speed.py`:

```python
def benchmark_my_new_observable() -> BenchmarkResult:
    """One-line description."""
    from nwt_substrate.<module> import <substrate_function>

    t0 = time.perf_counter_ns()
    value = substrate_function()
    elapsed_us = (time.perf_counter_ns() - t0) / 1e3

    PDG_VALUE = ...
    err_pct = abs(value - PDG_VALUE) / PDG_VALUE * 100

    return BenchmarkResult(
        name="…",
        substrate_time_us=elapsed_us,
        substrate_value=f"{value:.4f} vs PDG {PDG_VALUE:.4f}",
        substrate_accuracy=f"{err_pct:.2f}% vs PDG",
        traditional_method="…",
        traditional_cost="…",
        speedup_factor_str="~10^N× ()",
        notes="What this is, why it matters.",
    )
```

Then add to the `benches` list in `run_all()` and to the `__init__.py`
exports.

## Tests

The unit tests (`nwt_substrate/tests/test_benchmarks.py`) verify each
benchmark returns the expected structure and characteristic values
(incl. the O10 suite-DAG invariants and the sin²θ_W on-shell comparison).
Run with:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest nwt_substrate/tests/test_benchmarks.py -v
```

## References

- Paper 14: Newton's G from K_7 Wilson amplitude
- Paper 16: NWT Lagrangian (L_NWT, soliton sector)
- Paper 17: 5-loop α and electron mass ratio
- Paper 20: K_8 mass tower and neutrino sector
- Paper 21a/21b: cross-vendor / cross-architecture experimental verification
- VV's stage-7 substrate-EFT matching (substrate-DNA integers 18, 25)
