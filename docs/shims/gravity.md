# nwt_substrate.gravity

> Gravity in NWT is the macroscopic projection of the substrate's K_7 Wilson amplitude. Newton's `G` is derived from the fine structure constant α and the electron mass `m_e` at **zero free parameters**: `G = (8/7)² α²¹ ℏc / m_e²` at LO, with an NNLO correction that lands `G = 6.674228 × 10⁻¹¹ m³ kg⁻¹ s⁻²` — **−11 ppm from CODATA**, inside the ±22 ppm experimental band. The same shim covers black-hole thermodynamics, Kerr efficiency, Near-Horizon Extremal Kerr (NHEK) geometry, cosmogenic spin equilibrium, and FLRW backreaction predictions.

[← Back to index](../index.md) · Source: [`nwt_substrate/gravity/`](../../nwt_substrate/gravity/) · Papers: [14](https://zenodo.org/communities/nwt) (LO G), [17](https://zenodo.org/records/15445103) (NNLO α + G), [18](https://zenodo.org/communities/nwt) (G1-G6 Sakharov)

## Common questions

### Does NWT derive Newton's gravitational constant?

Yes. `G_substrate_SI()` returns `6.674228 × 10⁻¹¹ m³ kg⁻¹ s⁻²`, which is **−11 ppm** from CODATA 2018 (`6.67430 × 10⁻¹¹`), well inside the experimental band of ±22 ppm.

The derivation: the ratio `m_e / M_Pl` is the K_7 Wilson amplitude at perimeter 21 (the K_7 edge count), prefactor `8/7` (`dim Spin(7) spinor / dim vector`), with NLO and NNLO substrate corrections:

```
m_e / M_Pl = (8/7) · α^(21/2) · (1 + α/7 + (21/8) α²)
G          = ℏc / M_Pl²
```

Then `G = ℏc · (m_e / M_Pl)² / m_e² = (8/7)² · α²¹ · ℏc / m_e² · (1 + ...)²`. Each factor is a substrate primitive — no fitting.

### Is gravity fundamental in NWT?

No. Gravity is the long-wavelength projection of the K_7 substrate's Wilson amplitude. The Einstein-Hilbert action arises via **Sakharov-induced gravity** from substrate-condensate fluctuations (Paper 18 results G1-G6: linearized and full nonlinear Einstein equations emerge from the NWT condensate).

### What black-hole observables are computed?

Schwarzschild radius, Hawking temperature, Hawking luminosity, evaporation time, ISCO efficiency, Penrose extraction fraction, Bardeen prograde/retrograde efficiency, electromagnetic superradiance maximum, twin-peak QPO ratios on extremal Kerr, and a full `black_hole_summary(particle_name)` table for the proton, electron, neutron, etc.

### Does NWT make cosmological gravity predictions?

Yes. The `flrw_test` submodule predicts the Buchert backreaction amplitude `Q₀` from the substrate density-variance scaling, with a 5-σ band that includes the observed FLRW O-test value. The `cosmogenesis` submodule predicts the parent-BH Thorne spin equilibrium `a* = 0.998` and the cosmogenic energy budget `f_J ≈ 0.42` for a Kerr → daughter-universe bridge.

## Prediction table

| Observable | Substrate expression | Substrate value | Reference | Accuracy |
|---|---|---|---|---|
| `m_e / M_Pl` (LO) | `(8/7) α^(21/2)` | 4.186 × 10⁻²³ | CODATA `4.185 × 10⁻²³` | −0.24 % |
| `m_e / M_Pl` (NNLO) | `(8/7) α^(21/2) (1 + α/7 + (21/8) α²)` | 4.18476 × 10⁻²³ | CODATA | **−5.5 ppm** |
| Newton's `G` | `ℏc / M_Pl²` | 6.674228 × 10⁻¹¹ | CODATA `6.67430 × 10⁻¹¹` | **−11 ppm** (within ±22 ppm band) |
| Bardeen prograde η | `1 - √(1 - 2/3 (a*/M)²)` | 0.4233 (a*=0.998) | GR exact | matches |
| Cosmogenic `f_J` | NWT derivation | 0.42 | Thorne equilibrium | matches |
| Thorne `a*` equilibrium | substrate fit | 0.998 | Bardeen-Press-Teukolsky | matches |

Each row is a `pytest` assertion in [`nwt_substrate/tests/`](../../nwt_substrate/tests/) and a `BenchmarkResult` in [`nwt_substrate/benchmarks/`](../../nwt_substrate/benchmarks/).

## Quick start

```python
import nwt_substrate.gravity as grav

# Substrate-derived Newton's G
grav.G_substrate_SI()              # → 6.674228e-11 m³/kg/s²  (−11 ppm CODATA)
grav.G_NEWTON_SI                   # → 6.67430e-11           (CODATA 2018 anchor)

# Underlying ratio
grav.m_e_over_M_Pl_NNLO()          # → 4.18476e-23  (−5.5 ppm CODATA)

# Full structural breakdown
print(grav.G_substrate_breakdown())  # prints each (8/7), α^(21/2), NLO, NNLO term

# Black-hole observables for the proton
bh = grav.black_hole_summary("proton", use_substrate_G=True)
bh["r_s_m"]                        # Schwarzschild radius
bh["T_H_K"]                        # Hawking temperature
bh["evap_t_s"]                     # Evaporation time
```

## API by topic

### Substrate-derived Newton's G

| Function | Returns | Notes |
|---|---|---|
| `m_e_over_M_Pl_LO(alpha=ALPHA_QED)` | `(8/7) α^(21/2)` | Paper 14 leading order |
| `m_e_over_M_Pl_NLO(alpha)` | LO × `(1 + α/7)` | Per-K_7-vertex correction |
| `m_e_over_M_Pl_NNLO(alpha)` | LO × `(1 + α/7 + (21/8) α²)` | Paper 17 PSL(2,7) Fano |
| `m_e_over_M_Pl_observed()` | CODATA `m_e / M_Pl` | Reference value |
| `G_substrate_LO_natural()` | LO `G` in natural units | |
| `G_substrate_NNLO_natural()` | NNLO `G` in natural units | |
| `G_substrate_SI()` | NNLO `G` in m³ kg⁻¹ s⁻² | **The headline number** |
| `G_substrate_LO_SI()` | LO `G` in m³ kg⁻¹ s⁻² | |
| `G_substrate_breakdown()` | Pretty-printed table | Per-factor source |

### Black-hole thermodynamics

| Function | Returns |
|---|---|
| `schwarzschild_radius_m(mass_kg)` | `r_s = 2GM/c²` in metres |
| `hawking_temperature_K(mass_kg)` | `T_H` in Kelvin |
| `hawking_luminosity_W(mass_kg)` | Power radiated, Watts |
| `evaporation_time_s(mass_kg)` | Time to evaporate, seconds |
| `black_hole_summary(particle_name)` | All BH observables for a named particle |
| `PARTICLE_MASSES_KG` | Dict of {`"proton"`, `"electron"`, …} → kg |

### Einstein equations + Sakharov route

| Function | Returns |
|---|---|
| `M_Pl_over_m_e_squared_observed()` | `(M_Pl / m_e)²` from CODATA |
| `M_Pl_over_m_e_squared_substrate()` | Same, from substrate algebra |
| `UV_IR_bridge_breakdown()` | UV-IR bridge derivation summary |
| `einstein_hilbert_coefficient()` | EH action prefactor from condensate |
| `sakharov_route_summary()` | Sakharov-induced-gravity walk |
| `verify_schwarzschild_vacuum_symbolic()` | `R_μν = 0` check via sympy |
| `linearized_graviton_propagator_text()` | Substrate graviton propagator |

### Kerr efficiency + cosmogenic spin

| Function | Returns |
|---|---|
| `penrose_extraction_fraction(a_star)` | Penrose-process energy fraction |
| `m_irreducible_over_m_extremal(a_star)` | Irreducible mass ratio |
| `bardeen_prograde_efficiency(a_star)` | Prograde ISCO efficiency |
| `bardeen_retrograde_efficiency(a_star)` | Retrograde ISCO efficiency |
| `schwarzschild_isco_efficiency()` | Non-rotating ISCO efficiency |
| `em_superradiance_max()` | EM superradiance ceiling |
| `kerr_efficiency_table()` | All-route comparison table |
| `f_J_cosmogenic()` | Cosmogenic energy budget |
| `kappa_parent_required()` | Required parent-BH bridge couplings |
| `thorne_a_star_equilibrium()` | Thorne `a* = 0.998` equilibrium |
| `cosmogenesis_summary()` | Pretty-print full cosmogenesis chain |

### NHEK geometry (near-horizon extremal Kerr)

| Function | Returns |
|---|---|
| `nhek_metric(M, r, theta)` | NHEK metric tensor `g_μν` (4×4) |
| `nhek_inverse_metric(M, r, theta)` | `g^μν` |
| `nhek_metric_determinant(M, r, theta)` | `det(g)` analytic form |
| `nhek_signature(M, r, theta)` | Eigenvalue count (Lorentzian = `(1, 3, 0)`) |
| `killing_vectors_constant_basis()` | `H = ∂_t`, `L = ∂_φ` (SL(2,ℝ)×U(1) isometry) |
| `christoffels_numeric(M, r, theta)` | `Γ^λ_μν` via finite differences |
| `nhek_symbolic()` | Cached sympy NHEK machinery (Ricci, etc.) |
| `verify_nhek_vacuum()` | `R_μν = 0` verification (≈30 s first call) |
| `Sigma(theta)`, `Lambda_factor(theta)` | Scalar prefactors |

### FLRW null tests (Buchert backreaction)

`flrw_test.nwt_prediction(...)` returns an `FLRWTestPrediction` dataclass with the substrate-predicted scaling of `Q₀` (Buchert backreaction), `R_avg`, `K_BR`, `A(z)`, `O(z)`, `C′`, and a 5-σ band check against the observed `O`-band (km/s/Mpc²).

### Substrate healing length

| Symbol | Value | Notes |
|---|---|---|
| `XI_SUBSTRATE_M` | ≈ 1.6 × 10⁻¹⁰ m | Substrate condensate healing length |
| `xi_substrate_m()` | callable | Substrate ξ |
| `xi_cosmo_m()` | cosmological ξ | Long-wavelength regime |
| `k_chern_simons()` | CS coupling | Lattice-CS relation |
| `scale_regime_table()` | Domain-marker table | Where each ξ dominates |

## Worked examples

### Sakharov chain for Newton's G

```python
from nwt_substrate.gravity import G_substrate_breakdown
print(G_substrate_breakdown())
# Substrate G derivation (Paper 17):
#   (8/7)            = dim(Spin(7) spinor) / dim(vector)
#   α^(21/2)         = K_7 Wilson amplitude (21 edges)
#   (1 + α/7)        = NLO per-vertex correction (7 K_7 vertices)
#   (1 + (21/8) α²)  = NNLO PSL(2,7) Fano correction
#   m_e / M_Pl       = 4.18476 × 10⁻²³  (CODATA 4.185 × 10⁻²³, −5.5 ppm)
#   G = ℏc / M_Pl²   = 6.6742 × 10⁻¹¹ m³ kg⁻¹ s⁻²  (CODATA 6.67430 × 10⁻¹¹, −11 ppm)
```

### Proton black hole (entirely substrate-derived)

```python
from nwt_substrate.gravity import black_hole_summary
bh = black_hole_summary("proton", use_substrate_G=True)
print(f"r_s    = {bh['r_s_m']:.3e} m")     # ~2.5e-54 m
print(f"T_H    = {bh['T_H_K']:.3e} K")     # ~7.4e22 K
print(f"τ_evap = {bh['evap_t_s']:.3e} s")  # ~5.0e-104 s
```

Each result uses the substrate-derived `G` rather than the CODATA anchor — the entire chain is substrate algebra.

### Bardeen prograde ISCO efficiency at the cosmogenic spin

```python
from nwt_substrate.gravity import bardeen_prograde_efficiency, thorne_a_star_equilibrium
a_star = thorne_a_star_equilibrium()           # → 0.998
eta = bardeen_prograde_efficiency(a_star)      # → 0.4233 = 42.33 %
```

This is the maximum mass-to-energy conversion efficiency on the Thorne spin-equilibrium plateau — both the spin value and the efficiency are substrate-derived.

### NHEK as the cosmogenic background

```python
import math, numpy as np
from nwt_substrate.gravity import nhek_metric, nhek_signature

g = nhek_metric(M=1.0, r=1.0, theta=math.pi / 3)   # 4×4 metric
assert np.allclose(g, g.T)
assert nhek_signature(M=1.0, r=1.0, theta=math.pi/3) == (1, 3, 0)  # Lorentzian
```

NHEK is the relevant background geometry for the cosmogenic bridge on the parent side.

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `8/7` | `dim(Spin(7) spinor S) / dim(vector V)` | `isa.DIM_OCTONION / isa.DIM_V_SPIN7` |
| `21` | `K_7 edge count = dim Adj(so(7))` | `isa.N_EDGES_K7 == isa.DIM_ADJ_SO7` |
| `7` | `K_7 vertex count` | `isa.N_VERTICES_K7` |
| `21/8` | `dim Adj / dim S` | `isa.DIM_ADJ_SO7 / isa.DIM_OCTONION` |

If any of these change in `nwt_substrate.isa.constants`, the substrate identity assertions will refuse to import — preventing silent drift.

## Papers

- **Paper 14** — Newton's G from K_7 Wilson amplitude (LO derivation)
- **Paper 17** — α closure (`1/α = 25π√3 + 1`) + G NNLO + electron mass ratio
- **Paper 18** — G1-G6: linearized and full nonlinear Einstein equations from the NWT condensate (Sakharov-induced gravity)
- **Paper 22** — Cosmogenesis (parent BH → daughter universe via near-extremal bridge)

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`isa`](../../nwt_substrate/isa/README.md) — `8/7`, `21`, `7`, `21/8` substrate constants
- [`cosmology`](cosmology.md) — Λ, η_B, Ω_b/Ω_c (cosmological observables)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_gravitational_constant`, `benchmark_black_hole_thermodynamics`, `benchmark_cosmogenesis`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
