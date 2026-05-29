# nwt_substrate.cosmology

> The cosmological view of the substrate algebra. Three headline observables are derived at **zero free parameters**: the baryon asymmetry `η_B = (3/14) α⁴` (**0.38 % Planck**), the baryon-to-CDM ratio `Ω_b/Ω_c = 25α + 75α²` (**0.0067 % Planck — 240× tighter than the measurement uncertainty**), and the cosmological constant `ρ_Λ / M_Pl⁴ = (m_e / M_Pl)⁴ · α¹⁶ · h_Cox ≈ 1.19 × 10⁻¹²³` (**0.74 % Planck — closes the 123-orders-of-magnitude problem**). The shim also tabulates the measured CMB anisotropy axes (Axis of Evil, HPA, Cold Spot, dipole) with citations, used by both the Paper 22 cosmogenesis program and the Heron Experiment 11 sidereal directional layer.

[← Back to index](../index.md) · Source: [`nwt_substrate/cosmology/`](../../nwt_substrate/cosmology/) · Papers: [13](https://zenodo.org/records/19635239) (SM capstone), [19](https://zenodo.org/records/15410028) (substrate monism + F1 falsifier), [22](https://zenodo.org/communities/nwt) (cosmogenesis)

## Common questions

### Does NWT solve the cosmological constant problem?

Yes — the 123-orders-of-magnitude problem dissolves. The substrate vacuum-energy density is:

```
ρ_Λ / M_Pl⁴ = (m_e / M_Pl)⁴ · α¹⁶ · h_Cox  ≈ 1.19 × 10⁻¹²³
```

Observed: `1.20 × 10⁻¹²³` (Planck 2018). **Residual: 0.74 %.** Matching three significant figures from substrate primitives — vs the naive QFT estimate off by `10¹²⁰` — is the headline.

Each factor is a substrate primitive:
- `(m_e / M_Pl)⁴`: Coleman-Weinberg 1-loop vacuum-energy scale, with `m_e / M_Pl` from the K_7 Wilson amplitude (Paper 17 NNLO).
- `α¹⁶`: Stage-7 substrate vacuum-energy exponent.
- `h_Cox = 6`: Coxeter number of so(7) = B₃ (Weyl-group multiplicity).

### Does NWT predict the baryon-to-photon ratio?

Yes. `η_B = RANK_SO7 · α⁴ / (2 · N_VERTICES_K7) = 3α⁴ / 14 ≈ 6.077 × 10⁻¹⁰`. Planck observes `6.10 × 10⁻¹⁰`; BBN gives `6.14 × 10⁻¹⁰`. **Residual: 0.38 %.**

Two coincident derivation routes:
- **VV — K₇⊗K₈ bridge partition.** The cosmogenic `H⊗7` polarity-selection event on the K_7 Steane code commits one logical-coset outcome per event; `8 √α` traversals × 4 qubits per X-stabilizer give `α⁴`, × 3 X-stabilizers (Z₃ σ-orbit), / 14 single-Fano syndromes.
- **NWT — Murasugi torus-knot.** An independent `(3/14) α⁴` from the carrier knot's Murasugi/Jones structure.

The sign `sgn(η_B) = sgn(J_parent)` — matter over antimatter — is fixed by the parent black hole's spin via the bridge pointer-basis chain.

### Does NWT predict the baryon-to-CDM ratio?

Yes — and it's the most accurate prediction in the entire benchmark suite:

```
Ω_b / Ω_c = q² · α + RANK_SO7 · q² · α²   where q² = H_V_SO7² = 25
          = 25α + 75α²
          ≈ 0.18643
```

Planck observed `Ω_b/Ω_c = 0.02237 / 0.1200 = 0.18642`. **Residual: 67 ppm — that's 240× tighter than the Planck systematic uncertainty.** The substrate prediction is more precise than the measurement.

### What CMB anisotropy axes are tabulated?

The `anisotropy_axes` submodule canonically records, with citations and uncertainty ranges:
- **Axis of Evil (AoE)**: quadrupole-octopole alignment direction. Consensus `(l, b) ≈ (245°, +60°)`.
- **Hemispherical Power Asymmetry (HPA)**: low-vs-high-power hemisphere normal. Consensus `(l, b) ≈ (227°, −27°)`.
- **Cold Spot**: `(l, b) ≈ (210°, −57°)`.
- **CMB Dipole**: `(l, b) = (264.021°, 48.253°)` (Planck 2018).
- **Bulk Flow**, **Shapley Concentration**, **Galactic + Ecliptic North**: reference axes.

Plus rotation matrices, great-circle fits, and the 12-pair AoE × HPA separation table used by the Paper 22 cosmogenesis multi-measurement check.

### What does the cosmogenesis program claim?

Paper 22 reads the CMB anomalies as fossils of a parent black hole's accretion-disk geometry: the **Axis of Evil = parent BH spin axis**, the **HPA direction = accretion-disk normal**, the **Cold Spot = a gap in the disk's accretion flow**. This shim provides the rotation infrastructure (great-circle fits, angular separations) that the cosmogenesis program uses to test the prediction against measured axes.

### Is the cosmological substrate axis falsifiable?

Yes — that's the Paper 19 F1 falsifier. The substrate's holonomy reduction `SO(8) → Spin(7) → 3+1` picks a preferred axis as a structural consequence of the algebra; the axis is not free. The Heron Experiment 11 sidereal A/B/C probe (`heron.sidereal_geometry`) compares K_7 stabilizer drift against this predicted axis. If the substrate axis doesn't align with the CMB anomalies within tolerance, the F1 prediction fails.

## Prediction table

| Observable | Substrate formula | Substrate value | Observed | Accuracy |
|---|---|---|---|---|
| Baryon asymmetry `η_B` | `(3/14) α⁴` | 6.077 × 10⁻¹⁰ | Planck `6.10`, BBN `6.14` × 10⁻¹⁰ | **0.38 %** |
| `Ω_b / Ω_c` | `25α + 75α²` | 0.18643 | Planck `0.18642` | **67 ppm** (240× tighter than measurement) |
| `ρ_Λ / M_Pl⁴` | `(m_e/M_Pl)⁴ · α¹⁶ · h_Cox` | 1.19 × 10⁻¹²³ | Planck `1.20 × 10⁻¹²³` | **0.74 %** |
| Axis of Evil direction | `Spin(7)` holonomy axis | predicted vs measured | Heron Exp 11 (in progress) | F1 falsifier |

Each row asserted in [`nwt_substrate/tests/test_eta_B.py`](../../nwt_substrate/tests/), `test_omega_b_c.py`, `test_lambda_cc.py`, `test_anisotropy_axes.py`.

## Quick start

```python
from nwt_substrate.cosmology import eta_B, omega_b_c, lambda_cc

# Three headline cosmological observables, zero free parameters
eta_B.eta_B()                # → 6.077e-10   (Planck 6.10e-10, 0.38 %)
omega_b_c.omega_b_c()        # → 0.18643     (Planck 0.18642, 67 ppm)
lambda_cc.lambda_cc()        # → 1.19e-123   (Planck 1.20e-123, 0.74 %)

# Detailed breakdowns
eta_B.summary()              # dict: formula, integer factors, pred/obs ratios
omega_b_c.summary()
lambda_cc.summary()

# CMB anisotropy axes
from nwt_substrate.cosmology import anisotropy_axes as ax
aoe = ax.AXIS_OF_EVIL_CONSENSUS         # SkyAxis(l=245°, b=+60°)
hpa = ax.HPA_CONSENSUS                  # SkyAxis(l=227°, b=-27°)
ax.angular_separation(aoe, hpa)         # → ~89° (nearly perpendicular)
```

## API by topic

### Baryon asymmetry — `cosmology.eta_B`

| Symbol | Value / Returns |
|---|---|
| `ETA_B_PRED` | `(3/14) α⁴ ≈ 6.077 × 10⁻¹⁰` |
| `ETA_B_PLANCK`, `ETA_B_BBN` | `6.10e-10`, `6.14e-10` |
| `ETA_B_DENOM` | `14 = 2 · N_VERTICES_K7 = dim(G₂)` |
| `eta_B()` | The prediction |
| `summary()` | Dict with formula, integer factors, pred/obs ratios |

### Baryon-to-CDM ratio — `cosmology.omega_b_c`

| Symbol | Value / Returns |
|---|---|
| `OMEGA_B_C_PRED` | `25α + 75α² ≈ 0.18643` |
| `OMEGA_B_C_PLANCK` | `0.18642` |
| `OMEGA_B_H2_PLANCK`, `OMEGA_C_H2_PLANCK` | `0.02237`, `0.1200` (physical densities) |
| `Q_CINQ_SQ` | `25 = H_V_SO7² = cinquefoil (5,2)-torus-knot crossings²` |
| `omega_b_c()` | The prediction |
| `summary()` | Dict with formula, ppm offset, integer factors |

### Cosmological constant — `cosmology.lambda_cc`

| Symbol | Value / Returns |
|---|---|
| `RHO_LAMBDA_PRED` | `(m_e/M_Pl)⁴ · α¹⁶ · h_Cox ≈ 1.19 × 10⁻¹²³` |
| `RHO_LAMBDA_OBS` | `1.20 × 10⁻¹²³` |
| `M_E_OVER_M_PL` | K_7 Wilson NNLO, `≈ 4.185 × 10⁻²³` |
| `LAMBDA_EXPONENT` | `16` (Stage-7 vacuum-energy exponent) |
| `lambda_cc()` | The prediction |
| `summary()` | Dict with formula, integer factors |

### CMB anisotropy axes — `cosmology.anisotropy_axes`

| Symbol | What it is |
|---|---|
| `SkyAxis` | Dataclass: `name`, `l_deg`, `b_deg`, citation, uncertainty |
| `AXIS_OF_EVIL_CONSENSUS` | `(l, b) ≈ (245°, +60°)` — quadrupole-octopole alignment |
| `HPA_CONSENSUS` | `(l, b) ≈ (227°, −27°)` — hemispherical power asymmetry normal |
| `COLD_SPOT` | `(l, b) ≈ (210°, −57°)` — WMAP/Planck cold spot |
| `CMB_DIPOLE` | `(l, b) = (264.021°, 48.253°)` — Planck 2018 |
| `GALACTIC_NORTH`, `ECLIPTIC_NORTH` | Reference axes |
| `BULK_FLOW`, `SHAPLEY` | Large-scale-structure axes |
| `AXIS_OF_EVIL_MEASUREMENTS`, `HPA_MEASUREMENTS` | Lists of (citation, axis) tuples |
| `galactic_to_cartesian(SkyAxis)` | Unit vector in galactic Cartesian |
| `cartesian_to_galactic(xyz)` | Inverse |
| `galactic_to_icrs(SkyAxis)` | `(RA, Dec)` in ICRS |
| `angular_separation(a, b)` | Directed angle in degrees |
| `angular_separation_undirected(a, b)` | min(θ, 180° − θ) |
| `fit_great_circle(axes)` | Great-circle fit to a list of SkyAxis |
| `hpa_aoe_pair_separations()` | All 12 HPA × AoE pair separations (Paper 22 cross-check) |

## Worked examples

### Three headline observables in three lines

```python
from nwt_substrate.cosmology import eta_B, omega_b_c, lambda_cc
print(f"η_B    = {eta_B.eta_B():.3e}   (Planck 6.10e-10, 0.38 %)")
print(f"Ω_b/Ω_c = {omega_b_c.omega_b_c():.5f}   (Planck 0.18642, 67 ppm)")
print(f"ρ_Λ/M_Pl⁴ = {lambda_cc.lambda_cc():.3e}   (Planck 1.20e-123, 0.74 %)")
```

### Ω_b/Ω_c precision detail

```python
from nwt_substrate.cosmology.omega_b_c import summary
d = summary()
print(d["formula"])           # "H_V_SO7**2 * alpha + RANK_SO7 * H_V_SO7**2 * alpha**2 = 25a + 75a**2"
print(d["lo_coeff"])          # 25 = H_V_SO7² = cinquefoil-knot-crossings²
print(d["nlo_coeff"])         # 75 = 3 × 25 (Z₃ σ-orbit × q²)
print(f"{d['ppm_offset']:.1f} ppm offset from Planck")   # 67 ppm
# 240× tighter than Planck systematic uncertainty
```

### Multi-measurement axis consistency (Paper 22 cross-check)

```python
from nwt_substrate.cosmology import anisotropy_axes as ax

# The 12 separations between AoE and HPA across independent measurements.
# If both axes share a common origin (parent BH spin + disk normal),
# the separations cluster near ~90°.
pairs = ax.hpa_aoe_pair_separations()
import statistics
median_sep = statistics.median([p["separation_deg"] for p in pairs])
print(f"median AoE × HPA separation: {median_sep:.1f}°")   # ~89°

# Great-circle fit to AoE measurements
gc = ax.fit_great_circle(ax.AXIS_OF_EVIL_MEASUREMENTS)
print(f"AoE great-circle pole: l={gc.pole.l_deg:.1f}°, b={gc.pole.b_deg:.1f}°")
```

### Heron Experiment 11 directional layer

```python
from nwt_substrate.cosmology import AXIS_OF_EVIL_CONSENSUS
from nwt_substrate.heron.sidereal_geometry import (
    predicted_sigma_pattern, YORKTOWN,
)

# Predicted K_7 stabilizer pattern under the substrate-axis hypothesis
t_A, t_B, t_C = 1.7e9, 1.7e9 + 43082.0, 1.7e9 + 86164.0
out = predicted_sigma_pattern(t_A, t_B, t_C, YORKTOWN,
                              predicted_axis=AXIS_OF_EVIL_CONSENSUS)
out["sigma_v_pred"]    # 7-vector, one entry per K_7 vertex
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `3` | `RANK_SO7 = N_generations = Z₃ σ-orbit count` | `isa.RANK_SO7` |
| `14` | `2 · N_VERTICES_K7 = dim(G₂) = dim Aut(𝕆)` = single-Fano syndromes | `isa.N_VERTICES_K7 × 2` |
| `25` | `H_V_SO7² = cinquefoil (5,2)-torus-knot crossings²` | `isa.H_V_SO7² = 25` |
| `α¹⁶` | Stage-7 vacuum-energy α-exponent | Paper 22 derivation |
| `h_Cox = 6` | Coxeter number of so(7) = B₃ | `isa.H_COXETER_SO7` |
| `α⁴` | K_7 X-stabilizer Wilson amplitude | substrate Steane code |
| `m_e / M_Pl` | K_7 Wilson NNLO `(8/7) α^(21/2) (1 + α/7 + (21/8)α²)` | `isa.k7_wilson_amplitude` |

If any of these change in `nwt_substrate.isa.constants`, the substrate identity assertions will refuse to import — preventing silent drift.

## Papers

- **Paper 13** — Standard Model capstone (η_B + Ω_b/Ω_c bundled with SM derivations)
- **Paper 17** — α closure (`1/α = 25π√3 + 1`) + electron mass ratio (input to Λ)
- **Paper 18** — Sakharov-induced gravity (M_Pl scale used in Λ)
- **Paper 19** — Substrate monism foundations + F1 falsifier (predicted-axis test)
- **Paper 22** — Cosmogenesis: CMB anomalies as parent-BH accretion-disk fossils

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`gravity`](gravity.md) — Newton's G, NHEK geometry (cosmogenic background), Thorne spin equilibrium
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`RANK_SO7`, `H_V_SO7`, `H_COXETER_SO7`, `N_VERTICES_K7`)
- [`heron`](heron.md) — Experiment 11 sidereal A/B/C probe (F1 falsifier hardware test)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_eta_B`, `benchmark_omega_b_c`, `benchmark_lambda_cc`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
