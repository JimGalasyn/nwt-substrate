# nwt_substrate.neutrino

> The neutrino sector of NWT (Paper 20: *The Neutrino Sector from Spin(8) Triality on K_7 / K_8*). The substrate is extended to K_8 via the `Spin(7) ⊂ Spin(8)` stabilized eighth vertex, and the Wilson-amplitude machinery of Papers 17 & 19 is reused unchanged. Predictions, **all zero-fit except where noted**: three active masses `≈ (14.84, 17.16, 52.25) meV` (`Σ m_ν ≈ 84.3 meV`), three sterile partners `≈ (61.3, 70.8, 215.7) MeV`, universal active-sterile mixing `|U_α4|² = α^(9/2) ≈ 2.4 × 10⁻¹⁰`, a tri-bi-maximal PMNS matrix with `θ_13 = arcsin(√(3α)) ≈ 8.5°`, a Dirac phase `δ_CP = -2π/3 = -120°` from the `Z_3` winding of `π_1(PSU(3))`, and `|J_CP| ≈ 0.029`. The **only** non-substrate inputs are the experimental `Δm²_21` and `Δm²_31`, used solely to extrapolate `m_2, m_3` from the substrate-predicted lightest mass `m_1`.

[← Back to index](../index.md) · Source: [`nwt_substrate/neutrino/`](../../nwt_substrate/neutrino/) · Papers: [20](https://zenodo.org/communities/nwt) (neutrino sector), [17](https://zenodo.org/records/15445103) (α + Wilson amplitude), [19](https://zenodo.org/records/15410028) (substrate monism, W3.3-J/K)

## Common questions

### Does NWT predict the neutrino masses?

The lightest active mass `m_1` is a zero-fit substrate prediction. It comes from the unified Wilson amplitude on K_8 with `(N_v, N_e) = (8, 28)`:

```
m = (DIM_S_SPIN7 / N_v) · α^(N_e/2) · (1 + α/7) · (1 + 3α²) · m_Pl
```

The same recipe with `(7, 21)` reproduces the electron mass (Paper 17), so the formula is not bespoke to neutrinos.

```python
import nwt_substrate.neutrino as nu
nu.active_mass_lightest() * 1e3    # → 14.844 meV  (Paper 20 abstract: 14.84 meV)
```

`m_2` and `m_3` are then fixed by `m_i = √(m_1² + Δm²_i1)` using the measured `Δm²_21`, `Δm²_31` (see *What is fitted vs derived?*). Result: `(14.84, 17.16, 52.25) meV`, summing to `Σ m_ν ≈ 84.3 meV` — comfortably inside the Planck CMB+BAO bound `Σ m_ν < 120 meV`.

### Does NWT predict sterile neutrinos?

Yes. Three right-handed sterile partners follow from the same Wilson recipe on a K_8 sub-graph with `(N_v, N_e) = (8, 19)`, where `19 = 28 − 9` removes the 9 cross-orbit edges of the `Z_3 ⊂ G_2` decomposition. Equivalently, a seesaw scaling `m_N_i = m_ν_i / α^(9/2)`:

```python
nu.sterile_masses()    # ≈ (6.13e7, 7.08e7, 2.16e8) eV  =  (61.3, 70.8, 215.7) MeV
```

All three land in the νMSM 50–250 MeV window targeted by SHiP, PIONEER, NA62-HL and DUNE-near. The universal active-sterile mixing is `|U_α4|² = α^(9/2) ≈ 2.4 × 10⁻¹⁰` — right at the edge of current bounds.

### Where does `δ_CP = -2π/3` come from?

From the `Z_3` winding number of a closed loop in `π_1(PSU(3)) = Z(SU(3)) ≅ Z_3` (Paper 20 §7.5). The rephasing-invariant content is `|δ_CP| = 2π/3 = 120°`; the negative sign follows the Baez Fano octonion-product convention.

```python
import math
math.degrees(nu.delta_cp_winding())    # → -120.0°
```

Feeding `δ_CP = -2π/3` and the leading-order PMNS angles into the standard formula gives the Jarlskog invariant `|J_CP| ≈ 0.029`.

### What is fitted vs derived?

This is the one place the shim is **not** zero-fit, and it is honest about it. Derived from the substrate with no fitting: `m_1`, all sterile masses, `|U_α4|²`, all three PMNS angles at leading order, `δ_CP`, and `J_CP`. The only experimental inputs are the two oscillation mass-squared splittings, used **solely** to lift `m_1 → m_2, m_3`:

```python
nu.DELTA_M_SQ_21_eV2    # 7.41e-5  eV²  (PDG 2024 / NuFIT 6.0, solar)
nu.DELTA_M_SQ_31_eV2    # 2.51e-3  eV²  (PDG 2024 / NuFIT 6.0, atmospheric, NH)
```

Everything threads back to K_7 / K_8 edge counts, Spin(7) rep dimensions and `α_NWT`. The `substrate_breakdown()` function prints the full provenance table:

```python
print(nu.substrate_breakdown())
# nwt.neutrino — Paper 20 substrate-identity breakdown
# ====================================================
# K_7 vertex count             N_VERTICES_K7 = 7
# K_8 vertex count             N_VERTICES_K8 = 8
# K_8 edge count               N_EDGES_K8    = 28
# K_8 Z_3 edge partition       (6, 3, 12, 1, 6)
#                                (6 Lorentz + 3 internal-SU(2) + 12 SM-flavor + 1 Higgs + 6 Yukawa = 28)
# K_8 active edges             K8_ACTIVE_EDGE_COUNT    = 28
# K_8 sterile edges            K8_STERILE_EDGE_COUNT   = 19
# Seesaw edge difference       K8_SEESAW_EDGE_DIFFERENCE = 9  (= 12 SM-flavor − 3 internal)
# ...
```

### Why K_8 and not K_7?

The other shims live on K_7 (7 vertices, 21 edges = `so(7)`). The neutrino sector needs the `Spin(7) ⊂ Spin(8)` triality structure, which stabilizes an eighth vertex: K_8 has 8 vertices and `28 = 8·7/2` edges, partitioned as `(6, 3, 12, 1, 6)` = Lorentz + internal-SU(2) + SM-flavor + Higgs + Yukawa. The prefactor `DIM_S_SPIN7 / N_v` becomes `8/8 = 1` for the neutrino (vs `8/7` for the electron) — a direct consequence of moving onto the 8-vertex graph.

## Prediction table

| Observable | Substrate formula | Substrate value | Experiment | Accuracy |
|---|---|---|---|---|
| `m_1` (lightest active) | K_8 Wilson `(N_v, N_e)=(8, 28)` | 14.844 meV | — (zero-fit; abstract 14.84) | 0.04 % to Paper 20 target |
| `m_2` | `√(m_1² + Δm²_21)` | 17.160 meV | uses Δm²_21 anchor | by construction |
| `m_3` | `√(m_1² + Δm²_31)` | 52.253 meV | uses Δm²_31 anchor | by construction |
| `Σ m_ν` | `m_1 + m_2 + m_3` | 84.26 meV | Planck `< 120 meV` (CMB+BAO) | inside bound |
| `m_N1` | seesaw `m_ν1 / α^(9/2)` | 61.3 MeV | νMSM window | abstract 61.3 MeV |
| `m_N2` | seesaw `m_ν2 / α^(9/2)` | 70.8 MeV | νMSM window | abstract 70.8 MeV |
| `m_N3` | seesaw `m_ν3 / α^(9/2)` | 215.7 MeV | νMSM 50–250 MeV | abstract ≈ 218.8, `±5` tol |
| `\|U_α4\|²` | `α^(9/2)` | 2.42 × 10⁻¹⁰ | edge of current bounds | abstract 2.4 × 10⁻¹⁰ |
| `θ_12` (solar) | `arctan(1/√2)` | 35.26° | NuFIT ≈ 33.4° | ~5 % (LO) |
| `θ_23` (atmospheric) | `π/4` | 45.00° | NuFIT ≈ 49.0° | ~5 % (LO) |
| `θ_13` (reactor) | `arcsin(√(3α))` | 8.51° | NuFIT ≈ 8.6° | 0.2° (≈ 1 %) |
| `δ_CP` | `-2π/3` from `π_1(PSU(3))` `Z_3` | -120.0° | T2K ≈ -110° / NuFIT ≈ 177° | discriminating; DUNE/Hyper-K |
| `\|J_CP\|` | LO angles × `sin(δ_CP)` | 0.0295 | T2K ≈ 0.027 (~0.1σ) | NuFIT global ≈ 0.003 (~1.3σ) |

`θ_12` and `θ_23` are exact tri-bi-maximal at leading order; the `α`-suppressed `U_ℓ` rotation that lands NuFIT central values (`δθ_12 ≈ -1.8°`, `δθ_23 ≈ +2.7°`) is described in Paper 20 §7.6 but **not** exposed by this shim. Each row is asserted in [`nwt_substrate/tests/test_neutrino_shim.py`](../../nwt_substrate/tests/) (e.g. `test_m1_matches_paper20` at `abs=0.1` meV, `test_sterile_mass_n1_matches_paper20` at `abs=0.5` MeV, `test_theta_13_near_paper20_target` at `abs=0.2°`, `test_jarlskog_matches_paper20` at `abs=0.002`).

## Quick start

```python
import nwt_substrate.neutrino as nu

# Active 3-mass spectrum (eV); m_1 zero-fit, m_2/m_3 via measured Δm²
nu.active_masses()            # → (0.014844, 0.017160, 0.052253) eV
nu.sum_active_masses() * 1e3  # → 84.26 meV   (Planck bound < 120 meV)

# Sterile right-handed partners (eV) and universal mixing²
nu.sterile_masses()           # → (6.13e7, 7.08e7, 2.16e8) eV = (61.3, 70.8, 215.7) MeV
nu.sterile_active_mixing_sq() # → 2.4225e-10  = α^(9/2)
nu.seesaw_ratio()             # → 2.4225e-10  (same edge-count exponent)

# PMNS leading order + CP violation
nu.pmns_angles_leading_order()  # PMNSAngles(theta_12, theta_23, theta_13) in radians
nu.delta_cp_winding()           # → -2.0943951…  = -2π/3 = -120°
nu.jarlskog_invariant()         # → -0.02954

# Unified Wilson recipe — the electron, as a cross-check
nu.wilson_mass_eV(N_v=7, N_e=21)  # → 5.1105e5 eV ≈ 0.511 MeV (Paper 17 m_e)

# Full substrate-provenance table
print(nu.substrate_breakdown())
```

## API by topic

### Mass predictions

| Function | Returns |
|---|---|
| `wilson_mass_eV(N_v, N_e, alpha=ALPHA_NWT, order="NNLO")` | Unified Wilson amplitude mass in eV; `order` ∈ {`"LO"`, `"NLO"`, `"NNLO"`} |
| `active_mass_lightest(order="NNLO")` | Lightest active mass `m_1` in eV (K_8 Wilson, `(8, 28)`) |
| `active_masses(order="NNLO")` | `(m_1, m_2, m_3)` tuple in eV (normal hierarchy) |
| `sum_active_masses(order="NNLO")` | `Σ m_ν` in eV |
| `sterile_masses(order="NNLO")` | `(m_N1, m_N2, m_N3)` sterile masses in eV |
| `seesaw_ratio()` | `m_active / m_sterile = α^(9/2)` |

### Mixing & PMNS

| Symbol | Returns |
|---|---|
| `PMNSAngles` | NamedTuple: `theta_12`, `theta_23`, `theta_13` (radians) |
| `pmns_angles_leading_order()` | `PMNSAngles` at LO (tri-bi-maximal + `θ_13`) |
| `sterile_active_mixing_sq()` | `\|U_α4\|² = α^(9/2) ≈ 2.4 × 10⁻¹⁰` |

### CP violation

| Function | Returns |
|---|---|
| `delta_cp_winding()` | Dirac phase `δ_CP = -2π/3` (radians) from `π_1(PSU(3))` `Z_3` winding |
| `jarlskog_invariant()` | `J_CP` from LO PMNS angles × `sin(δ_CP)` → ≈ -0.0295 |

### Experimental anchors (the only non-substrate inputs)

| Symbol | Value |
|---|---|
| `DELTA_M_SQ_21_eV2` | `7.41 × 10⁻⁵` eV² (solar; PDG 2024 / NuFIT 6.0) |
| `DELTA_M_SQ_31_eV2` | `2.51 × 10⁻³` eV² (atmospheric, NH; PDG 2024 / NuFIT 6.0) |

### Substrate breakdown

| Function | Returns |
|---|---|
| `substrate_breakdown()` | Multi-line string: K_8 edge decomposition + every prediction with its substrate origin |

## Worked examples

### Active and sterile mass spectra side by side

```python
import nwt_substrate.neutrino as nu

active = nu.active_masses()                      # eV
sterile = nu.sterile_masses()                    # eV
for i, (a, s) in enumerate(zip(active, sterile), start=1):
    print(f"  ν_{i}: {a*1e3:6.2f} meV   N_{i}: {s/1e6:6.1f} MeV")
# ν_1:  14.84 meV   N_1:   61.3 MeV
# ν_2:  17.16 meV   N_2:   70.8 MeV
# ν_3:  52.25 meV   N_3:  215.7 MeV
print(f"  Σ m_ν = {nu.sum_active_masses()*1e3:.2f} meV  (< 120 meV Planck)")
# Σ m_ν = 84.26 meV  (< 120 meV Planck)
```

### PMNS angles in degrees

```python
import math
import nwt_substrate.neutrino as nu

p = nu.pmns_angles_leading_order()
print(f"  θ_12 = {math.degrees(p.theta_12):.2f}°   (= arctan(1/√2),  NuFIT ≈ 33.4°)")
print(f"  θ_23 = {math.degrees(p.theta_23):.2f}°   (= π/4,           NuFIT ≈ 49.0°)")
print(f"  θ_13 = {math.degrees(p.theta_13):.2f}°    (= arcsin(√(3α)), NuFIT ≈ 8.6°)")
print(f"  δ_CP = {math.degrees(nu.delta_cp_winding()):.1f}°  |J_CP| = {abs(nu.jarlskog_invariant()):.4f}")
# θ_12 = 35.26°   (= arctan(1/√2),  NuFIT ≈ 33.4°)
# θ_23 = 45.00°   (= π/4,           NuFIT ≈ 49.0°)
# θ_13 = 8.51°    (= arcsin(√(3α)), NuFIT ≈ 8.6°)
# δ_CP = -120.0°  |J_CP| = 0.0295
```

### Tracing a number to the substrate

```python
import nwt_substrate.neutrino as nu
print(nu.substrate_breakdown())
# Active neutrino masses (eV):
#   m_1 = 1.4844e-02  (K_8 Wilson, N_v=8, N_e=28)
#   m_2 = 1.7160e-02  (extrapolated via Δm²_21 = 7.410e-05 eV²)
#   m_3 = 5.2253e-02  (extrapolated via Δm²_31 = 2.510e-03 eV²)
#   Σ m_ν = 8.4257e-02
# Sterile neutrino masses (eV):
#   m_N1 = 6.1277e+07
#   m_N2 = 7.0835e+07
#   m_N3 = 2.1570e+08
#   Seesaw ratio          α^(9/2) = 2.4225e-10
#   |U_α4|² (mixing²)             = 2.4225e-10
# PMNS angles (leading order from Spin(8) triality):
#   θ_12 = 35.2644°  = arctan(1/√2)
#   θ_23 = 45.0000°  = π/4
#   θ_13 = 8.5087°   = arcsin(√(3α))   [3 = RANK_SO7]
# Dirac CP phase from π_1(PSU(3)) Z_3 winding:
#   δ_CP   = -120.0000°  = -2π/3 (Baez Fano octonion-product convention)
#   |J_CP| = 0.0295
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `DIM_S_SPIN7 = 8` | Spinor rep of Spin(7) = `dim(octonions)`; the `8/N_v` Wilson prefactor | `isa.DIM_S_SPIN7` |
| `N_VERTICES_K8 = 8` | K_8 vertex count (`Spin(7) ⊂ Spin(8)` stabilized 8th vertex); `N_v` for ν & N | `isa.N_VERTICES_K8` |
| `N_EDGES_K8 = 28` | `8·7/2`; total K_8 edges = `(6,3,12,1,6)` partition sum | `isa.N_EDGES_K8` |
| `K8_ACTIVE_EDGE_COUNT = 28` | Active-ν edge count `N_e` (full K_8); `α^(28/2)` power | `isa.K8_ACTIVE_EDGE_COUNT` |
| `K8_STERILE_EDGE_COUNT = 19` | Sterile edge count `N_e = 28 − 9`; `α^(19/2)` power | `isa.K8_STERILE_EDGE_COUNT` |
| `K8_SEESAW_EDGE_DIFFERENCE = 9` | `28 − 19 = 12 (SM-flavor) − 3 (internal)`; seesaw exponent `α^(9/2)` | `isa.K8_SEESAW_EDGE_DIFFERENCE` |
| `RANK_SO7 = 3` | Rank of so(7) = `Z_3` σ-orbit count; the `3` in `θ_13 = arcsin(√(3α))` | `isa.RANK_SO7` |
| `NLO_VERTEX_COEFFICIENT = 1/7` | `1 / N_VERTICES_K7`; the `(1 + α/7)` NLO bracket | `isa.NLO_VERTEX_COEFFICIENT` |
| `3α²` | NNLO Wilson bracket `(1 + 3α²)` | Paper 17 → 20 NNLO |
| `α^(9/2)` | Seesaw ratio = `\|U_α4\|²` (same edge-count exponent) | `isa.ALPHA_NWT ** 4.5` |
| `3α` | `RANK_SO7 · α`; `sin²θ_13` | `isa.RANK_SO7 × ALPHA_NWT` |

These constants are asserted in `isa.constants` at import (e.g. `K8_PARTITION` sums to `N_EDGES_K8`, `K8_SEESAW_EDGE_DIFFERENCE == 28 − 19 == 12 − 3`), so any drift refuses to import.

## Papers

- **Paper 20** — *The Neutrino Sector from Spin(8) Triality on K_7 / K_8*: active + sterile mass spectra, PMNS, `δ_CP`, Jarlskog. ([NWT community on Zenodo](https://zenodo.org/communities/nwt))
- **Paper 17** — α closure (`1/α = 25π√3 + 1`) + the Wilson-amplitude mass machinery reused here. ([Zenodo](https://zenodo.org/records/15445103))
- **Paper 19** — Substrate monism foundations; W3.3-J II (active mass) and W3.3-K (tri-bi-maximal PMNS + `θ_13`). ([Zenodo](https://zenodo.org/records/15410028))

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`electroweak`](electroweak.md) — CKM Wolfenstein + `δ_CP` for the quark sector (the charged-current partner of this leptonic mixing)
- [`particles`](particles.md) — charged-lepton and hadron masses (same Wilson-amplitude lineage)
- [`cosmology`](cosmology.md) — `Σ m_ν` joins the CMB+BAO neutrino-mass bound; baryon asymmetry shares the `Z_3` σ-orbit
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`DIM_S_SPIN7`, `N_VERTICES_K8`, `N_EDGES_K8`, `K8_*`, `RANK_SO7`)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_neutrino_sector`, `benchmark_pmns_angles`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
