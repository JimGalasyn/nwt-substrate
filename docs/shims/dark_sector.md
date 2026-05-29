# nwt_substrate.dark_sector

> Dark-sector phenomenology from the substrate **K_8 mass tower** `m_n = c_n · α^(N_e/2) · M_Pl`, with the rung index `N_e` set by K_8 closed-walk edge enumeration. The **98 GeV WIMP** is the `N_e = 16` rung (`m_χ = α^8 · M_Pl = 98.18 GeV`). The shim computes the Higgs-portal spin-independent direct-detection cross section `σ_SI`, the off-shell LHC production cross section, and the comparison against the **LZ-2024** limit. **Honest status:** at the default substrate coupling `g_Hχχ = α²`, the 98 GeV WIMP is **EXCLUDED** by LZ-2024 (σ_SI ≈ 6.8 × 10⁻⁴² cm², ~1.6 × 10⁶× over the limit); it becomes **ALLOWED** only if the portal coupling is `α⁴`-suppressed (σ_SI ≈ 1.9 × 10⁻⁵⁰ cm²). The mass scale is fixed; the portal coupling is the open input.

[← Back to index](../index.md) · Source: [`nwt_substrate/dark_sector/`](../../nwt_substrate/dark_sector/) · Papers: [16](https://zenodo.org/communities/nwt) (L_NWT Higgs portal), [20](https://zenodo.org/communities/nwt) (K_8 neutrino/DM tower), [19](https://zenodo.org/records/15410028) (substrate monism), [13](https://zenodo.org/records/19635239) (SM capstone)

## Common questions

### Does NWT predict a dark matter candidate?

Yes. The substrate's K_8 mass tower `m_n = c_n · α^(N_e/2) · M_Pl` spans ~20 rungs from sub-eV to `M_Pl`, with the rung index `N_e` set by K_8 closed-walk edge enumeration. The `N_e = 16` rung lands at `m_χ = α^8 · M_Pl ≈ 98.18 GeV`, squarely in the WIMP band of direct-detection experiments. The same tower also predicts warm-DM and sterile-ν species at other rungs (the `benchmark_wimp_tower` lists `N_e=22 → 38.2 eV` warm DM and `N_e=26 → 2.03 keV` sterile ν).

### Why 98 GeV?

The exponent is `N_e / 2 = 16 / 2 = 8`, so `m_χ = α^8 · M_Pl`. With the substrate `α = 1/(25π√3 + 1)` and `M_Pl = 1.22089 × 10¹⁹ GeV`, this evaluates to **98.18 GeV** — no fitting. The integer `N_e = 16` is the "Higgs-sector" rung of the K_8 enumeration; it is the same scale flagged by `benchmark_higgs_mass_vs_98gev` as a possible second Higgs-sector scalar near the LEP/ATLAS ~95 GeV hint.

### Is it ruled out by LZ-2024?

At the **default** substrate coupling `g_Hχχ = α²`, **yes — it is EXCLUDED.** Running `compare_with_lz_limit` on the default WIMP gives:

```python
{'sigma_si_cm2': 6.773623377251768e-42,
 'lz_limit_cm2': 4.2e-48,
 'ratio_pred_over_limit': 1612767.4707742305,
 'verdict': 'EXCLUDED'}
```

The predicted `σ_SI ≈ 6.8 × 10⁻⁴² cm²` is **~1.6 × 10⁶** times above the LZ-2024 limit of `4.2 × 10⁻⁴⁸ cm²` at `m_χ ≈ 100 GeV`. Because `σ_SI ∝ g_eff² ∝ g_Hχχ²`, the candidate survives only if the Higgs-portal coupling is far more suppressed: at `g_Hχχ = α⁴` the cross section drops to `σ_SI ≈ 1.9 × 10⁻⁵⁰ cm²` (ratio ≈ 0.0046), which is **ALLOWED** (below LZ-G3, possibly reachable next generation). The mass scale `α^8 · M_Pl` is fixed by the substrate; the portal coupling `g_Hχχ` is the genuinely open input the program still has to derive from K_7/K_8 portal algebra.

### How does it relate to the 95–98 GeV LHC/LEP excess?

The source flags it as a candidate match. `benchmark_higgs_mass_vs_98gev` notes that the `N_e = 16` rung at 98.18 GeV is a *second* Higgs-sector scalar distinct from the SM Higgs (the substrate gets the 125 GeV Higgs separately via `λ_H = 18α → m_h = 126.20 GeV`, ~0.9%), and tags it with "(95 GeV LEP/ATLAS hint?)". This is presented as a suggestive coincidence, not a derived signal.

### Is the portal coupling derived or assumed?

Assumed (parametrically), as the source is explicit about. The Lagrangian `L_portal = -g_Hχχ · |ψ|² · χ²` (Paper 16 L_NWT, the L_NWT complex scalar ψ) is substrate-motivated, but the dimensionless coupling `g_Hχχ` is left as a free input. The module's docstring lists candidate substrate values — `O(1)` (excluded), `α` (excluded), `α²` (default, excluded by LZ), `α⁴` (allowed) — and states plainly that fixing `g_Hχχ` from K_7/K_8 portal algebra is the substrate program's remaining task.

## Prediction table

| Observable | Substrate formula | Substrate value | Reference / limit | Status |
|---|---|---|---|---|
| WIMP mass `m_χ` | `α^8 · M_Pl` (K_8 `N_e=16`) | 98.18 GeV | WIMP band; ~95 GeV LEP/ATLAS hint | mass fixed, no fit |
| Higgs-portal coupling `g_Hχχ` | `α²` (default) | 5.325 × 10⁻⁵ | substrate input (open) | parametric |
| `g_eff = g_Hχχ · v` | `α² · 246.22 GeV` | 1.311 × 10⁻² | tree-level | — |
| `σ_SI` (at `g_Hχχ = α²`) | `g_eff² f_N² μ²/(π m_h⁴)` | 6.77 × 10⁻⁴² cm² | LZ-2024 `< 4.2 × 10⁻⁴⁸ cm²` | **EXCLUDED** (1.6 × 10⁶×) |
| `σ_SI` (at `g_Hχχ = α⁴`) | same, `α⁴` portal | 1.92 × 10⁻⁵⁰ cm² | LZ-2024 `< 4.2 × 10⁻⁴⁸ cm²` | **ALLOWED** (0.0046×) |
| LHC off-shell `σ(pp→χχ)` | off-shell Higgs, `√s=13.6 TeV` | 1.82 × 10⁻⁵ fb | — | rough estimate |
| LZ-2024 limit @ 100 GeV | external | 4.2 × 10⁻⁴⁸ cm² | arXiv:2410.17036 | reference |

Each row asserted in [`nwt_substrate/tests/test_dark_sector_wimp.py`](../../nwt_substrate/tests/test_dark_sector_wimp.py): `test_wimp_mass_is_alpha8_mpl` (90 < m_χ < 110, `n_e=16`), `test_sigma_si_excluded_at_g_alpha2` (ratio > 10⁵, verdict EXCLUDED), `test_sigma_si_allowed_at_g_alpha4` (1e-52 < σ < 1e-48, verdict ALLOWED), `test_sigma_si_scales_as_g_squared` (σ_SI ∝ g_eff²), `test_lhc_cross_section_nonnegative`, `test_predict_all_returns_dict`.

## Quick start

```python
import nwt_substrate.dark_sector as ds

# The default 98 GeV WIMP (g_Hχχ = α²)
w = ds.WIMP_98GeV()
w.mass_gev          # → 98.18060691931561   (α^8 · M_Pl, K_8 N_e=16)
w.n_e               # → 16
w.g_hxx             # → 5.325216806387437e-05   (= α²)
w.g_eff             # → 0.013111748820687148    (= g_Hχχ · v_Higgs)

# Direct-detection cross section + LZ-2024 verdict
sigma = ds.sigma_si_higgs_portal(w)         # → 6.773623377251768e-42 cm²
ds.compare_with_lz_limit(sigma)["verdict"]  # → 'EXCLUDED'  (1.6e6× over limit)

# A more suppressed portal coupling survives the LZ limit
from nwt_substrate.isa.constants import ALPHA_SUBSTRATE
w4 = ds.WIMP_98GeV(g_hxx=ALPHA_SUBSTRATE ** 4)
ds.compare_with_lz_limit(ds.sigma_si_higgs_portal(w4))["verdict"]  # → 'ALLOWED'

# LHC off-shell production estimate at √s = 13.6 TeV
ds.lhc_production_cross_section(w)          # → 1.8216221299670693e-05 fb

# Bundled summary (prints a table + returns a dict)
ds.predict_all()
```

## API by topic

### The candidate — `WIMP_98GeV`

A frozen dataclass (all fields substrate-fixed except `g_hxx`).

| Field / property | Type | Value (default) |
|---|---|---|
| `mass_gev` | `float` (field) | `ALPHA_SUBSTRATE ** 8 * M_PL_GEV` = 98.18060691931561 GeV |
| `n_e` | `int` (field) | `16` (K_8 closed-walk edge count; `N_e/2 = 8` is the α-exponent) |
| `k8_sector` | `str` (field) | `"DM-tower (carrier knot in K_8 not K_7)"` |
| `g_hxx` | `float` (field) | `ALPHA_SUBSTRATE ** 2` = 5.325216806387437e-05 (open input; override at construction) |
| `g_eff` | `property` | `g_hxx · V_HIGGS_GEV` — coefficient of `h·χ²` after EWSB |
| `reduced_mass_with_nucleon` | `property` | `m_χ m_N / (m_χ + m_N)` = 0.9301044744081993 GeV |

### Direct detection

| Function | Returns |
|---|---|
| `sigma_si_higgs_portal(wimp)` | Spin-independent χ-nucleon cross section in **cm²**: `σ_SI = g_eff² f_N² μ² / (π m_h⁴)` (tree-level) |
| `compare_with_lz_limit(sigma_si_cm2, lz_limit_at_100gev_cm2=4.2e-48)` | Dict: `sigma_si_cm2`, `lz_limit_cm2`, `ratio_pred_over_limit`, `verdict` (`'ALLOWED'` / `'EXCLUDED'`) |

### Collider

| Function | Returns |
|---|---|
| `lhc_production_cross_section(wimp, root_s_gev=13600.0)` | Rough off-shell `σ(pp → h* → χχ)` in **fb**. For `2 m_χ ≈ 196 GeV > m_h`, on-shell `h → χχ` is kinematically forbidden, so production runs through an off-shell Higgs (heavily suppressed). HONEST: a real number needs MadGraph + a substrate UFO model |

### Lagrangian

| Function | Returns |
|---|---|
| `higgs_portal_lagrangian_text()` | The L_NWT Higgs-portal Lagrangian as a multi-line string (`L_portal = -g_Hχχ·|ψ|²·χ²` → after EWSB → `-g_eff·h·χ²`) |

### Summary

| Function | Returns |
|---|---|
| `predict_all(wimp=None, show_table=True)` | Structured dict with nested `wimp`, `direct_detection` (incl. `lz_comparison`), and `lhc` blocks; optionally prints a formatted table. Defaults to the `g_Hχχ = α²` WIMP |

### Module-level constants (in `wimp_98gev`)

| Symbol | Value | Meaning |
|---|---|---|
| `V_HIGGS_GEV` | 246.22 | Higgs VEV |
| `M_HIGGS_GEV` | 125.10 | observed Higgs mass |
| `M_NUCLEON_GEV` | 0.939 | average nucleon mass |
| `F_N_HIGGS` | 0.30 | Higgs-nucleon form factor (lattice + chiral PT) |
| `GEV_INV2_TO_CM2` | (0.19732698e-13)² ≈ 3.894e-28 | GeV⁻² → cm² conversion |

## Worked examples

### Full prediction summary (default `g_Hχχ = α²`)

```python
import nwt_substrate.dark_sector as ds
ds.predict_all()
# Prints:
#   m_χ                =    98.181 GeV    (α^8 · M_Pl, K_8 N_e=16)
#   g_Hχχ              = 5.325e-05      (substrate input)
#   g_eff = g·v        = 1.311e-02
#   σ_SI (Higgs portal) = 6.774e-42 cm²
#   LZ-2024 @ 100 GeV  < 4.200e-48 cm²    EXCLUDED
#   ratio (pred/limit) = 1.613e+06
#   σ_LHC (off-shell)  = 1.822e-05 fb @ √s=13.6 TeV
```

### LZ verdict flips with the portal coupling

```python
import nwt_substrate.dark_sector as ds
from nwt_substrate.isa.constants import ALPHA_SUBSTRATE

for g in (ALPHA_SUBSTRATE ** 2, ALPHA_SUBSTRATE ** 4):
    w = ds.WIMP_98GeV(g_hxx=g)
    cmp = ds.compare_with_lz_limit(ds.sigma_si_higgs_portal(w))
    print(f"g_Hχχ = α^{8 if g < 1e-5 else 2}:  σ_SI = {cmp['sigma_si_cm2']:.2e} cm²  "
          f"ratio = {cmp['ratio_pred_over_limit']:.2e}  → {cmp['verdict']}")
# g_Hχχ = α^2:  σ_SI = 6.77e-42 cm²  ratio = 1.61e+06  → EXCLUDED
# g_Hχχ = α^8:  σ_SI = 1.92e-50 cm²  ratio = 4.57e-03  → ALLOWED
# (σ_SI ∝ g_eff², so dropping α²→α⁴ cuts the cross section by α⁴ ≈ 2.8e-9)
```

### The Higgs-portal Lagrangian text

```python
import nwt_substrate.dark_sector as ds
print(ds.higgs_portal_lagrangian_text())
#
#     L_portal  =  - g_Hχχ · |ψ|² · χ²
#
#     After EWSB:  |ψ|² → v² + 2 v h + h²
#
#     L_portal  →  - g_eff · h · χ²    with  g_eff = g_Hχχ · v_Higgs
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `N_e = 16` | K_8 closed-walk edge count for the DM rung; α-exponent `N_e/2 = 8` | K_8 mass-tower enumeration (Paper 20 sector) |
| `α^8 = α^(N_e/2)` | half-edge-count power → DM mass-tower exponent | `isa.ALPHA_SUBSTRATE = 1/(25π√3 + 1)` |
| `M_Pl` | Planck-mass anchor of the tower | `isa.M_PL_GEV = 1.22089 × 10¹⁹ GeV` |
| `K_8` (28 edges, 8 vertices) | K_7 + 1 stabilized Spin(7) ⊂ Spin(8) vertex; carrier knot lives in K_8 not K_7 | `isa.N_VERTICES_K8 = 8`, `isa.N_EDGES_K8 = 28` |
| `α²` (default `g_Hχχ`) | K_7/K_8 portal double-EM suppression (open input) | `isa.ALPHA_SUBSTRATE ** 2` |
| `v_Higgs = 246.22 GeV` | Higgs VEV (also derived in [`electroweak`](electroweak.md) at 28 ppm) | `wimp_98gev.V_HIGGS_GEV` |
| `g_eff = g_Hχχ · v` | coefficient of `h·χ²` after EWSB | `WIMP_98GeV.g_eff` |

The mass scale (`α^8 · M_Pl`) is asserted against `ALPHA_SUBSTRATE ** 8 * M_PL_GEV` in `test_dark_sector_wimp.py`; if the substrate `α` or `M_Pl` change in `nwt_substrate.isa.constants`, the test will refuse to pass, preventing silent drift.

## Papers

- **Paper 16** — L_NWT Lagrangian (`paper16_nwt_lagrangian.tex`); supplies the complex scalar ψ and the Higgs-portal coupling `L_portal = -g_Hχχ·|ψ|²·χ²`.
- **Paper 20** — K_8 neutrino/dark-sector tower; closed-walk edge enumeration that fixes the `N_e` rungs (active ν `N_e=28`, sterile ν `N_e=19`, DM `N_e=16`).
- **Paper 19** — substrate monism foundations (the K_7/K_8 substrate from which the tower descends).
- **Paper 13** — Standard Model capstone (SM-sector context the dark scalar couples into).

External limits cited in source: **LZ-2024** (arXiv:2410.17036) and **XENONnT** (arXiv:2303.14729).

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`cosmology`](cosmology.md) — `η_B`, `Ω_b/Ω_c`, `ρ_Λ` (the relic-abundance / cosmological context for the DM tower)
- [`electroweak`](electroweak.md) — Higgs VEV `v_EW = 246.21 GeV`, Higgs mass via `λ_H = 18α`, the EW sector the portal couples through
- [`particles`](particles.md) — substrate mass/width spectrum the K_8 tower extends
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`ALPHA_SUBSTRATE`, `M_PL_GEV`, `N_VERTICES_K8`, `N_EDGES_K8`)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_wimp_tower`, `benchmark_higgs_mass_vs_98gev`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
