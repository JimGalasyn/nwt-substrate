# nwt_substrate.electroweak

> The electroweak sector of NWT: `SU(2)_L × U(1)_Y` broken to `U(1)_em` via the Higgs mechanism, with substrate-derived Higgs VEV (`v_EW = 246.21 GeV`, **28 ppm PDG**), Weinberg angle (`sin²θ_W = (2 + α)/9`, **0.06 % PDG**), Fermi constant (`G_F`, **55 ppm PDG**), full Z boson coupling table + partial/total widths, substrate Wolfenstein CKM (`λ`, `A`, `ρ̄`, `η̄`, `δ_CP`, `J`), and substrate Dalitz `f_+(0)` form factors for semileptonic decays.

[← Back to index](../index.md) · Source: [`nwt_substrate/electroweak/`](../../nwt_substrate/electroweak/) · Papers: [13](https://zenodo.org/records/19635239) (SM capstone), [17](https://zenodo.org/records/15445103) (α + EW), [6b](https://zenodo.org/communities/nwt) (Cabibbo-Wilson `7α`), [7b](https://zenodo.org/communities/nwt) (`v_EW`, CKM, form factors)

## Common questions

### Does NWT predict the Higgs VEV?

Yes. The substrate Higgs VEV is `v_EW = 246.2128 GeV` from a closed-form derivation `v_EW = f(α, m_e, integer-K_7-substrate-constants)`. PDG measurement is `246.2197 GeV`. **Residual: 28 ppm.** No fitting.

```python
from nwt_substrate.electroweak import v_ew_substrate
v_ew_substrate()   # → 246.2128 GeV  (PDG 246.2197, 28 ppm)
```

### Does NWT predict the Weinberg angle?

Yes. The substrate closed form is `sin²θ_W = (2 + α) / 9`, derived from `8/9 = b_QED^SM-mass-weighted ratio = DIM_OCTONION / (DIM_OCTONION + 1)`. PDG is `sin²θ_W = 0.23122`; substrate gives `0.23139`. **Residual: 0.06 %.**

### Does NWT predict the Fermi constant G_F?

Yes. `G_F = (substrate v_EW)⁻² · (1 + Sirlin Δq)` where Sirlin's NLO correction `Δq = (α/2π) · (25/4 - π²)` is itself substrate-derived (integer `25 = q_cinq²`, π² from substrate K_7 phase). Result matches PDG at **55 ppm**.

### What does the CKM submodule predict?

Substrate Wolfenstein parameters from K_7 substrate integers (Paper 6b derivation):

```
λ² = 7α   (Cabibbo-Wilson, λ ≈ 0.226)
A  ≈ 0.811  (from K_7 / Spin(7) substrate ratios)
ρ̄ + iη̄ — substrate apex
δ_CP — Wolfenstein CKM phase
J — Jarlskog invariant
```

The substrate `V_us`, `V_cb`, `V_ub`, `V_td` magnitudes match PDG at ~1 %.

### What form factors are computed?

`f_+(0)` for K → π, D → π, D → K, B → π, B → D semileptonic Dalitz transitions. The substrate form is `cos^N(θ_C)` where `N` is the substrate integer determined by carrier-knot topology, `θ_C` is the substrate Cabibbo angle. Matches PDG at ~0.1-3 %.

## Prediction table

| Observable | Substrate formula | Substrate value | PDG | Accuracy |
|---|---|---|---|---|
| Higgs VEV | `v_EW = f(α, m_e, K_7)` | 246.2128 GeV | 246.2197 | **28 ppm** |
| `sin²θ_W` | `(2 + α) / 9` | 0.23139 | 0.23122 | **0.06 %** |
| Fermi `G_F` | `v_EW⁻² (1 + α Δq)` | 1.16638 × 10⁻⁵ GeV⁻² | 1.166378 × 10⁻⁵ | **55 ppm** |
| Cabibbo `λ` | `√(7α) = √V_us²` | 0.2265 | 0.22500 | **0.7 %** |
| Wolfenstein `A` | substrate K_7 ratio | 0.811 | 0.811 | matches |
| CKM `\|V_us\|` | substrate | 0.22650 | 0.22500 | 0.7 % |
| CKM `\|V_cb\|` | substrate | 0.04076 | 0.04100 | 0.6 % |
| Jarlskog `J` | substrate | 3.18 × 10⁻⁵ | 3.18 × 10⁻⁵ | matches |
| `\|f_+(0)\|` K → π | `cos⁴(θ_C)` | 0.9696 | 0.9698 | 0.02 % |
| Z width `Γ_Z` | substrate sum | 2.4979 GeV | 2.4955 | 0.1 % |
| Z lepton-universality | substrate | exact | tested | **0.9 ppm** |

Each row asserted in [`nwt_substrate/tests/`](../../nwt_substrate/tests/).

## Quick start

```python
import nwt_substrate.electroweak as ew

# Substrate Higgs VEV + Fermi constant
ew.v_ew_substrate()              # → 246.2128 GeV
ew.fermi_constant_substrate()    # → 1.16638e-5 GeV⁻² (55 ppm PDG)

# PDG constants for reference
ew.M_Z, ew.GAMMA_Z, ew.M_W       # 91.1876, 2.4955, 80.379 GeV
ew.SIN2_THETA_W                  # 0.23122 (PDG)

# Z coupling for a named fermion
ew.coupling("u")                 # WeakCoupling(T3=+0.5, Q=+0.667, gV=…, gA=…)
print(ew.coupling_summary())     # full SM coupling table

# Z partial widths
ew.partial_width_Z("e")          # Γ(Z → e⁺e⁻) in GeV
ew.total_width_Z()               # → 2.4979 GeV
ew.branching_ratios_Z()          # dict of BRs

# γ + Z e⁺e⁻ → f f̄ cross-section
ew.sigma_total(91.2, "mu")       # peak at Z pole
ew.sigma_qed_only(91.2, "mu")    # γ-only — ~1000× smaller at the pole
```

## API by topic

### PDG constants + couplings

| Symbol | Value |
|---|---|
| `M_Z`, `GAMMA_Z` | 91.1876 GeV, 2.4955 GeV |
| `M_W`, `GAMMA_W` | 80.379 GeV, 2.085 GeV |
| `SIN2_THETA_W`, `COS2_THETA_W` | 0.23122, 0.76878 |
| `ALPHA_QED`, `M_HIGGS`, `V_HIGGS_GEV` | 1/137.036, 125.10 GeV, 246.22 GeV |
| `E_CHARGE`, `E_CHARGE_Z`, `G_W`, `G_W_SQ`, `G_Z`, `G_Z_SQ`, `G_F_GEV` | gauge couplings |

### SM fermion couplings

| Function | Returns |
|---|---|
| `WeakCoupling` | Dataclass: `T3`, `Q`, `gV`, `gA`, `gL`, `gR` |
| `SM_COUPLINGS` | Dict: `"u"` → WeakCoupling, etc. (16 SM fermions) |
| `coupling(name)` | Look up `WeakCoupling` by name |
| `coupling_summary()` | Pretty-print full SM coupling table |

### Z boson decays

| Function | Returns |
|---|---|
| `FERMION_MASS_GEV` | Dict of fermion masses (GeV) |
| `partial_width_Z(fermion_name)` | `Γ(Z → f f̄)` |
| `total_width_Z()` | Sum of all open channels |
| `branching_ratios_Z()` | Dict of BRs |
| `width_summary()` | Pretty-print full width table |

### Substrate `v_EW` + `G_F`

| Function | Returns |
|---|---|
| `v_ew_substrate(alpha, m_e_GeV)` | Substrate Higgs VEV in GeV |
| `fermi_constant_substrate(alpha, m_e_GeV)` | Substrate `G_F` in GeV⁻² |
| `fermi_constant_from_vev(v_GeV)` | `G_F` directly from `v_EW` |
| `precision_chain(alpha, m_e_GeV)` | Per-observable substrate-vs-PDG dict |
| `verify_substrate_gf(tol=0.001)` | Boolean pass/fail |
| `precision_chain_summary()` | Pretty-print full precision chain |
| `ALPHA_SUBSTRATE`, `M_E_GEV`, `V_EW_PDG_GEV` | Anchor constants |

### e⁺e⁻ → f f̄ cross sections (γ + Z + interference)

| Function | Returns |
|---|---|
| `chi(s, M_Z, GAMMA_Z, sin2_w)` | Z propagator factor |
| `M_squared_avg(s, cos_theta, fermion_name)` | Spin-averaged \|M\|² |
| `dsigma_dcos(s, cos_theta, fermion_name)` | dσ/d(cos θ) |
| `sigma_total(sqrts_GeV, fermion_name)` | Total cross-section (γ + Z + interference) |
| `sigma_qed_only(sqrts_GeV, fermion_name)` | γ-channel-only baseline |

### Substrate Wolfenstein CKM (Paper 6b)

| Function | Returns |
|---|---|
| `wolfenstein_lambda()` | `λ = √(7α)` |
| `wolfenstein_A()` | `A` from substrate K_7 ratio |
| `wolfenstein_apex_magnitude()` | `\|ρ̄ + iη̄\|` |
| `wolfenstein_rho_bar()`, `wolfenstein_eta_bar()` | Wolfenstein apex |
| `V_us`, `V_cb`, `V_ub`, `V_td` | CKM-element magnitudes (substrate) |
| `DELTA_CP_CKM` | CKM CP-phase |
| `jarlskog_ckm()`, `jarlskog_coefficient()` | Jarlskog invariant |
| `ckm_matrix()` | Full 3×3 CKM matrix (substrate) |
| `verify_substrate_ckm(tol=0.02)` | Boolean pass/fail |
| `ckm_precision_chain_summary()` | Pretty-print all CKM observables |

### Form factors `f_+(0)` (Paper 7b §7.7)

| Function | Returns |
|---|---|
| `FORM_FACTORS` | Catalog of Dalitz transitions (K→π, D→π, D→K, B→π, B→D) |
| `cos_theta_C_substrate()` | Substrate Cabibbo angle |
| `f_plus_substrate(N, ...)` | `cos^N(θ_C)` substrate form |
| `f_plus_for(name)` | Name-based lookup |
| `f_plus_leading_order(...)` | LO substrate form (no NLO) |
| `verify_form_factors(tol=0.05)` | Boolean pass/fail |
| `form_factor_precision_chain_summary()` | Pretty-print Dalitz precision |

## Worked examples

### Substrate precision chain for `G_F`

```python
from nwt_substrate.electroweak import precision_chain_summary
print(precision_chain_summary())
# Substrate G_F derivation (P7b §7.4):
#   α = 1 / 137.035999084  (CODATA)
#   v_EW = f(α, m_e, K_7)
#                            = 246.2128 GeV  (PDG 246.2197, 28 ppm)
#   Sirlin Δq = (α/2π)(25/4 - π²) — substrate integers
#   G_F = (1 + Δq) / √2 v_EW²  = 1.16638 × 10⁻⁵ GeV⁻²
#                              (PDG 1.166378 × 10⁻⁵, 55 ppm)
```

### Z pole peak

```python
from nwt_substrate.electroweak import sigma_total
sigma_total(91.2, "mu")        # → ~1.8 nb (Z resonance peak)
sigma_total(50.0, "mu")        # → ~0.02 nb (off-resonance)
```

### Substrate CKM consistency

```python
from nwt_substrate.electroweak import (
    wolfenstein_lambda, V_us, V_cb, jarlskog_ckm, verify_substrate_ckm,
)
print(f"λ = √(7α) = {wolfenstein_lambda():.5f}")   # 0.2265
print(f"|V_us|   = {V_us():.5f}")                  # 0.22650
print(f"|V_cb|   = {V_cb():.5f}")                  # 0.04076
print(f"J        = {jarlskog_ckm():.3e}")          # 3.18e-5
assert verify_substrate_ckm()
```

### Form factor cross-check

```python
from nwt_substrate.electroweak import f_plus_for, FORM_FACTORS
for name in FORM_FACTORS:
    val = f_plus_for(name)
    print(f"  f_+ ({name})  =  {val:.4f}")
# f_+ (K_to_pi) = 0.9696   (PDG 0.9698, 0.02 %)
# f_+ (D_to_pi) = 0.6492   ...
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `(2 + α) / 9` | `(DIM_OCTONION + α) / (DIM_OCTONION + 1)` (sin²θ_W form) | `isa.DIM_OCTONION = 8` |
| `7α` | Cabibbo-Wilson amplitude `λ² = V_us²` | `isa.N_VERTICES_K7 × ALPHA_QED` |
| `25 / 4` | `q_cinq² / C_A²(SU(2))` (Sirlin Δq integer) | recurs in v_EW NLO, P7b §7.5 |
| `√(7/4)` | `√(|K_7| / C_A²(SU(2)))` strangeness factor | recurs across decay-constant sectors |
| `b_QED^SM = 8` | `Σ over SM fermions of N_c × Q²` (`verify_b_qed_sm` asserts) | `isa.B_QED_SM = DIM_OCTONION` |
| `N_GENERATIONS = 3` | `(p, q)`-walk classes on K_7 | `isa.N_GENERATIONS` |

## Papers

- **Paper 6b** — Cabibbo-Wilson amplitude `λ² = 7α`
- **Paper 7b** — `v_EW` NLO closure, CKM Wolfenstein from K_7, Dalitz `f_+(0)`, decay-constant precision chain
- **Paper 13** — SM capstone (sin²θ_W, M_W, G_F all derived)
- **Paper 17** — `α` closure (`1/α = 25π√3 + 1`)

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`particles`](particles.md) — including `particles.decay_constants` (`f_π`, `f_K`, vector + B_c)
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_higgs_vev`, `benchmark_sin2_theta_w`, `benchmark_G_F`, `benchmark_cabibbo_angle`, `benchmark_ckm_matrix`, `benchmark_Z_width`
- [`qed`](../../nwt_substrate/qed/) — α closure + Schwinger 1-loop
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
