# nwt_substrate.particles

> The particle catalog of NWT. Every particle is represented by a `Particle` object with substrate quantum numbers `(p, q, m, n_q)` from the Paper 6 carrier-knot mass formula `m / m_e = ((p² + q²) / 5) · (β/β_e) · (ln 8β / ln 8β_e) · n_q^q`. The catalog covers the **80-entry Standard Model mass spectrum** at **<1 % median PDG residual**, plus decay constants (`f_π`, `f_K`, `f_D`, `f_B`, vector mesons + `B_c`), Gell-Mann-Nishijima charge consistency, and the substrate pattern-stability ratio `ρ = m / Γ`.

[← Back to index](../index.md) · Source: [`nwt_substrate/particles/`](../../nwt_substrate/particles/) · Papers: [6](https://zenodo.org/records/15376291) (carrier-knot mass), [13](https://zenodo.org/records/19635239) (SM capstone), [21a/21b](https://zenodo.org/communities/nwt) (P7b decay synthesis)

## Common questions

### How are particles represented in NWT?

Each particle has a **carrier knot** characterized by an integer triple `(p, q)` (the torus knot indices, e.g. trefoil = T(2,3)) plus an integer `m` (substrate winding) and `n_q ∈ [0, 6]` (crossing number, one of 7 = `N_VERTICES_K7` carrier topologies):

| `n_q` | Carrier | Physics |
|---|---|---|
| 0 | unknotted | leptons (electron, muon, …) |
| 1 | extended unknot | charged leptons (lepton tower extension) |
| 2 | Hopf link | mesons (π, K, …) |
| 3 | trefoil | baryons (proton, neutron, …) |
| 4 | figure-8 | tetraquarks |
| 5 | cinquefoil | nucleons + pentaquarks (P_c, …) |
| 6 | hexa-knot | hexaquarks, dibaryons |

The integer `7 = N_CARRIER_TYPES = N_VERTICES_K7` is forced by the substrate — one carrier per K_7 vertex.

### What is the Paper 6 mass formula?

```
m / m_e = ((p² + q²) / 5) · (β / β_e) · (ln 8β / ln 8β_e) · n_q^q
```

where `5 = (p² + q²)|_electron = 2² + 1²` (the electron carrier-knot anchor), `β = (p² + q²) / (4 p q)` is the carrier-knot eccentricity, and `n_q^q` is the q-fold link enhancement. The formula reproduces the entire 24-particle compendium at **0.76 % median residual** (Paper 6) and extends to the 80-entry SM mass spectrum at **<1 % median** in Paper 13.

### What decay constants does it predict?

The unified `particles.decay_constants` module covers:
- **Light pseudoscalars** (π, K, η): Cabibbo law `f = m · √(7α)` for K, η and Fibonacci anomaly `f_π = m_π / 5^(1/4)` for the pion. **1.1-4.7 %** vs PDG.
- **Heavy pseudoscalars** (D, D_s, B, B_s): `f_X² = f_π² · R_s · N_X / m_X` with `R_s ∈ {1, √(7/4)}`. **1.1-2.6 %** vs PDG.
- **Vector mesons + B_c** (ρ, ω, K*, φ, J/ψ, Υ, D*, D_s*, B*, B_s*, B_c): `f_X* · m_X* = 7α · m_τ² · C(X*)`. **0.2-3.6 %** on the C-ratio.

See P7b §2-3, §7.5, §7.6 for the derivations.

### What is the pattern-stability ratio?

`ρ = m / Γ` quantifies how "topologically stable" each particle is — i.e., how much its substrate carrier-knot resists the gauge-current bath that drives decay. The classification: **passive** (proton, electron — `Γ` ≈ 0) versus **BPS** (`ρ` sits on a topologically protected ledge) versus **active** (heavy, `ρ` controlled by phase-space alone). The compendium claim is that **every K_7 walk is either passive or BPS** — no active outliers.

## Prediction table

| Observable | Substrate route | Accuracy |
|---|---|---|
| 24-particle compendium | Paper 6 `m / m_e` formula | 0.76 % median |
| 80-entry SM mass spectrum | Paper 13 capstone | <1 % median PDG |
| `m_p` (proton) | `(p,q,m,n_q) = (2,3,?,3)` trefoil | 0.05 % PDG |
| `f_π` decay constant | Fibonacci anomaly | −2.3 % PDG |
| `f_K` decay constant | Cabibbo scale `m · √(7α)` | +1.4 % PDG |
| `f_D`, `f_B` | `f_π² · N / m` (Paper 6 N) | 1.1 %, 1.2 % |
| `f_J/ψ` | `7α m_τ² C(X*)/m` (C=15/2) | ~6 % (combined S_V + C) |
| Gell-Mann-Nishijima `Q = T₃ + Y/2` | substrate identity | exact (integer-pair) |

Each row asserted in `nwt_substrate/tests/test_compendium.py` and the matching `nwt_substrate/tests/test_particles_decay_constants.py`.

## Quick start

```python
import nwt_substrate as nwt

# Particle by name → substrate prediction
p = nwt.particle("p")              # proton
p.mass_pred                        # → 937.24 MeV   (PDG 938.27, 0.11 %)
p.carrier                          # → "trefoil" (n_q = 3)
p.quantum_numbers                  # → (p=2, q=3, m=?, n_q=3)

# Compendium walk
for name in ["e", "mu", "tau", "p", "n", "Lambda_b"]:
    pred = nwt.particle(name).mass_pred
    print(f"{name:<10} {pred:.3f} MeV")

# Pattern-stability classification
from nwt_substrate.particles import stability_ratio_for, classify_regime
classify_regime(stability_ratio_for("p"))   # → "passive"  (proton is stable)
classify_regime(stability_ratio_for("pi+")) # → "BPS"       (long-lived charged pion)
```

## API by topic

### Particle + factory

| Function | Returns |
|---|---|
| `Particle` | Dataclass: name, `(p, q, m, n_q)`, carrier, predicted mass, PDG mass |
| `particle(name)` | Look up by canonical name ("p", "n", "mu", "Lambda_b", …) |
| `list_particles()` | All compendium keys |

### Paper 6 mass formula

| Function | Returns |
|---|---|
| `paper6_mass_ratio(p, q, m, n_q)` | `m / m_e` from the closed form |
| `paper6_mass_mev(p, q, m, n_q)` | Mass in MeV |
| `ME_MEV` | Electron mass anchor (CODATA, MeV) |

### Charge

| Function | Returns |
|---|---|
| `gell_mann_nishijima(T3, Y)` | `Q = T₃ + Y/2` — substrate-consistent electric charge |

### Decay constants ([details](#decay-constants))

| Function | Returns |
|---|---|
| `f_pi_substrate()` | `f_π = m_π⁰ / 5^(1/4)` in GeV |
| `m_tau_substrate()` | `m_τ = 25 m_e / [α(1-α)²]` in GeV |
| `vector_meson_binding_scale()` | `S_V = 7α m_τ²` (= `f_ρ m_ρ`) in GeV² |
| `cabibbo_scale_fX(m)` | `m · √(7α)` for light pseudoscalars |
| `fibonacci_anomaly_fX(m, walk_length=5)` | `m / 5^(1/4)` (pion) |
| `light_meson_fX_for("pi"|"K"|"eta")` | Name-based lookup |
| `heavy_meson_fX(N, m, strange)` | `f_π · √(R_s · N / m)` |
| `heavy_meson_fX_for("D"|"D_s"|"B"|"B_s")` | Name-based lookup |
| `vector_meson_fX(C, m)` | `7α m_τ² · C / m` |
| `vector_meson_fX_for("rho"|"K*"|"J/psi"|…|"B_c")` | Name-based lookup |
| `verify_decay_constants(tol=0.07)` | All-sector pass/fail |
| `precision_chain_summary()` | Pretty-print all 18 mesons + PDG comparison |

### Stability-ratio diagnostics

| Function | Returns |
|---|---|
| `stability_ratio(m, Gamma)` | `ρ = m / Γ` |
| `stability_ratio_for(name)` | `ρ` from compendium |
| `log10_stability_ratio_for(name)` | `log₁₀ ρ` |
| `classify_regime(rho)` | `"passive" / "BPS" / "active"` |
| `all_k7_walks_are_passive_or_BPS()` | Bool: substrate-consistency check |
| `stability_summary()` | Pretty-print full classification |

### Substrate identity

| Function | Returns |
|---|---|
| `substrate_breakdown()` | Pretty-print carrier table + Paper 6 formula derivation |
| `substrate` namespace | `N_CARRIER_TYPES`, `MAX_CROSSING_NUMBER`, `CARRIER_NAMES`, `N_VERTICES_K7` |

## Worked examples

### Compendium accuracy

```python
import nwt_substrate as nwt

# Proton — trefoil carrier
p = nwt.particle("p")
print(f"  pred = {p.mass_pred:.4f} MeV")   # 937.24
print(f"  PDG  = {p.mass_pdg:.4f} MeV")    # 938.27
print(f"  err  = {(p.mass_pred / p.mass_pdg - 1) * 100:.3f} %")  # −0.110 %

# Walk the whole compendium
errs = []
for name in nwt.particles.list_particles():
    par = nwt.particle(name)
    if par.mass_pdg:
        errs.append(abs(par.mass_pred / par.mass_pdg - 1.0) * 100)
print(f"median residual = {sorted(errs)[len(errs)//2]:.3f} %")   # < 1 %
```

### Decay constants

```python
from nwt_substrate.particles.decay_constants import (
    f_pi_substrate, light_meson_fX_for,
    heavy_meson_fX_for, vector_meson_fX_for,
)

# Substrate pion decay constant (Fibonacci-anomaly closure)
print(f"f_π   = {f_pi_substrate() * 1e3:.2f} MeV")              # 90.24
print(f"f_K   = {light_meson_fX_for('K') * 1e3:.2f} MeV")       # 111.5
print(f"f_D   = {heavy_meson_fX_for('D') * 1e3:.2f} MeV")       # 209.6
print(f"f_B_s = {heavy_meson_fX_for('B_s') * 1e3:.2f} MeV")     # 224.0
print(f"f_J/ψ = {vector_meson_fX_for('J/psi') * 1e3:.2f} MeV")  # 391.4
```

### Substrate-identity verification

```python
from nwt_substrate.particles import substrate_breakdown, all_k7_walks_are_passive_or_BPS

print(substrate_breakdown())   # carrier table + Paper 6 formula derivation

# The substrate forbids "active" walks — every K_7 carrier walk must
# be either passive or BPS-protected.
assert all_k7_walks_are_passive_or_BPS()
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `7` | `N_CARRIER_TYPES = N_VERTICES_K7` (one carrier per K_7 vertex) | `isa.N_VERTICES_K7` |
| `5` | `(p² + q²)\|_electron = 2² + 1²` (mass-formula anchor) | Paper 6 |
| `7α` | Cabibbo-Wilson amplitude (`λ² = V_us²`) | `isa.ALPHA_QED * isa.N_VERTICES_K7` |
| `√(7/4)` | `√(|K_7| / C_A²(SU(2)))` strangeness factor | recurs in v_EW NLO, Sirlin Δq |
| `5^(1/4)` | Pion Fibonacci-anomaly cinquefoil walk-length | F_5 = 5 |
| `25` | `q_cinq² = 5²` (cinquefoil signature) | also `B_s` `N_X = 25`, `m_τ = 25 m_e/[α(1-α)²]` |

## Papers

- **Paper 6** — Carrier-knot mass formula (24-particle compendium, 0.76 % median)
- **Paper 13** — Standard Model capstone (80-entry mass spectrum, <1 % median)
- **Paper 7b (P7b)** — Decay-constant synthesis (light + heavy + vector + B_c, 0.2-3.6 % C-ratio)

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`electroweak`](electroweak.md) — Higgs VEV, sin²θ_W, G_F, Z widths, CKM, form factors
- [`neutrino`](neutrino.md) — K_8 extension: 3 active + 3 sterile + PMNS
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (7, 21, 8, 3, …)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_composite_particles`, `benchmark_mass_spectrum`, `benchmark_decay_constants`, `benchmark_vector_meson_decay`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
