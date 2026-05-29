# nwt_substrate.atomic

> The atomic-physics view of the substrate algebra: the substrate Coulomb potential `V(r) = -α ℏc / r` from the Paper 8 closure (polar-photon vertex sustained at frequency `ω`, substrate distance `r = c/ω` per light-cone tick), which yields the full Bohr-Rydberg spectrum plus precision EM for hydrogen — Rydberg (**+15 ppm PDG**), Bohr radius (**−8 ppm**), hydrogen `R_H` (**+15 ppm**), 21 cm hyperfine (**+0.056 %**), Lyman α (**−0.059 %**), fine-structure 2P (**−0.154 %**), electron `a_e` one-loop / Schwinger (**+7 ppm**), and the Lamb-shift scale. Zero free EW parameters; `m_e` is the dimensional anchor and `α` comes from Paper 17 (K_7 Wilson).

[← Back to index](../index.md) · Source: [`nwt_substrate/atomic/`](../../nwt_substrate/atomic/) · Papers: [8](https://zenodo.org/communities/nwt) (substrate Coulomb closure), [17](https://zenodo.org/records/15445103) (α + K_7 Wilson)

## Common questions

### Where does the Coulomb potential come from in NWT?

From the Paper 8 substrate closure (Galasyn 2026-05-23). The potential is

```
V(r) = -α ℏc / r       [substrate Coulomb potential]
```

derived from three substrate pieces, not postulated:
- the **P1 polar-photon vertex** carries Wilson amplitude `α` per vertex application;
- **substrate distance** is the D2 light-cone identification `r = c / ω_vertex` — each polar-photon vertex application is one "light-cone tick" at frequency `ω`, so two patches exchanging photons at rate `ω` sit at substrate-distance `r = c/ω`;
- **substrate α** is the Paper 17 K_7 Wilson value `α = 1/(25π√3 + 1)`.

The `1/r` Coulomb form and the `α` coefficient both emerge from the vertex structure, so the entire hydrogen spectrum inherits substrate-`α` precision (`+7.6 ppm` in `α` propagates to `+15 ppm` in `R_y`, which scales as `α²`).

### Does NWT predict the Rydberg constant and the Bohr radius?

Yes — both from the substrate Coulomb potential plus substrate `α`, with `m_e` as the only dimensional anchor:

```python
import nwt_substrate.atomic as at
at.rydberg()        # → 13.605901 eV   (PDG 13.605693, +15 ppm)
at.bohr_radius()    # → 52.91732 pm    (PDG 52.91772,  −8 ppm)
```

`R_y = m_e c² α²/2` and `a_0 = ℏc/(m_e c² α)`. The signs differ because `R_y ∝ α²` (over-predicts) while `a_0 ∝ 1/α` (under-predicts), both tracking the `+7.6 ppm` substrate-`α` offset.

### Does NWT predict the 21 cm line?

Yes. The 1S hyperfine splitting is `ΔE_HF = (4/3) g_p (m_e/M_p) α⁴ m_e c²`, giving `5.878 μeV` vs PDG `5.874 μeV` — **+0.056 %**, corresponding to a 21 cm wavelength near `1.421 GHz`. The proton g-factor `g_p = 5.5857` (PDG) is an input, as are `m_e` and `M_p`; `α` is substrate.

### Does NWT reproduce the electron g−2?

At one loop, yes — exactly Schwinger's 1948 result `a_e^{(1)} = α/(2π)`, evaluated at substrate `α`:

```python
at.electron_a_e_one_loop()   # → 0.0011614186   (PDG 1-loop 0.00116141, +7 ppm)
```

The `+7 ppm` gap is the substrate-`α` offset (`a_e^{(1)} ∝ α`). Higher-loop QED terms are not part of this closure; the one-loop coefficient is reproduced from substrate `α` alone.

### What about the Lamb shift?

The substrate gives the **scale** `α⁵ m_e c² ≈ 10.57 μeV` (`lamb_shift_scale`), not the full 2S–2P splitting. The PDG Lamb shift (`4.37 μeV`, `1057.845 MHz`) is recovered by multiplying this scale by the standard QED Bethe-log + cancellation factor `≈ 0.41` — a Bethe-Salpeter one-loop calculation using substrate `α`, not derived from the substrate alone at this closure level. This is flagged honestly in the docstring as "scale ✓" rather than a precision match.

### How many free parameters are there?

Zero in the atomic EM sector. `α` is fixed by Paper 17; `m_e` is the dimensional anchor (it sets the eV scale); `M_p` and `g_p` are PDG inputs that enter only the reduced-mass and hyperfine corrections. No quantity is tuned to atomic spectroscopy data.

## Prediction table

| Observable | Substrate formula | Substrate value | PDG | Accuracy |
|---|---|---|---|---|
| Rydberg `R_y` | `m_e c² α² / 2` | 13.605901 eV | 13.6056931 | **+15 ppm** |
| Bohr radius `a_0` | `ℏc / (m_e c² α)` | 52.91732 pm | 52.91772 | **−8 ppm** |
| Hydrogen `R_H` | `R_y · M_p/(M_p + m_e)` | 13.598495 eV | 13.5982886 | **+15 ppm** |
| Lyman α | `(3/4) R_H` | 10.19887 eV | 10.20486 | **−0.059 %** |
| 21 cm hyperfine | `(4/3) g_p (m_e/M_p) α⁴ m_e c²` | 5.878 μeV | 5.87433 | **+0.056 %** |
| Fine structure 2P | `m_e c² α⁴ / 32` | 45.284 μeV | 45.354 | **−0.154 %** |
| Electron `a_e` 1-loop | `α / (2π)` (Schwinger) | 0.0011614186 | 0.00116141 | **+7 ppm** |
| Lamb shift scale | `α⁵ m_e c²` | 10.57 μeV | 4.37 (× Bethe-log 0.41) | scale ✓ |

The ppm/% figures above are the signed values returned by `precision_chain()` (positive = substrate > PDG). Each row is asserted in [`nwt_substrate/tests/test_atomic_hydrogen.py`](../../nwt_substrate/tests/test_atomic_hydrogen.py).

## Quick start

```python
import nwt_substrate.atomic as at

# Bohr-Rydberg spectrum
at.rydberg()              # → 13.605901 eV   (PDG 13.605693, +15 ppm)
at.bohr_radius()          # → 52.91732 pm    (PDG 52.91772,  −8 ppm)
at.hydrogen_R_H()         # → 13.598495 eV   (reduced-mass, +15 ppm)
at.lyman_alpha()          # → 10.19887 eV    (−0.059 %)

# Precision EM
at.hyperfine_21cm()       # → 5.8776e-06 eV  (21 cm line, +0.056 %)
at.fine_structure_2P()    # → 4.5284e-05 eV  (2P splitting, −0.154 %)
at.electron_a_e_one_loop()# → 0.0011614186   (Schwinger 1-loop, +7 ppm)
at.lamb_shift_scale()     # → 1.0575e-05 eV  (α⁵ m_e c² scale)

# Physical constants (PDG)
at.M_E_EV, at.M_P_EV      # 510998.9461 eV, 938272088.16 eV
at.HBAR_C_EV_PM, at.G_P   # 197326.9804 eV·pm, 5.5856947

# Full precision chain + verification
print(at.precision_chain_summary())
assert at.verify_substrate_hydrogen()['pass']
```

## API by topic

### Physical constants (PDG)

| Symbol | Value | What it is |
|---|---|---|
| `M_E_EV` | `510998.9461` eV | Electron rest mass (dimensional anchor) |
| `M_P_EV` | `938272088.16` eV | Proton rest mass (reduced-mass + hyperfine) |
| `HBAR_C_EV_PM` | `197326.9804` eV·pm | `ℏc` (sets the length scale in `a_0`) |
| `G_P` | `5.5856947` | Proton g-factor (21 cm hyperfine input) |

### Spectrum formulas

| Function | Formula | Returns |
|---|---|---|
| `rydberg(alpha, m_e_eV)` | `m_e c² α² / 2` | Rydberg `R_y` in eV |
| `bohr_radius(alpha, m_e_eV, hbar_c_eV_pm)` | `ℏc / (m_e c² α)` | Bohr radius `a_0` in pm |
| `hydrogen_R_H(alpha, m_e_eV, M_p_eV)` | `R_y · M_p/(M_p + m_e)` | Reduced-mass `R_H` in eV |
| `lyman_alpha(alpha, m_e_eV, M_p_eV)` | `(3/4) R_H` | Lyman-α energy in eV |
| `hyperfine_21cm(alpha, m_e_eV, M_p_eV, g_p)` | `(4/3) g_p (m_e/M_p) α⁴ m_e c²` | 1S hyperfine (21 cm) in eV |
| `fine_structure_2P(alpha, m_e_eV)` | `m_e c² α⁴ / 32` | 2P₃/₂−2P₁/₂ splitting in eV |

All functions default `alpha = ALPHA_SUBSTRATE` (Paper 17 K_7 Wilson) and PDG masses, so calling them with no arguments yields the substrate prediction.

### Precision EM

| Function | Formula | Returns |
|---|---|---|
| `electron_a_e_one_loop(alpha)` | `α / (2π)` (Schwinger 1948) | Electron anomalous moment, one loop |
| `lamb_shift_scale(alpha, m_e_eV)` | `α⁵ m_e c²` | Lamb-shift energy scale in eV (no Bethe-log) |

### Precision chain + verification

| Function | Returns |
|---|---|
| `precision_chain(alpha, m_e_eV, M_p_eV, g_p, hbar_c_eV_pm)` | Dict: each of 7 observables → `{substrate, pdg, ppm}` (signed) |
| `verify_substrate_hydrogen(alpha, m_e_eV, M_p_eV, ppm_tol_ppm=100, percent_tol_percent=0.5)` | `precision_chain` dict augmented with `pass` (bool) and `per_observable_pass` (dict) |
| `precision_chain_summary(alpha, m_e_eV, M_p_eV)` | Pretty-printed substrate-vs-PDG table (string) |

The four ppm-level observables (`rydberg`, `bohr_radius`, `hydrogen_R_H`, `a_e_one_loop`) are checked against the 100 ppm tolerance; the three percent-level observables (`hyperfine_21cm`, `fine_structure_2P`, `lyman_alpha`) against the 0.5 % tolerance.

## Worked examples

### The four headline values

```python
import nwt_substrate.atomic as at
print(f"R_y    = {at.rydberg():.6f} eV        (+15 ppm)")
print(f"a_0    = {at.bohr_radius():.5f} pm        (−8 ppm)")
print(f"21 cm  = {at.hyperfine_21cm()*1e6:.4f} μeV     (+0.056 %)")
print(f"a_e^1  = {at.electron_a_e_one_loop():.10f}   (+7 ppm)")
# R_y    = 13.605901 eV        (+15 ppm)
# a_0    = 52.91732 pm        (−8 ppm)
# 21 cm  = 5.8776 μeV     (+0.056 %)
# a_e^1  = 0.0011614186   (+7 ppm)
```

### Full precision chain

```python
import nwt_substrate.atomic as at
print(at.precision_chain_summary())
# Substrate hydrogen precision chain (substrate α vs PDG):
#
#   Rydberg R_y          = 13.605901 eV    vs PDG 13.605693  →   +15.3 ppm
#   Bohr radius a_0      = 52.91732 pm   vs PDG 52.91772    →    -7.6 ppm
#   Hydrogen R_H         = 13.598495 eV    vs PDG 13.598289  →   +15.2 ppm
#   21-cm hyperfine      = 5.8776 μeV     vs PDG 5.8743   →  +0.056 %
#   2P fine structure    = 45.2840 μeV    vs PDG 45.3540  →  -0.154 %
#   Electron a_e^(1-loop) = 0.0011614186     vs PDG 0.0011614100  →    +7.4 ppm
#   Lyman-α              = 10.19887 eV     vs PDG 10.20486     →  -0.059 %
#
#   Substrate Coulomb:   V(r) = -α ℏc / r
#   Substrate distance:  r = c / ω_vertex (light-cone tick)
#   Inputs:              α (Paper 17 K_7 Wilson) + m_e (scale anchor)
#   Free params:         0 (atomic EM sector); m_e is the dimensional anchor
```

### Verification asserts all observables in tolerance

```python
import nwt_substrate.atomic as at
result = at.verify_substrate_hydrogen()
assert result['pass'] is True
assert all(result['per_observable_pass'].values())
# per_observable_pass = {
#   'rydberg': True, 'bohr_radius': True, 'hydrogen_R_H': True,
#   'hyperfine_21cm': True, 'fine_structure_2P': True,
#   'a_e_one_loop': True, 'lyman_alpha': True,
# }

# A 1 ppm tolerance fails: Rydberg sits at +15 ppm.
strict = at.verify_substrate_hydrogen(ppm_tol_ppm=1.0)
assert strict['pass'] is False
assert strict['per_observable_pass']['rydberg'] is False
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `α` | `1 / (25π√3 + 1)` — K_7 Wilson amplitude (`+7.6 ppm` CODATA) | Paper 17; imported as `electroweak.substrate_gf.ALPHA_SUBSTRATE` |
| `V(r) = -α ℏc / r` | P1 polar-photon vertex (amplitude `α`) + D2 light-cone distance | Paper 8 closure |
| `r = c / ω_vertex` | substrate distance = one light-cone tick at vertex frequency `ω` | Paper 8 / Paper 22 spacetime emergence |
| `m_e` | dimensional anchor (sets the eV energy scale); `M_E_EV = 510998.9461` | PDG |
| `ℏc` | `HBAR_C_EV_PM = 197326.9804` eV·pm (length scale in `a_0`) | PDG |
| `M_p` | `M_P_EV = 938272088.16` eV (reduced-mass + hyperfine) | PDG |
| `g_p` | `G_P = 5.5856947` (21 cm hyperfine coupling) | PDG |

`α` is shared with `electroweak.substrate_gf`, so the atomic shim and the EW shim agree on the fine-structure constant by construction (`test_uses_same_alpha_as_ew_substrate_gf`).

## Papers

- **Paper 8** — substrate Coulomb closure `V(r) = -α ℏc / r` (Galasyn 2026-05-23); hydrogen Bohr-Rydberg spectrum + precision EM.
- **Paper 17** — `α` closure (`1/α = 25π√3 + 1`), the K_7 Wilson amplitude that anchors every observable in this shim.

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`electroweak`](electroweak.md) — `α`-derived EW sector (Higgs VEV, `sin²θ_W`, `G_F`); shares `ALPHA_SUBSTRATE`
- [`particles`](particles.md) — Paper 6 mass spectrum (`m_e`, `m_p` carrier-knot masses)
- [`qed`](qed.md) — α closure + Schwinger 1-loop building block
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (single source of truth)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_atomic_hydrogen` (a₀, Lyman α, 21 cm, Lamb, Rydberg)
- [`docs/FAQ.md`](../FAQ.md) — one-liner atomic / QED Q&A
