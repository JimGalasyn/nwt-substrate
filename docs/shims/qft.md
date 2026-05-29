# nwt_substrate.qft

> The QFT view of the substrate algebra: Lagrangian densities, field types (scalar / Dirac / gauge / Faddeev-Popov ghost), Feynman-rule extraction, and 1-loop β-functions — all as a **relabelling** of substrate-algebraic constructions that `nwt.qed` and `nwt.qcd` already compute. The substrate `Cl(1,3)` Dirac spinors *are* the 8-dimensional octonion-basis spinors; the `SU(3)` Gell-Mann generators *are* the substrate color algebra; the running couplings come from the *same* 1-loop vacuum-polarisation diagrams. Different vocabulary, same primitives. Continuum field theory for the NWT Lagrangian (**Paper 16**).

[← Back to index](../index.md) · Source: [`nwt_substrate/qft/`](../../nwt_substrate/qft/) · Papers: [16](https://zenodo.org/communities/nwt) (NWT Lagrangian, `L_NWT`)

## Common questions

### How is QFT the same thing as the substrate?

`nwt.qft` does not re-implement QFT separately from the rest of the library. Each `Lagrangian` is a QFT-vocabulary *view* whose fields and vertices point back at substrate primitives. The thesis (`qft/__init__.py` docstring): textbook quantum field theory and the substrate algebra are not separate things —

- the `Cl(1,3)` Dirac spinor is the **8-dimensional** spinor built from `Cl(0,7)` octonion left-multiplication (Wick-rotated to Lorentzian signature);
- the `SU(3)` Gell-Mann generators are the substrate's color algebra (8 gluons = `DIM_OCTONION`);
- the 3-gluon and 4-gluon vertices are encoded in standard non-Abelian Yang-Mills form;
- the running couplings come from the same 1-loop vacuum-polarisation diagrams (`amplitudes.vacuum_polarization`).

Every `Field` carries a `.substrate_primitive` string, and every `Lagrangian` exposes `.substrate_view()` that surfaces the underlying algebra.

```python
import nwt_substrate.qft as qft
print(qft.qed)                  # textbook L_QED expression
print(qft.qed.substrate_view()) # the substrate primitives behind it
```

### What β-functions does it derive?

1-loop β-function coefficients, delegated to `amplitudes.vacuum_polarization` so they are *the same numbers* used elsewhere in the library — not a parallel calculation:

```python
qft.qed.beta_0()                  # → 5.333... = 16/3   (QED, 9 charged SM fermions)
qft.qcd.beta_0(n_f_dirac=6)       # → 7.0    = 11 − 2·6/3
qft.qcd.beta_0(n_f_dirac=5)       # → 7.666... = 23/3   (below top)
qft.qcd.beta_0(n_f_dirac=0)       # → 11.0   = 11·C_A/3, pure glue
```

`β_0(QED) = (2/3)·b_QED` with `b_QED = Σ_f N_c Q²` over the 9 charged SM Dirac fermions `= 8 = DIM_OCTONION`, giving `16/3`. `β_0(QCD) = (11·C_A − 4·T_F·n_f)/3 = 11 − 2·n_f/3` for `SU(3)`. The sign tells the story: QED's `+16/3` is screening (Landau pole), QCD's `+7` (for `n_f=6`) is anti-screening (asymptotic freedom).

`Lagrangian.beta_0()` also recognises any bare non-Abelian `SU(N)` tag from `yang_mills(N)`, applying the same `11·N/3 − 2·n_f/3` structure. So `qft.yang_mills(N=3).beta_0(n_f_dirac=6)` → **`7.0`** (identical to `qft.qcd.beta_0(n_f_dirac=6)`); with no `n_f` it defaults to pure glue (`n_f = 0`): `qft.yang_mills(N=3).beta_0()` → `11.0`, `qft.yang_mills(N=2).beta_0()` → `22/3`.

### How do I compose Lagrangians?

`Lagrangian.__add__` concatenates field and interaction lists and unions the gauge-symmetry tags:

```python
L_sm = qft.qed + qft.qcd          # name "QED + QCD"
len(L_sm.fields)                  # → 18  (9 QED + 9 QCD)
len(L_sm.interactions)            # → 18  (9 QED vertices + 9 QCD vertices)
L_sm.gauge_symmetries             # → ['U(1)_em', 'SU(3)_color']
print(L_sm.fields_summary())      # pretty-print every field
```

`qft.standard_model_matter()` is exactly `qed + qcd` (Higgs / EW pieces deferred to [`electroweak`](electroweak.md)).

### What is `multiview`?

A single-object multi-lens inspector. `qft.multiview(name)` shows one substrate object through:

- **QFT view** — the `Field` dataclass (`Dirac spinor e (m=…, Q=…)`);
- **Substrate primitive** — the `.substrate_primitive` string (`Cl(1,3)` 8-dim spinor, etc.);
- **Particle compendium** — the Paper 6 `(p, q, m, n_q)` tuple + predicted mass, *when the name resolves in the compendium*.

```python
print(qft.multiview("electron"))  # QFT view + 8-dim Cl(1,3) substrate primitive
print(qft.multiview("xyzzy"))     # → "Unknown particle 'xyzzy'.  Try: [...]"
```

Accepts aliases: `e`/`electron`, `mu`/`muon`/`μ`, `tau`/`τ`, `u`/`up`, `d`/`down`, `s`, `c`, `b`, `t`, `photon`/`γ`/`A`, `gluon`/`g`/`A_g`.

### Where do the Feynman rules come from?

`Lagrangian.feynman_rules()` returns `{'propagators': {...}, 'vertices': {...}}` where the values are *references to the substrate primitives* used by the amplitude module — e.g. QED vertices point at `nwt_substrate.amplitudes.vertices.qed_vertex`, gluon self-interactions at `three_gluon_vertex_lorentz` / `four_gluon_vertex_lorentz`, propagators at `amplitudes.propagators.*`. Building a Lagrangian and building substrate amplitudes are two views of one calculation.

## Prediction table

| Fact | Substrate route | Value (verified) |
|---|---|---|
| `qft.qed.beta_0()` | `(2/3)·b_QED`, `b_QED = Σ N_c Q² = 8 = DIM_OCTONION` | **16/3 = 5.3333…** |
| `qft.qcd.beta_0(n_f_dirac=6)` | `11 − 2·6/3` | **7.0** |
| `qft.qcd.beta_0(n_f_dirac=5)` | `11 − 2·5/3` | **23/3 = 7.6667…** |
| `qft.qcd.beta_0(n_f_dirac=0)` | `11·C_A/3`, `C_A = 3` (pure glue) | **11.0** |
| `qft.yang_mills(N=3).beta_0(n_f_dirac=6)` | `11·N/3 − 2·n_f/3`, `N=3` | **7.0** |
| `qft.yang_mills(N=2).beta_0()` | `11·N/3`, `N=2` (pure glue) | **22/3 = 7.3333…** |
| QED charged-fermion count | 3 leptons + 6 quarks | **9** Dirac fermions, 9 QED vertices |
| QCD field content | 6 quarks + 1 gluon + 1 FP ghost | **8** fields, vertices: qqg ×6, 3-gluon, 4-gluon, ghost-gluon |
| `len((qed + qcd).fields)` | concatenation | **18** fields, **18** interactions |
| `SU(N)` gauge generators | `N² − 1` | 3 (N=2), 8 (N=3), 24 (N=5) |
| Dirac spinor dimension | `Cl(1,3)` octonion-basis spinor | **8** = `DIM_OCTONION` |
| `klein_gordon(charge=0)` | real scalar | `real=True`; `charge≠0` → complex |

Each row asserted in [`nwt_substrate/tests/test_qft_shim.py`](../../nwt_substrate/tests/test_qft_shim.py) (`test_qed_beta_0_matches_pdg`, `test_qcd_beta_0_matches_running_couplings_module`, `test_qed_has_9_charged_fermions`, `test_qcd_has_6_quarks_one_gluon_one_ghost`, `test_lagrangian_sum_concatenates_fields_and_vertices`, `test_qft_beta_0_consistent_with_vacuum_polarization`, …).

## Quick start

```python
import nwt_substrate.qft as qft

# Textbook Lagrangian display
print(qft.qed)
# Lagrangian QED: L_QED = − (1/4) F^{μν} F_{μν} + Σ_f ψ̄_f (i γ^μ D_μ − m_f) ψ_f, ...

# Look behind the QFT vocabulary at the substrate primitives
print(qft.qed.substrate_view())   # fields → Cl(1,3) spinors / U(1)_em generators

# Feynman rules delegate to substrate amplitude primitives
rules = qft.qed.feynman_rules()
rules['vertices']['QED vertex (e)']   # → 'nwt_substrate.amplitudes.vertices.qed_vertex'
rules['propagators']['e']             # → 'nwt_substrate.amplitudes.propagators.fermion_propagator'

# 1-loop β-functions (same numbers as vacuum_polarization)
qft.qed.beta_0()                  # → 5.333... (= 16/3)
qft.qcd.beta_0(n_f_dirac=6)       # → 7.0
qft.qcd.beta_0(n_f_dirac=5)       # → 7.666... (= 23/3)

# Compose theories
L_sm = qft.qed + qft.qcd
len(L_sm.fields), len(L_sm.interactions)   # → (18, 18)

# Generic Yang-Mills SU(N) (kinetic + 3-/4-gauge vertices)
qft.yang_mills(N=3).fields[0].n_generators # → 8 (= N² − 1)

# Free scalar
print(qft.klein_gordon("φ", mass=1.0))     # real Klein-Gordon Lagrangian

# Multi-lens inspection of one object
print(qft.multiview("electron"))
```

## API by topic

### Field types

| Symbol | Kind | Substrate primitive |
|---|---|---|
| `ScalarField(name, mass, charge, real, spin)` | complex/real scalar φ | 1-d C-valued field; prop `i/(p² − m²)` |
| `DiracSpinor(name, mass, charge, n_color, spin)` | Dirac fermion ψ | `Cl(1,3)` **8-dim** octonion-basis spinor; vertex `−i e γ^μ` |
| `GaugeField(name, gauge_group, coupling, n_generators, massless, spin)` | Yang-Mills `A^a_μ` | `g` Lie-algebra generators on fundamental rep |
| `FaddeevPopovGhost(name, gauge_group, coupling)` | FP ghost `c, c̄` | anti-commuting scalar; prop `i δ^{ab}/k²`; BRST closure |

Every field exposes a `.substrate_primitive` property and a readable `__str__`.

### Lagrangian / InteractionTerm

| Symbol | Returns |
|---|---|
| `Lagrangian(name, fields, interactions, gauge_symmetries, lorentz, text)` | a QFT Lagrangian-density view |
| `Lagrangian.feynman_rules()` | dict `{'propagators': {...}, 'vertices': {...}}` → substrate primitives |
| `Lagrangian.substrate_view()` | multi-line description of underlying substrate algebra |
| `Lagrangian.fields_summary()` | pretty-print field list |
| `Lagrangian.beta_0(n_f_dirac=None)` | 1-loop β₀ for `U(1)_em`, `SU(3)_color`, or any `SU(N)` factor (`11·N/3 − 2·n_f/3`); `None` otherwise |
| `L1 + L2` | combine — concatenate fields/interactions, union gauge tags |
| `InteractionTerm(name, fields, expression, coupling, substrate_vertex_fn)` | one vertex, pointing at a substrate vertex factor |

### Prebuilt Lagrangians

| Symbol | Theory |
|---|---|
| `qed` | `L_QED`: 9 charged Dirac fermions + photon; gauge `U(1)_em`; 9 QED vertices |
| `qcd` | `L_QCD`: 6 quarks + gluon + FP ghost; gauge `SU(3)_color`; qqg + 3-gluon + 4-gluon + ghost-gluon vertices |
| `yang_mills(N, coupling=1.0)` | pure `L_YM(SU(N))`, `N²−1` generators, 3-/4-gauge vertices |
| `klein_gordon(name="φ", mass=0.0, charge=0.0)` | free scalar (real if `charge==0`, else complex) |
| `standard_model_matter()` | `qed + qcd` (EW/Higgs deferred to `electroweak`) |

### Prebuilt particles / fields

| Symbol | Object |
|---|---|
| `electron`, `muon`, `tau` | charged-lepton `DiracSpinor`s (`n_color=1`, `Q=−1`) |
| `up_quark`, `down_quark`, `strange`, `charm`, `bottom`, `top` | quark `DiracSpinor`s (`n_color=3`) |
| `photon` | `GaugeField("A", "U(1)_em", E_CHARGE, n_generators=1)` |
| `gluon` | `GaugeField("A_g", "SU(3)_color", G_S, n_generators=8)` |
| `ghost_qcd` | `FaddeevPopovGhost("c", "SU(3)_color", G_S)` |

### multiview

| Symbol | Returns |
|---|---|
| `multiview(particle_name)` | multi-line string: QFT view + substrate primitive (+ compendium lens when resolvable) |

### Constants

| Symbol | Value |
|---|---|
| `ALPHA_QED` | `1/137.035999084` ≈ 0.00729735 |
| `E_CHARGE` | 0.30282212 (`= √(4π α)` at `q² → 0`) |
| `ALPHA_S` | 0.1179 (at `M_Z`) |
| `G_S` | 1.21720 (`= √(4π α_s)`) |

## Worked examples

### Print `L_QED` and look behind it

```python
import nwt_substrate.qft as qft

print(qft.qed.text)
# L_QED = − (1/4) F^{μν} F_{μν} + Σ_f ψ̄_f (i γ^μ D_μ − m_f) ψ_f, with D_μ = ∂_μ + i e Q_f A_μ

print(qft.qed.substrate_view())
# Lagrangian QED as substrate algebra:
#
# Fields (→ substrate primitives):
#   e   →  Cl(1,3) 8-dimensional spinor (octonion left-mult, Wick-rotated); ...
#   ...
#   A   →  U(1)_em Lie-algebra generators on fundamental rep
#
# Interaction vertices (→ substrate factors):
#   QED vertex (e)  →  nwt_substrate.amplitudes.vertices.qed_vertex
#   ...
# Gauge symmetries: U(1)_em
# Lorentz: SO(1,3)
```

### Derive the 1-loop β-functions

```python
from nwt_substrate.amplitudes import vacuum_polarization as vp

# QED: (2/3) · b_QED, b_QED = Σ N_c Q² over 9 charged SM fermions = 8
qft.qed.beta_0()                              # → 5.333... = 16/3
qft.qed.beta_0() == vp.qed_beta_0_total(vp.standard_qed_species())   # → True

# QCD: 11 − 2·n_f/3 (asymptotic freedom while n_f ≤ 16)
qft.qcd.beta_0(n_f_dirac=0)   # → 11.0     pure glue (11·C_A/3, C_A=3)
qft.qcd.beta_0(n_f_dirac=5)   # → 7.666... = 23/3
qft.qcd.beta_0(n_f_dirac=6)   # → 7.0

# Generic Yang-Mills SU(N) is recognised too (11·N/3 − 2·n_f/3)
qft.yang_mills(N=3).beta_0(n_f_dirac=6)       # → 7.0
qft.yang_mills(N=2).beta_0()                  # → 7.333... = 22/3  (pure glue)
```

### Compose the Standard-Model matter Lagrangian

```python
L_sm = qft.qed + qft.qcd
print(L_sm.name)                  # QED + QCD
print(len(L_sm.fields))           # 18
print(len(L_sm.interactions))     # 18
print(L_sm.gauge_symmetries)      # ['U(1)_em', 'SU(3)_color']
print(L_sm.fields_summary())      # one line per field

# Equivalent helper
assert len(qft.standard_model_matter().fields) == 18
```

### Multi-lens inspection

```python
print(qft.multiview("electron"))
# === Multi-view of 'electron' ===
#
# QFT view:
#   Dirac spinor e (m=0.000510998928, Q=-1.0)
#
# Substrate primitive:
#   Cl(1,3) 8-dimensional spinor (octonion left-mult, Wick-rotated); ...
```

## Substrate constants used here

| Magic number | Substrate identity | Source |
|---|---|---|
| `8` (Dirac spinor dim) | `DIM_OCTONION = DIM_S_SPIN7`; same 8 as the 8 gluons & octonion algebra | `isa.DIM_OCTONION` |
| `b_QED = 8` | `Σ_f N_c Q²` over 9 charged SM fermions `= DIM_OCTONION` | `isa.B_QED_SM = DIM_OCTONION` |
| `16/3` (β₀ QED) | `(2/3)·b_QED = (2/3)·8` | `vacuum_polarization.qed_beta_0_total` |
| `11·C_A/3` (β₀ QCD glue) | `C_A = N_c = 3 = RANK_SO7` | `isa.C_A_SU3 = N_C_SU3` |
| `4/3` (`C_F`, in `qcd_beta_0`) | `(N_c² − 1)/(2 N_c) = DIM_OCTONION/(2·N_c)` | `isa.C_F_SU3` |
| `T_F = 1/2` (quark loop index) | `Tr(T^a T^b) = T_R δ^{ab}` | `isa.T_R_SU3` |
| `8` (`SU(3)` gluons) | `N_c² − 1 = DIM_OCTONION` (gluon = adjoint rep) | `isa.N_GLUONS = DIM_OCTONION` |
| `1` (QED loop structural factor) | `(substrate trace 8 / Dirac trace 4) × U(1)_em projector (1/2) = 1` | `vacuum_polarization.QED_LOOP_STRUCTURAL_FACTOR` |

## Papers

- **Paper 16** — NWT Lagrangian (`L_NWT`, three-field / soliton sector). The `nwt.qft` shim is the continuum-field-theory view of that Lagrangian, with QED / QCD / Yang-Mills as named sub-Lagrangians. (Source references: `condensate/abelian_higgs.py` §L_2, `dark_sector/wimp_98gev.py` `L_NWT`, `benchmarks/README.md`.)

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`electroweak`](electroweak.md) — `SU(2)_L × U(1)_Y`, Higgs VEV, `sin²θ_W`, `G_F`, Z widths, CKM (the EW pieces deferred by `standard_model_matter`)
- [`particles`](particles.md) — the Paper 6 carrier-knot mass spectrum that `multiview` cross-references
- [`qed`](qed.md) — α closure + Schwinger 1-loop (the QED vertex substrate primitives)
- [`qcd`](qcd.md) — `SU(3)` color algebra, gluon vertices, running `α_s`
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`DIM_OCTONION = 8`, `N_C_SU3 = 3`, `C_F_SU3 = 4/3`, `B_QED_SM = 8`)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `L_NWT` (Paper 16) → scattering amplitudes, vacuum polarisation
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
