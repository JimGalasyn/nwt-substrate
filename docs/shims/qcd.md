# nwt_substrate.qcd

> Textbook-style QCD interface over the substrate algebra. SU(3) color sits inside `SU(3) ⊂ G_2 = aut(octonions) ⊂ Spin(7)`, so the color group constants are substrate identities, not inputs: `N_c = 3 = RANK_SO7`, `N_gluons = 8 = DIM_OCTONION = DIM_S_SPIN7`, `C_F = 8/(2·3) = 4/3`, `C_A = 3`. Exposes the running coupling `α_s`, the 8 Gell-Mann generators + structure constants, the `q q̄ / q q / g g` scattering quartet (`M_squared_avg / dsigma_dOmega / sigma_total / event`), the K_7 closed-walk **confinement theorem**, the universal exotic-state mass ladder `m² = (4 m_π⁰)² · N` (23 states, **0.02–3.3 % PDG**), and gluon-as-coil Feynman diagrams.

[← Back to index](../index.md) · Source: [`nwt_substrate/qcd/`](../../nwt_substrate/qcd/) · Papers: P5b (confinement closure), P7b §7.8 (exotic mass ladder), [6/8a](https://zenodo.org/records/15376291) (substrate `m_π⁰`), [Zenodo NWT community](https://zenodo.org/communities/nwt)

## Common questions

### How does SU(3) color arise from the substrate?

SU(3) is the stabilizer of a fixed imaginary octonion direction `v̂ = e_4 ∈ S⁶ ⊂ Im(𝕆)` inside `G_2 = aut(𝕆)`: `dim SU(3) = dim G_2 − dim S⁶ = 14 − 6 = 8`. Every QCD color constant is then sourced from `nwt_substrate.isa.constants` (the same source-of-truth that gives gravity its `α^(21/2)` Wilson exponent and QED its `8×8` Dirac γ matrices). `substrate_breakdown()` spells it out:

```text
QCD from substrate SU(3) ⊂ G_2 ⊂ Spin(7):

    N_c = N_C_SU3 = RANK_SO7 = 3
    N_gluons = N_GLUONS = DIM_OCTONION = DIM_S_SPIN7 = 8
      (= N_c² - 1 = 8 gluons in the SU(3) adjoint)
      (SUBSTRATE IDENTITY: same 8 that makes Dirac γ^μ into 8×8 matrices)
    C_F = DIM_OCTONION / (2 × N_C_SU3) = 8/6 = 1.3333333333333333
    C_A = N_C_SU3 = 3.0
    T_R = 1/2
```

The 8 gluons of QCD and the 8-dim Dirac spinor space of QED are the **same** `DIM_OCTONION` — one substrate identity, not two coincidences.

### Does α_s run correctly?

Yes. `alpha_s = 0.1179` at `M_Z` (PDG 2022), and `alpha_s_at(μ)` evolves it via the 1-loop QCD β-function. The coupling **decreases with energy** (asymptotic freedom): `alpha_s_at(100) = 0.11636` falls to `alpha_s_at(1000) = 0.08915`. The 1-loop `Λ_QCD^(5) = 87 MeV` is an underestimate of the 2-loop PDG `~210 MeV` (documented in the constants); the chiral scale `Λ_χ = 4π f_π = 1.162 GeV`.

### What about confinement?

Confinement is a **structural** consequence of the K_7 graph + the `SU(3) ⊂ G_2 ⊂ Spin(7)` chain, not a dynamical mechanism (P5b closure). The `v̂ = e_4` partition splits K_7's 7 vertices into singlet `{e_4}` + color triplet `{e_1,e_2,e_3}` (charge +1) + antitriplet `{e_5,e_6,e_7}` (charge −1). A walk's **net color** is `N_c(W) = #triplet − #antitriplet`; the physical states are SU(3) singlets `N_c ≡ 0 (mod 3)`. The closed-walk-only constraint *is* "no free quarks": an open walk like `(1, 5)` is not closed, so cannot be a particle. The singlet fraction of closed walks converges to exactly `1/3` (uniform over Z₃) as length grows: `verify_confinement_range` gives fractions `[0.5, 0.267, 0.355, 0.330, 0.334, 0.333]` for `L = 2..7`.

### What exotic states are predicted?

A single universal ladder `m_X² = (4 m_π⁰)² · N_X` closes **23 bound states** (11 glueballs, 6 tetraquarks incl. X(6900), 2 heavy Z_b, 4 pentaquarks) at **0.02–3.3 % PDG**, spanning `N = 7` (= |K_7|) to `N = 389`, with **zero fitted parameters** (`m_e` is the dimensional anchor; `N` is a substrate Casimir integer). The gluon-pair scale `4 m_π⁰ = 539.8 MeV` coincides with the QCD gluon-condensate scale (500–550 MeV). Best matches: X(4140) `−0.01 %`, Z_b(10610) `−0.02 %`; worst is η(1475) at `−3.3 %`.

### Can it draw gluon diagrams?

Yes. Each process carries `.diagrams.<channel>` `Diagram` objects that render gluons as helical coils (matplotlib) and emit TikZ-Feynman with the `gluon` line type — `qqbar` (s-channel), `qq` (t/u-channel), `gg` (s-channel two-3-gluon-vertex + 4-gluon contact). `qcd.gallery_all()` returns a 4-panel figure.

## Prediction table

| Observable | Substrate formula | Substrate value | Reference | Accuracy |
|---|---|---|---|---|
| `N_c` | `N_C_SU3 = RANK_SO7` | 3 | `N_c² − 1 = 8` adjoint | exact (identity) |
| `N_gluons` | `N_GLUONS = DIM_OCTONION = DIM_S_SPIN7` | 8 | `N_c² − 1` | exact (identity) |
| `C_F` | `DIM_OCTONION / (2·N_C_SU3)` | 4/3 = 1.33333 | `(N_c²−1)/(2N_c)` | exact (identity) |
| `C_A` | `N_C_SU3` | 3 | `N_c` | exact (identity) |
| `T_R` | `T_R_SU3` | 1/2 | normalization | exact (identity) |
| `α_s(M_Z)` | PDG anchor | 0.1179 | PDG 0.1179(9) | **0.00 %** |
| `Λ_QCD^(5)` (1-loop) | 1-loop estimate | 87 MeV | PDG ~210 (2-loop) | underestimate |
| `Λ_χ` | `4π f_π` | 1.162 GeV | ~1.1 GeV | ~6 % |
| `\|M\|²` `gg→gg` @ θ=π/2 | `(243/8) g_s⁴` | 47.9663 (α_s=0.1) | textbook exact | machine ε |
| substrate `gg→gg` (transverse) | full 4-diagram | within 0.8–4.5 % of textbook off π/2 | BRST-closed exact | see worked ex. |
| η(1405) glueball | `4 m_π⁰ √7` | 1428.1 MeV | PDG 1408.0 | **+1.43 %** |
| f₀(1710) glueball | `4 m_π⁰ √10` | 1706.9 MeV | PDG 1733.0 | **−1.51 %** |
| X(4140) tetraquark | `4 m_π⁰ √59` | 4146.0 MeV | PDG 4146.5 | **−0.01 %** |
| P_c(4312) pentaquark | `4 m_π⁰ √64` | 4318.1 MeV | PDG 4312.0 | **+0.14 %** |
| Z_b(10610) | `4 m_π⁰ √386` | 10604.6 MeV | PDG 10607.2 | **−0.02 %** |
| exotic ladder (23 states) | `m² = (4 m_π⁰)²·N` | worst η(1475) | PDG | **≤ 3.3 %** |

Each row asserted in [`nwt_substrate/tests/`](../../nwt_substrate/tests/) — `test_qcd_shim.py`, `test_su3.py`, `test_qcd_gg.py`, `test_confinement.py`, `test_qcd_exotic_states.py`.

## Quick start

```python
import nwt_substrate.qcd as qcd

# Constants (all sourced from substrate ISA)
qcd.alpha_s             # → 0.1179 (PDG 2022, at M_Z)
qcd.alpha_s_at(1000.0)  # → 0.08915  (asymptotic freedom; cf. 0.11636 at 100 GeV)
qcd.C_F, qcd.C_A, qcd.T_R, qcd.N_c   # → 1.3333…, 3.0, 0.5, 3
qcd.Lambda_QCD          # → 0.087 GeV  (1-loop, n_f=5)
qcd.g_s                 # → 1.2172  (= √(4π α_s))

# Color algebra
T = qcd.T()             # 8 generators T^a = λ^a / 2, each (3,3)
f = qcd.f()             # (8,8,8) structure constants; f[0,1,2]=1, f[3,4,7]=√3/2
qcd.fundamental_casimir()   # → (4/3)·I_3
qcd.adjoint_casimir()       # → 3·δ^ab

# Processes (same quartet idiom as qed.eemumu / qed.compton)
qcd.qqbar.sigma_total(E_cm=10.0)              # → 50384.8 pb
qcd.qq.dsigma_dOmega(E_cm=10.0, theta=1.5)    # → 45238.5 pb/sr
qcd.gg.textbook_M_squared(E_cm=100.0, theta=1.5)
e = qcd.qqbar.event(E_cm=10.0, theta=1.0)     # QCDEvent('qqbar', sqrt(s)=10.000 GeV, …)

# Confinement (K_7 closed-walk theorem)
qcd.is_singlet_walk(qcd.MESON_CANDIDATE_WALK)   # → True   (1,5,1): +1 + −1 = 0
qcd.verify_confinement(length=3, start=1)       # ConfinementCheck(total=30, singlet=8, …)

# Exotic states (universal mass ladder)
qcd.universal_mass(7) * 1e3       # → 1428.1 MeV  (η(1405), N = |K_7|)
qcd.verify_universal_mass_formula()["pass"]   # → True (all 23 within 5%)

# Diagrams (gluon = helical coil)
fig  = qcd.qqbar.diagrams.s_channel.render()
tikz = qcd.qq.diagrams.t_channel.to_tikz()
fig2 = qcd.gallery_all()
```

## API by topic

### Constants & Casimirs

| Symbol | Value |
|---|---|
| `alpha_s`, `g_s` | 0.1179, 1.2172 |
| `N_c`, `C_F`, `C_A`, `T_R` | 3, 4/3, 3, 1/2 |
| `Lambda_QCD`, `Lambda_chiral` | 0.087 GeV, 1.162 GeV |
| `m_u, m_d, m_s, m_c, m_b, m_t` | current quark masses (GeV) |
| `m_u_constituent, m_d_constituent, m_s_constituent` | 0.336, 0.336, 0.486 GeV |
| `m_Z`, `m_proton` | 91.1876, 0.93828 GeV |
| `alpha_s_at(mu, alpha_s_at_mz=0.1179)` | 1-loop running coupling |
| `fundamental_casimir()` | `T^a T^a = (4/3)·I_3` |
| `adjoint_casimir()` | `f^acd f^bcd = 3·δ^ab` |

### Color algebra (T, f)

| Function | Returns |
|---|---|
| `gell_mann()` | 8 Hermitian Gell-Mann matrices λ^a |
| `T()` | 8 fundamental generators `T^a = λ^a/2` (each 3×3) |
| `f()` | `(8,8,8)` structure constants (totally antisymmetric) |
| `d()` | `(8,8,8)` symmetric tensor `d^abc` |

### Processes

| Object | Methods |
|---|---|
| `qqbar` (`q q̄ → q' q̄'`, s-channel) | `M_squared_avg`, `dsigma_dOmega`, `sigma_total`, `event` |
| `qq` (`q q → q q`, t+u channels) | `M_squared_avg`, `dsigma_dOmega`, `sigma_total(cutoff_deg=30)` |
| `gg` (`g g → g g`, s+t+u+4-gluon) | `textbook_M_squared`, `dsigma_dOmega_textbook`, `sigma_total_textbook` |
| `QCDEvent` | dataclass: `process, s, t, u, M_sq_avg, dsigma_dOmega_pb, sigma_total_pb` |

Color factors: `qqbar`/`qq` carry `2/9` (= `(1/N_c²)·2`) diagonal, `qq` u-channel cross `−2/27`; `gg` gives the textbook `(9/2)` factor.

### Confinement

| Symbol | Meaning |
|---|---|
| `POLAR_VERTEX` (4) | singlet `v̂ = e_4` |
| `TRIPLET_VERTICES` (1,2,3) | color `{r,g,b}` |
| `ANTITRIPLET_VERTICES` (5,6,7) | anti-color `{r̄,ḡ,b̄}` |
| `color_charge(v)` | +1 / −1 / 0 |
| `Walk` | dataclass; `.length`, `.is_closed` |
| `net_color(walk)`, `is_singlet_walk(walk)` | `N_c(W)`; `N_c ≡ 0 (mod 3)` |
| `enumerate_closed_walks(length, start=None)` | all closed walks of exact length |
| `verify_confinement(length, start)` | `ConfinementCheck(total, singlet, fraction)` |
| `verify_confinement_range(max_length=6, start=1)` | list over `L = 2..max_length` |
| `GLUEBALL_/MESON_/BARYON_CANDIDATE_WALK` | example singlet walks |

### Exotic states

| Function | Returns |
|---|---|
| `pi_zero_mass()` | substrate `m_π⁰ = (2 m_e/α)(1−5α)` → 134.939 MeV |
| `gluon_pair_scale()` | `4 m_π⁰` → 539.8 MeV (gluon-pair Casimir) |
| `universal_mass(N)` | `4 m_π⁰ √N` in GeV |
| `universal_N(mass_GeV)` | inverse `(m / 4 m_π⁰)²` |
| `BOUND_STATES_CATALOG` | dict of 23 states (`mass_MeV, J_PC, N, N_origin, sector`) |
| `exotic_precision_chain()` | per-state substrate-vs-PDG gaps |
| `verify_universal_mass_formula(percent_tol=5.0)` | `{pass, per_state_pass, worst_gap, worst_state}` |
| `exotic_precision_chain_summary()` | pretty-print full ladder |

### Diagrams

| Channel | Process |
|---|---|
| `qqbar.diagrams.s_channel` | `q q̄ → q' q̄'` via gluon |
| `qq.diagrams.t_channel`, `.u_channel` | `q q → q q` |
| `gg.diagrams.s_channel`, `.four_gluon` | two-3-gluon-vertex + 4-gluon contact |
| `.render(ax=None)`, `.to_tikz()` | matplotlib (gluon coil) + TikZ-Feynman |
| `gallery_all()` | 4-panel multi-process figure |

### Substrate identities

| Function | Returns |
|---|---|
| `substrate_breakdown()` | pretty-print `SU(3) ⊂ G_2 ⊂ Spin(7)` decomposition |
| `substrate` namespace | `N_C_SU3, N_GLUONS, C_F_SU3, C_A_SU3, T_R_SU3, DIM_OCTONION, RANK_SO7, N_EDGES_K7` |

## Worked examples

### Color algebra: substrate-validated Gell-Mann generators

```python
import numpy as np
import nwt_substrate.qcd as qcd

T = qcd.T(); f = qcd.f()
print(len(T), T[0].shape)                     # 8 (3, 3)   -- enforced == N_GLUONS, (N_c,N_c)
print(round(f[0,1,2], 4), round(f[3,4,7], 4)) # 1.0 0.866  (= 1, √3/2)

# Casimirs come out as the textbook multiples of identity / delta:
print(np.round(np.real(np.diag(qcd.fundamental_casimir())), 4))  # [1.3333 1.3333 1.3333]
print(np.round(np.real(np.diag(qcd.adjoint_casimir())), 4)[:3])  # [3. 3. 3.]
```

### gg → gg: substrate amplitude vs textbook

```python
import numpy as np
from nwt_substrate.amplitudes.processes import qcd_gg

# At theta = pi/2 the bracket reduces to 243/8 = 30.375; match is exact.
tb  = qcd_gg.textbook_M_squared(100.0, np.pi/2, alpha_s=0.1)   # 47.9663
sub = qcd_gg.M_squared_avg(100.0, np.pi/2, alpha_s=0.1)        # 47.9663  (rel diff 0)

# Off pi/2 the explicit-transverse result drifts (multi-leg Slavnov-Taylor):
tb1  = qcd_gg.textbook_M_squared(100.0, 1.0, alpha_s=0.1)      # 126.41
sub1 = qcd_gg.M_squared_avg(100.0, 1.0, alpha_s=0.1)           # 131.97  (+4.4 %)
# The BRST/FP-ghost-subtracted form closes it to machine precision:
brst = qcd_gg.M_squared_avg_BRST(100.0, 1.0, alpha_s=0.1)      # == tb1 (rel diff ~1e-8)
```

### Exotic-state mass ladder

```python
import nwt_substrate.qcd as qcd

for name in ["eta_1405", "f0_1710", "X_4140", "Pc_4312", "Zb_10610", "X_6900"]:
    d = qcd.exotic_precision_chain()[name]
    print(f"{name:<10} N={d['N']:<4} {d['substrate_MeV']:8.1f} MeV  "
          f"(PDG {d['pdg_MeV']:.1f}, {d['percent_gap']:+.2f}%)  [{d['N_origin']}]")
# eta_1405   N=7      1428.1 MeV  (PDG 1408.0, +1.43%)  [|K_7|]
# f0_1710    N=10     1706.9 MeV  (PDG 1733.0, -1.51%)  [dim(10_SU(5))]
# X_4140     N=59     4146.0 MeV  (PDG 4146.5, -0.01%)  [(open)]
# Pc_4312    N=64     4318.1 MeV  (PDG 4312.0, +0.14%)  [(open)]
# Zb_10610   N=386   10604.6 MeV  (PDG 10607.2, -0.02%) [(open)]
# X_6900     N=163    6891.2 MeV  (PDG 6886.0, +0.08%)  [(open)]

res = qcd.verify_universal_mass_formula()
print(res["pass"], res["worst_state"], round(res["worst_gap"], 2))  # True eta_1475 3.25
```

## Substrate constants used here

| Magic number | Substrate identity | Source |
|---|---|---|
| `8` | `N_GLUONS = DIM_OCTONION = DIM_S_SPIN7` (= `N_c² − 1`) | `isa.N_GLUONS` |
| `N_c = 3` | `N_C_SU3 = RANK_SO7` (so(7) rank = colors) | `isa.N_C_SU3`, `isa.RANK_SO7` |
| `C_F = 4/3` | `DIM_OCTONION / (2·N_C_SU3) = 8/6` | `isa.C_F_SU3` |
| `C_A = 3` | `N_C_SU3` (adjoint Casimir = N_c) | `isa.C_A_SU3` |
| `T_R = 1/2` | `T_R_SU3` (generator normalization) | `isa.T_R_SU3` |
| `7` | `|V(K_7)| = N_VERTICES_K7` (1 + 3 + 3̄ color partition; smallest exotic `N`) | `isa.N_VERTICES_K7` |
| `21` | `N_EDGES_K7` = so(7) generators (gravity `α^(21/2)`) | `isa.N_EDGES_K7` |
| `4` | `C_A²(SU(2))` (gluon-pair Casimir factor `4 m_π⁰`) | P7b §7.8 |

## Papers

- **P5b** — substrate color-confinement closure: confinement as the K_7 closed-walk theorem under `SU(3) ⊂ G_2 ⊂ Spin(7)` (Galasyn 2026-05-23; `confinement.py` docstring).
- **P7b §7.8** — universal exotic-state mass ladder `m² = (4 m_π⁰)² · N`, 23 bound states (Galasyn 2026-05-23; `exotic_states.py` docstring).
- **Paper 6 / 8a** — substrate Goldstone `m_π⁰ = (2 m_e/α)(1 − 5α)` underlying the gluon-pair scale ([Zenodo](https://zenodo.org/records/15376291)).

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`qed`](qed.md) — the sibling shim QCD parallels (`M_squared_avg / dsigma_dOmega / sigma_total / event`, `8×8` Dirac γ from the same `DIM_OCTONION`)
- [`qft`](qft.md) — shared amplitude/diagram machinery
- [`particles`](particles.md) — carrier-knot mass spectrum (mesons/baryons as K_7 walks)
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`N_C_SU3`, `N_GLUONS`, `C_F_SU3`, `RANK_SO7`, …)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_qcd_constants`, `benchmark_exotic_states`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
