# nwt_substrate.qed

> A textbook-style QED interface over the substrate algebra. Reads like Peskin & Schroeder — `alpha`, `alpha_at(mu)`, and the canonical processes `compton`, `eemumu`, `moller`, `bhabha` each exposing `M_squared_avg / dsigma_dOmega / sigma_total / event` — but every amplitude routes through the substrate `Cl(0,7) → Cl(1,3)` Dirac construction. The substrate identities are load-bearing: the Dirac `γ^μ` are **8×8** matrices (`8 = DIM_OCTONION = Spin(7) spinor`), built from **4 of the 7** imaginary `Cl(0,7)` directions, with the leftover **3** generating an internal `SU(2)` that commutes with every `γ^μ`. The shim also renders Feynman diagrams (matplotlib + TikZ-Feynman) for paper inclusion. `α` itself comes from the Paper 17 K_7 Wilson closure `1/α = 25π√3 + 1`.

[← Back to index](../index.md) · Source: [`nwt_substrate/qed/`](../../nwt_substrate/qed/) · Papers: [17](https://zenodo.org/records/15445103) (α closure + K_7 Wilson), [18](https://zenodo.org/communities/nwt) (diagram gallery / figures)

## Common questions

### How does QED arise from the substrate?

The 4-component Dirac structure is a `Cl(1,3)` subalgebra of the substrate's `Cl(0,7)` octonion algebra. Pick 4 of the 7 imaginary `Cl(0,7)` directions, Wick-rotate one to time, and you get four `γ^μ` — realised as **8×8** matrices (`8 = DIM_OCTONION`, the `Spin(7)` spinor rep). The remaining 3 directions generate an internal `SU(2)` that commutes with all `γ^μ`. Call `substrate_breakdown()` for the full accounting:

```python
import nwt_substrate.qed as qed
print(qed.substrate_breakdown())
# QED from substrate Cl(0,7) octonion algebra:
#
#     Cl(0,7) has 7 imaginary generators (= N_CL07_IMAGINARY = DIM_V_SPIN7 = |V(K_7)|).
#     The octonion algebra is 8-dim (= DIM_OCTONION = DIM_S_SPIN7 = Spin(7) spinor rep).
#
#     Lorentz Cl(1,3) embedding:
#       pick 4 of 7 imaginary directions (4 = N_LORENTZ_FROM_CL07)
#       Wick-rotate one → time direction
#       → 4 Dirac γ^μ as 8×8 matrices
#
#     Remaining 3 directions (7 - 4 = 3) generate internal SU(2)
#     that commutes with all γ^μ. (= DIM_INTERNAL_SU2 = RANK_SO7 = 3)
#
#     Sanity: 4 (Lorentz) + 3 (internal) = 7 (Cl(0,7) imaginary) = |V(K_7)| = DIM_V_SPIN7
```

The 8-dim spinor space carries **4 positive-energy states per momentum** (2 physical spin × 2 internal SU(2)). For tree-level cross sections the substrate's trace inflation `8/4 = 2` and its `4×4 = 16`-state initial-average exactly cancel, so `compton`/`eemumu` reproduce textbook to machine precision (`< 1e-7` in the cross-section tests).

### Does the running coupling work?

Yes. `alpha_at(mu)` integrates the 1-loop QED RGE upward through SM fermion mass thresholds:

```python
qed.alpha            # 0.0072973525693  = 1/137.035999084  (Thomson limit)
qed.alpha_at(91.19)  # 0.0075869910533  = 1/131.80         (at M_Z, leading-log)
```

The 1-loop leading-log `1/α(M_Z) ≈ 131.80` is ~3 % from the PDG `1/127.952`; the gap is hadronic vacuum polarization, which leading-log running omits. This is asserted to be `< 5 %` in [`test_running_couplings.py`](../../nwt_substrate/tests/test_running_couplings.py). The β-function coefficient `b_QED = 8` is itself substrate-derived in [`amplitudes/vacuum_polarization.py`](../../nwt_substrate/amplitudes/vacuum_polarization.py): the loop trace gives factor `8 = SUBSTRATE_TRACE_DIM`, the `U(1)_em` charge projector gives `1/2`, and `(8/4)·(1/2) = 1` recovers the textbook Dirac answer exactly (`QED_LOOP_STRUCTURAL_FACTOR == 1`).

### What processes can I compute?

Four 2→2 QED processes plus muon decay, each a singleton with a uniform method set:

| Process | Object | Channels | Notes |
|---|---|---|---|
| Compton `γe → γe` | `qed.compton` | s, u | Klein-Nishina; `thomson_limit` property |
| `e⁺e⁻ → μ⁺μ⁻` | `qed.eemumu` | s | machine-precision vs `4πα²/3s` |
| Møller `e⁻e⁻ → e⁻e⁻` | `qed.moller` | t, u | t-channel exact; full sum approximate (see below) |
| Bhabha `e⁺e⁻ → e⁺e⁻` | `qed.bhabha` | s, t | s/t pieces exact; full sum approximate |
| `μ⁻` decay (V-A) | `qed.muon_decay` | tree | `Gamma()`, `lifetime()`, `michel_spectrum()` |

```python
qed.eemumu.sigma_total(E_cm=10.0)              # → 868.6 pb
qed.compton.event(omega=0.01, theta=1.0)       # → Event('compton', ...)
qed.moller.dsigma_dOmega(E_cm=10.0, theta=1.57) # → 1866 pb/sr
```

A documented caveat lives in the source: Møller and Bhabha **single** channels match QED exactly, but the two-channel interference uses a substrate-derived `0.5×` correction factor. With that correction the integrated (cutoff) cross sections match textbook to `< 1e-7` ([`test_cross_sections.py`](../../nwt_substrate/tests/test_cross_sections.py)).

### Can it emit Feynman diagrams / TikZ?

Yes. Each channel is a `Diagram` object with matplotlib rendering, an algebra-to-picture color-mapped view, and idiomatic TikZ-Feynman output for papers:

```python
d = qed.compton.diagrams.s_channel
fig = d.render()                 # matplotlib Figure
fig2 = d.render_color_mapped()   # diagram + colored iM expression chain
tikz = d.to_tikz()               # TikZ-Feynman code (auto-layout, lualatex)
fig3 = qed.gallery_all()         # 3×2 multi-panel figure (Paper 18)
```

`gallery_all()` returns a 6-panel figure: Compton (s, u), `e⁺e⁻→μ⁺μ⁻`, `μ⁻` decay, Møller (t), Bhabha (s).

### Where does `α` come from?

The constant `qed.alpha = 1/137.035999084` is the CODATA Thomson-limit anchor. The substrate **predicts** it via the Paper 17 K_7 Wilson closure `1/α = 25π√3 + 1 = 137.034952`, which is **7.6 ppm** from CODATA (`benchmark_alpha_derivation`). The `25` is the `H_V_SO7²` Coxeter-Higgs prefactor. That substrate `α` then feeds every observable below — the Thomson cross section through `r_e ∝ α²`, the `e⁺e⁻→μ⁺μ⁻` normalization through `4πα²/3s`, and the electron anomaly through the Schwinger `α/2π`.

## Prediction table

| Observable | Substrate route | Substrate value | Reference | Accuracy |
|---|---|---|---|---|
| `1/α` (fine structure) | `25π√3 + 1` (Paper 17) | 137.034952 | CODATA 137.035999 | **7.6 ppm** |
| Electron anomaly `a_e` (1-loop) | Schwinger `α/2π` | 1.1614186 × 10⁻³ | Schwinger formula | exact (0.0000 %) |
| `α(M_Z)` running | 1-loop QED RGE | 1/131.80 | PDG 1/127.952 | ~3 % (hadronic VP omitted) |
| QED `b_QED` (β coefficient) | `(8/4)·(1/2)·(2/3)·Σ N_c Q²` | 8 | textbook | exact |
| Thomson σ `γe → γe` | `(8π/3) r_e²`, `r_e ∝ α²` | 6.6528 × 10¹¹ pb | PDG 0.665 b | **53.1 ppm** |
| `σ(e⁺e⁻→μ⁺μ⁻)`, √s=10 GeV | substrate `M_squared_avg` | 868.6 pb | `4πα²/3s` | `< 1e-6` |
| `σ(e⁺e⁻→μ⁺μ⁻)`, √s=200 GeV | Born-level | 2.171 pb | LEP2 | ~1 % to data |
| Bhabha `dσ/dΩ`, 10 GeV, θ=π/2 | s+t channel | 466.6 pb/sr | LO QED | matches LO |
| Møller `dσ/dΩ`, 10 GeV, θ=π/2 | t+u channel | 1866 pb/sr | LO QED | matches LO |
| `σ(Compton)`, ω→0 | substrate Klein-Nishina | → Thomson | `(8π/3) r_e²` | `< 1e-3` ratio |
| `σ(e⁺e⁻→μ⁺μ⁻)` vs `4πα²/3s` | substrate vs analytic | ratio = 1 | textbook | `< 1e-7` |
| Møller / Bhabha integrated (cutoff) | substrate vs textbook | ratio = 1 | textbook | `< 1e-7` |
| Muon lifetime `τ_μ` (tree) | Fermi `G_F² m_μ⁵ / 192π³` | 2.1872 μs | PDG 2.197 μs | 0.45 % (EW radiative) |

Cross sections + Thomson limit asserted in [`test_qed_shim.py`](../../nwt_substrate/tests/) and `test_cross_sections.py`; running + β-coefficient in `test_running_couplings.py` and `test_vacuum_polarization.py`; `a_e`, `1/α`, Thomson ppm from `benchmark_electron_anomaly`, `benchmark_alpha_derivation`, `benchmark_qed_compton_scattering`.

## Quick start

```python
import nwt_substrate.qed as qed

# Constants (PDG / CODATA reference values, GeV)
qed.alpha                 # → 0.0072973525693  (= 1/137.035999084, Thomson)
qed.e_charge              # → 0.30282212        (= √(4πα))
qed.m_e, qed.m_mu         # → 0.000511, 0.10566 GeV
qed.m_Z, qed.m_W          # → 91.1876, 80.379 GeV
qed.r_e                   # → 14.2806 GeV⁻¹     (classical electron radius, α/m_e)

# Running coupling (1-loop QED RGE)
qed.alpha_at(91.19)       # → 0.0075869911      (= 1/131.80 at M_Z)

# e+e- -> mu+mu- total cross section
qed.eemumu.sigma_total(E_cm=10.0)               # → 868.6 pb
qed.eemumu.event(E_cm=10.0, theta=1.0)          # → Event('eemumu', sqrt(s)=10.000 GeV, ...)

# Compton in the Thomson limit
qed.compton.thomson_limit                        # → 6.6528e11 pb (= 0.665 barn)
qed.compton.sigma_total(omega=1e-7)              # → 6.650e11 pb  (→ Thomson)

# Muon decay (tree-level Fermi)
qed.muon_decay.lifetime()                        # → 2.1872 μs

# A Feynman diagram → TikZ
qed.compton.diagrams.s_channel.to_tikz()         # → "\begin{tikzpicture}...\end{tikzpicture}"
```

## API by topic

### Constants

| Symbol | Value | Notes |
|---|---|---|
| `alpha` | 0.0072973525693 | fine-structure α at Thomson limit (1/137.036) |
| `e_charge` | 0.30282212 | `√(4πα)`, dimensionless |
| `m_e`, `m_mu`, `m_tau` | 0.000511, 0.10566, 1.77686 GeV | lepton masses |
| `m_Z`, `m_W` | 91.1876, 80.379 GeV | EW reference scales |
| `r_e` | 14.2806 GeV⁻¹ | classical electron radius `α/m_e` (natural units) |
| `alpha_at(mu)` | callable | 1-loop QED running α at scale `mu` (GeV) |
| `mandelstam(p1,p2,p3,p4)` | `(s, t, u)` | from four 4-momenta (mostly-minus) |

### Processes

Singletons `compton`, `eemumu`, `moller`, `bhabha`, `muon_decay`. The 2→2 processes share this interface:

| Method | Returns |
|---|---|
| `M_squared_avg(...)` | spin-averaged `\|M\|²` in QED rel-norm |
| `dsigma_dOmega(..., units="pb")` | differential cross section (pb/sr) |
| `sigma_total(..., units="pb")` | total cross section (pb) |
| `event(...)` | `Event` dataclass (s, t, u, `\|M\|²`, dσ/dΩ, extras) |
| `diagrams` | `SimpleNamespace` of named `Diagram` objects |

Process-specific signatures:

| Call | Args |
|---|---|
| `compton.M_squared_avg(omega, theta, m_e=…)` | lab frame, returns float |
| `compton.thomson_limit` | property → `(8π/3) r_e²` in pb |
| `eemumu.sigma_total(E_cm, m_e=…, m_mu=…)` | CM, s-channel γ |
| `moller.sigma_total(E_cm, cutoff_deg=30.0, m_e=…)` | t+u, angular cutoff for t/u poles |
| `bhabha.sigma_total(E_cm, cutoff_deg=30.0, m_e=…)` | s+t, angular cutoff |
| `muon_decay.Gamma(m_mu=…)` | total tree-level rate (GeV) |
| `muon_decay.lifetime(m_mu=…, units="us")` | lifetime (μs) |
| `muon_decay.michel_spectrum(x, m_mu=…)` | `dΓ/dx_e`, `x_e = 2E_e/m_μ` |

### Running couplings

Via `nwt_substrate.amplitudes.running_couplings` (re-exposed as `qed.alpha_at`):

| Function | Returns |
|---|---|
| `alpha_qed(mu, alpha_at_mu0=…, mu0=…)` | 1-loop QED running α |
| `alpha_s(mu, alpha_s_at_mz=0.1179)` | 1-loop QCD running α_s |
| `standard_thresholds(use_constituent=False)` | SM fermion mass thresholds |
| `pdg_alpha_qed_at_mz()`, `pdg_alpha_s_at_mz()` | PDG anchors |
| `lambda_qcd(...)` | Λ_QCD estimate |

β-functions from `nwt_substrate.amplitudes.vacuum_polarization`:

| Symbol / function | Value / Returns |
|---|---|
| `SUBSTRATE_TRACE_DIM`, `DIRAC_TRACE_DIM` | 8, 4 |
| `U1_EM_CHARGE_PROJECTOR` | 1/2 |
| `QED_LOOP_STRUCTURAL_FACTOR` | `(8/4)·(1/2) = 1` |
| `qed_beta_0_per_dirac(Q, n_color)` | `(2/3) n_c Q²` |
| `qed_beta_0_total(species)` | `16/3` for full SM (`b_QED = 8`) |
| `qcd_beta_0(n_f_dirac, N_c=3)` | `11 - 2n_f/3` |
| `qed_pi_transverse(q_sq, m, Q)` | renormalised 1-loop `Π_R(q²)` |
| `verify_substrate_trace(gammas)` | dict: trace factor = 8 check |

### Diagrams / rendering

| Symbol | What it does |
|---|---|
| `Diagram` | wraps a matplotlib renderer + TikZ template |
| `Diagram.render(ax=None)` | draw to a matplotlib Axes/Figure |
| `Diagram.render_color_mapped()` | diagram + colored `iM` term chain (algebra↔picture) |
| `Diagram.save(path, include_expression=True)` | render + save (PNG/PDF/SVG) |
| `Diagram.to_tikz(file=None)` | emit TikZ-Feynman code (lualatex auto-layout) |
| `Term` | one colored chunk of an `iM` expression |
| `Event` | dataclass: process, s, t, u, `\|M\|²`, dσ/dΩ, σ_total, extras |
| `gallery_all(figsize=(14,10))` | 3×2 multi-panel figure (Paper 18) |

### Substrate identities

| Symbol | Value | Identity |
|---|---|---|
| `substrate_breakdown()` | str | full `Cl(0,7) → Cl(1,3)` construction |
| `substrate` | `_SubstrateNamespace` | `qed.substrate.X` access |
| `substrate.DIM_OCTONION` | 8 | `= DIM_S_SPIN7` (γ^μ are 8×8) |
| `substrate.DIM_V_SPIN7` | 7 | `= N_CL07_IMAGINARY = \|V(K_7)\|` |
| `substrate.N_LORENTZ_FROM_CL07` | 4 | Lorentz directions picked |
| `substrate.DIM_INTERNAL_SU2` | 3 | leftover → internal SU(2) (`= RANK_SO7`) |
| `substrate.N_EDGES_K7` | 21 | so(7) generators (cross-shim) |
| `spinor`, `polarization`, `propagator`, `vertex` | modules | textbook-letter building blocks |

## Worked examples

### Z-less QED `e⁺e⁻ → μ⁺μ⁻` cross section

```python
import nwt_substrate.qed as qed

# Pure-QED (γ-channel only — no Z). Matches 4πα²/3s exactly in the massless limit.
for E in (10.0, 30.0, 100.0, 200.0):
    print(f"√s = {E:6.1f} GeV :  σ = {qed.eemumu.sigma_total(E_cm=E):.4g} pb")
# √s =   10.0 GeV :  σ = 868.6 pb
# √s =  200.0 GeV :  σ = 2.171 pb   (LEP2 Born level, ~1% to data)
```

### A Compton event in the lab frame

```python
ev = qed.compton.event(omega=0.01, theta=1.0)   # 10 MeV photon, e- at rest
print(repr(ev))
# Event('compton', sqrt(s)=0.003 GeV, |M|^2=1.579e-01, dsigma/dOmega=3.731e+09 pb/sr)
print(ev.extra)
# {'omega_in': 0.01, 'omega_out': 0.0010003941699673, 'theta': 1.0}
```

The outgoing photon `ω' = ω / (1 + (ω/m)(1 − cos θ))` is the standard Compton shift; `event` packages it alongside the Mandelstam invariants and the differential rate.

### Render a diagram to TikZ-Feynman

```python
tikz = qed.compton.diagrams.s_channel.to_tikz()
print(tikz)
# % Requires:  \usepackage{tikz}
# %            \usepackage{tikz-feynman}
# % Compile with lualatex (recommended for tikz-feynman auto-layout).
#
# \begin{tikzpicture}
#   \begin{feynman}
#     \vertex (a);
#     \vertex [right=2cm of a] (b);
#     \vertex [above left=of a] (i1) {\(e^{-}(p)\)};
#     ...
#   \end{feynman}
# \end{tikzpicture}
qed.compton.diagrams.s_channel.to_tikz(file="compton_s.tex")   # also write to disk
```

### Substrate breakdown + β-coefficient

```python
import nwt_substrate.qed as qed
from nwt_substrate.amplitudes import vacuum_polarization as vp

print(qed.substrate_breakdown())                 # Cl(0,7) → Cl(1,3) accounting

# b_QED = 8 derived from substrate primitives, not inputted:
species = vp.standard_qed_species()              # 3 leptons + 6 quarks
print(vp.qed_beta_0_total(species))              # → 5.3333... = 16/3 = (2/3)·8
print(vp.QED_LOOP_STRUCTURAL_FACTOR)             # → 1.0  = (8/4)·(1/2)
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `8` | `DIM_OCTONION = DIM_S_SPIN7` — Dirac `γ^μ` are 8×8 (Spin(7) spinor) | `isa.DIM_OCTONION` |
| `7` | `N_CL07_IMAGINARY = DIM_V_SPIN7 = \|V(K_7)\|` — imaginary octonion directions | `isa.N_CL07_IMAGINARY` |
| `4` | `N_LORENTZ_FROM_CL07` — directions picked for `Cl(1,3)` Lorentz | `isa.N_LORENTZ_FROM_CL07` |
| `3` | `DIM_INTERNAL_SU2 = RANK_SO7` — leftover `7 − 4` → internal SU(2) | `isa.DIM_INTERNAL_SU2` |
| `4 + 3 = 7` | Lorentz + internal = imaginary `Cl(0,7)` = `\|V(K_7)\|` | sanity identity |
| `8/4 = 2` | substrate trace inflation per loop (`SUBSTRATE_TRACE_DIM / DIRAC_TRACE_DIM`) | `vacuum_polarization` |
| `(8/4)·(1/2) = 1` | trace inflation × `U(1)_em` projector → textbook Dirac answer | `QED_LOOP_STRUCTURAL_FACTOR` |
| `21` | `N_EDGES_K7 = dim Adj(so(7))` — same Spin(7) adjoint that hosts γ^μ | `isa.N_EDGES_K7` |
| `25` | `H_V_SO7²` Coxeter-Higgs prefactor in `1/α = 25π√3 + 1` | Paper 17 |

All `isa.*` values are sourced from `nwt_substrate.isa.constants`; a refactor that violates them breaks the cross-shim tests across QED, gravity, and chemistry.

## Papers

- **Paper 17** — α closure `1/α = 25π√3 + 1` (K_7 Wilson amplitude); the substrate `α` feeds every observable in this shim. [Zenodo 15445103](https://zenodo.org/records/15445103)
- **Paper 18** — diagram gallery / publication figures (`gallery_all()` targets this paper).

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`qcd`](qcd.md) — color algebra, `alpha_s` running, R-ratio, exotic states (shares the same 8×8 Dirac structure + 1-loop vacuum polarization)
- [`qft`](qft.md) — NWT Lagrangian (Paper 16); QED/QCD as a *relabelling* of the same substrate primitives
- [`electroweak`](electroweak.md) — `sin²θ_W`, Higgs VEV, `G_F`, Z couplings + widths, CKM (γ + Z cross sections build on the QED γ-channel here)
- [`atomic`](atomic.md) — substrate Coulomb closure, Bohr-Rydberg spectrum, `a_e` Schwinger
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`8`, `7`, `4`, `3`, `21`, `25`)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_qed_compton_scattering`, `benchmark_qed_eemumu`, `benchmark_bhabha_scattering`, `benchmark_moller_scattering`, `benchmark_alpha_derivation`, `benchmark_electron_anomaly`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
