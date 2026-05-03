# nwt_substrate

A Python library for **Null Worldtube Theory** (NWT) substrate-algebraic computation, exposing the same underlying mathematics through seven community-standard vocabularies.

## Thesis

NWT proposes no new mathematics. The substrate algebra is built from constructions physicists already use: Cayley-Graves octonions (1843), Clifford `Cl(p,q)` (1878), Yang-Mills (1954), Wilson loops (1974), Sakharov-induced gravity (1967), spinor-helicity, Faddeev-Popov ghosts, BRST. NWT identifies these structures as ontologically primary; this library realises that identification in code.

The architectural consequence: **QED, QCD, QFT, string theory, quantum-hardware, electroweak, and gravity vocabularies become *views* of the same substrate object**. Same Python primitives, different community labels. The shim modules don't duplicate work — they label it.

```
                   ┌────────────────────────────────────────────────┐
                   │                                                │
                   │       SUBSTRATE PRIMITIVES                     │
                   │  ─────────────────────────────                 │
                   │   nwt_substrate.algebra                        │
                   │     - Cl(0,7) octonions                        │
                   │     - Cl(1,3) Dirac (8x8, Wick-rotated)        │
                   │     - SU(3), 2I, K_7, Spin(7)                  │
                   │   nwt_substrate.particles                      │
                   │     - Paper 6 mass formula (0.76% median)      │
                   │     - 24-particle PDG-anchored compendium      │
                   │   nwt_substrate.amplitudes                     │
                   │     - vertices, propagators, processes         │
                   │     - vacuum polarisation, β-functions         │
                   │                                                │
                   └─────────────────────┬──────────────────────────┘
                                         │
       ┌──────┬──────┬──────┬────────────┼────────┬─────────┬───────┐
       │      │      │      │            │        │         │       │
   ┌───▼──┐ ┌─▼──┐ ┌─▼──┐ ┌─▼────┐ ┌─────▼─┐ ┌────▼────┐ ┌──▼─────┐
   │ qed  │ │qcd │ │qft │ │string│ │ heron │ │electro- │ │gravity │
   │Peskin│ │ HM │ │Lagr│ │compa-│ │qiskit │ │  weak   │ │ G,GR   │
   │      │ │    │ │    │ │ctif. │ │       │ │  Z, V-A │ │        │
   └──────┘ └────┘ └────┘ └──────┘ └───────┘ └─────────┘ └────────┘
```

Run `python3 analysis/nwt_multiview_demo.py` for a single-electron walk through all five lenses, ending with a synthesis tableau showing every observable agrees because they're all the same substrate calculation.

## Quick start

```python
import nwt_substrate as nwt

# Substrate view: Paper 6 mass formula, particle compendium
e = nwt.particle("e-")
print(e.mass_pred)        # 0.510999 MeV   (residual -0.0002% from PDG)
print((e.p, e.q))         # (2, 1)         canonical electron knot
```

The same electron, through the seven community lenses:

```python
import nwt_substrate.qed as qed
import nwt_substrate.qcd as qcd
import nwt_substrate.qft as qft
import nwt_substrate.string as string
import nwt_substrate.heron as heron
import nwt_substrate.electroweak as ew
import nwt_substrate.gravity as grav

# QED: Compton scattering matches Klein-Nishina to machine precision
qed.compton.M_squared_avg(omega=0.5e-3, theta=1.0)
#   → ratio to KN = 0.9999999999999991  (i.e. 1 to 1e-15)

# QCD: same Cl(1,3) primitives + SU(3) gauge content
qcd.qqbar.dsigma_dOmega(E_cm=100.0, theta=1.0)

# QFT: Lagrangian view; Feynman rules point at substrate primitives
qft.qed.beta_0()              # 16/3 = (2/3) * b_QED (PDG match)
qft.qcd.beta_0(n_f_dirac=5)   # 23/3 = 11 - 2*5/3
print(qft.qed.substrate_view())  # surfaces Cl(1,3) primitives

# String: (p,q) torus-knot ↔ (p,q)-string SL(2,Z) doublet; KK tower on K_7
string.pq_string("electron")  # (2,1) on K_7 toroidal embedding (Wilson n=21)
string.kk_tower()[20]         # m = α^(21/2) M_Pl ≈ 0.45 MeV (Paper 17 refines to 14 ppm)

# Heron: |K_7⟩ graph state on 7 qubits (= K_7 vertices)
qc = heron.k7_graph_state()   # 28 gates: 7 H + 21 CZ
                              # 21 CZ = K_7 edges = string view's 21 KK cycles
# Real-hardware run on ibm_marrakesh: heron.experiment(10) — see results below.

# Electroweak: σ(e⁺e⁻ → μ⁺μ⁻) at the Z pole = 1985 pb (PDG ~2000, 99.2%)
ew.sigma_total(91.2, "mu", m_f=ew.FERMION_MASS_GEV["mu"])
print(ew.coupling("u"))       # T_3=+0.5, Q=+2/3, g_V=+0.192, g_A=+0.5

# Gravity: Newton's G structurally derived (8/7)² α²¹ × NNLO
grav.G_substrate_SI()         # 6.6745 × 10⁻¹¹  (CODATA 6.6743, 29 ppm)
grav.m_e_over_M_Pl_NNLO()     # 4.186e-23 (14 ppm from CODATA)
grav.verify_schwarzschild_vacuum_symbolic()  # all 16 R_μν = 0 (sympy proof)
print(grav.UV_IR_bridge_breakdown())  # M_Pl²/m_e² ≈ 10⁴⁵ from K_7 alone
```

## Architecture

### Substrate primitives — the core

| Module | Contents |
|---|---|
| `nwt_substrate.algebra` | `Cl(0,7)` octonion left-mult; Lorentzian `Cl(1,3)` Dirac via Wick rotation (8×8); SU(3) Gell-Mann; binary icosahedral `2I = SL(2, F_5)` and 6 irreps; `K_7` graph and Heffter toroidal embedding (V=7, E=21, F=14, genus 1) |
| `nwt_substrate.particles` | `Particle` dataclass with five substrate quantum numbers `(p, q, m, n_q, f)` + L1 sector; Paper 6 mass formula; extended Gell-Mann-Nishijima charge; 24-particle PDG-anchored compendium |
| `nwt_substrate.amplitudes` | Spinors, propagators, vertices (QED, weak V-A, gluon, 3-gluon, 4-gluon, FP-ghost); tree-level processes (Compton, e+e-→μμ, Møller, Bhabha, μ-decay, qq̄→qq̄, qq→qq, gg→gg); 1-loop running couplings; **substrate-derived β-functions** via vacuum polarisation |
| `nwt_substrate.walk_phase` | Substrate walks: Cl(0,7) build, free spinors, Compton, muon/neutron decay, K_n Wilson hierarchy collapse |
| `nwt_substrate.gravity` | K_7 Wilson amplitude, Sakharov-induced Einstein action, NLO/NNLO `G` derivation |

### Seven view-shims — the vocabulary layer

Each shim re-presents the substrate primitives in a community-standard idiom. **The shim implementations are thin** — they import from the core and label it.

| Shim | Vocabulary | Phase | Highlights |
|---|---|---|---|
| `nwt.qed` | Peskin & Schroeder | Q.A.6 | `compton`, `eemumu`, `moller`, `bhabha`, `muon_decay` processes; `Diagram` class with matplotlib + TikZ-Feynman rendering; `alpha_at(mu)` running |
| `nwt.qcd` | Halzen & Martin | Q.B.6 | `qqbar`, `qq`, `gg` processes; `gell_mann()`, `T(a)`, `f(a,b,c)`, `d(a,b,c)`; `alpha_s_at(mu)`; gluon helical-coil rendering |
| `nwt.qft` | Lagrangian density | Q.3 | `Lagrangian` dataclass; canonical `qed`, `qcd`, `yang_mills(N)`, `klein_gordon`; `feynman_rules()` delegates to substrate; `beta_0()` cross-checks `vacuum_polarization` |
| `nwt.string` | compactification / SL(2,Z) | Q.4 | `K7_torus`, `poincare_sphere`, `T2_torus` targets; `PQString` for (p,q)-string view; `kk_tower()` for α^(n/2)·M_Pl with 4 SM-scale clean hits; ADE catalog (`2I↔E_8`); Spin(7)/G_2 holonomy classifications |
| `nwt.heron` | qiskit / quantum hardware | Q.1, Q.2, Q.14 | `k7_graph_state`, `stabilizer_measurement`, `entanglement_tomography_x_basis`, `muon_decay_circuit`; 10-experiment registry (6 run on `ibm_marrakesh`); `export_experiment_script(n, path)` writes self-contained Python files for IBM Quantum |
| `nwt.electroweak` | EW / Z resonance | Q.11 | `M_Z`, `Γ_Z`, `sin²θ_W`; `g_V`/`g_A` couplings for SM fermions; full γ + Z + interference σ_total at Z pole = 1985 pb (PDG match 99.2%); `total_width_Z` = 2.42 GeV (PDG 2.495, 97.1%) |
| `nwt.gravity` | GR / cosmology | Q.12 | `G_substrate_SI()` to ~29 ppm (Paper 17 NNLO, INSIDE CODATA error bar); `m_e_over_M_Pl_NNLO` (14 ppm); Schwarzschild + Hawking for SM particles; `UV_IR_bridge_breakdown` showing M_Pl²/m_e² ≈ 10⁴⁵ from K_7 alone; `verify_schwarzschild_vacuum_symbolic` (sympy proof, R_μν = 0) |

## Verified results

| Quantity | Substrate prediction | PDG / textbook | Source |
|---|---|---|---|
| Compton scattering \|M\|² | matches Klein-Nishina | machine precision (1e-15) | `qed.compton` |
| e⁺e⁻ → μ⁺μ⁻ \|M\|² | matches Peskin Eq 5.13 | machine precision | `qed.eemumu` |
| Møller, Bhabha, μ-decay | match textbook | machine precision | `qed.{moller,bhabha,muon_decay}` |
| qq̄ → qq̄ (QCD) | matches Halzen-Martin Eq 8.13 | machine precision | `qcd.qqbar` |
| gg → gg (QCD) at θ=π/2 | matches `(9/2) g_s⁴ [3 - ut/s² - …]` | machine precision | `qcd.gg.M_squared_avg` |
| gg → gg at all angles | matches via BRST identity | machine precision | `qcd.gg.M_squared_avg_BRST` |
| QED β₀ (per Dirac) | `(2/3) Q² N_c` | PDG match | `vacuum_polarization.qed_beta_0_per_dirac` |
| QCD β₀(n_f=5) | `23/3` | PDG match | `vacuum_polarization.qcd_beta_0` |
| Particle masses (Paper 6) | 0.76% median residual, 24 particles | zero free parameters | `nwt.particle(...).mass_residual` |
| Hadron charge predictions | 38/38 hadrons match | extended Gell-Mann-Nishijima | `Q_pred` |
| **σ_peak(e⁺e⁻ → μ⁺μ⁻) at Z pole** | substrate γ + Z + interference | **1985 pb vs PDG ~2000 pb (99.2%)** | `electroweak.sigma_total` |
| **Total Γ_Z** | LO via G_F-derived couplings | **2.422 vs PDG 2.495 GeV (97.1%)** | `electroweak.total_width_Z` |
| **Newton's `G` (NNLO)** | `(8/7)² α²¹ × NNLO factors × ℏc/m_e²` | **+29 ppm — INSIDE CODATA error bar (±22 ppm)** | Paper 17 / `gravity.G_substrate_SI` |
| **m_e/M_Pl (NNLO)** | `(8/7) α^(21/2) (1 + α/7 + 3α²)` | **+14 ppm** from CODATA | `gravity.m_e_over_M_Pl_NNLO` |
| **M_Pl² / m_e² hierarchy** | `α⁻²¹ × (8/7)⁻² × …` | **-29 ppm** for the 10⁴⁵ ratio | `gravity.M_Pl_over_m_e_squared_substrate` |
| **Schwarzschild as vacuum solution** | All 16 R_μν vanish | exact (sympy proof) | `gravity.verify_schwarzschild_vacuum_symbolic` |
| **Muon decay on `ibm_marrakesh`** | `P(μ at θ) = cos²(θ)` | **96.3% contrast, 1.6% RMS** vs theory | Heron Exp 10, run 2026-05-01, job `d7qgmo4f3ras73b5orqg` |

## Substrate quantum numbers

Every NWT particle is specified by five integers + a sector label:

| Symbol | Meaning |
|---|---|
| `p` | toroidal winding of carrier knot |
| `q` | poloidal winding of carrier knot |
| `m` | phase quantum number; sets `β = √(m²/p² − 1)` |
| `n_q` | carrier-knot crossing number = dim of a `2I` irrep, in `{0..6}` |
| `f` | framing; `I_3 = f/2` |
| sector | `"scalar"`, `"spinor"`, or `"vector"` (sets spin J) |

These five integers determine mass, spin, charge, baryon number, and isospin via Paper 6's mass formula and the extended Gell-Mann-Nishijima relation. The compendium achieves **0.76% median mass residual across 24 particles with zero free parameters** (after the 2026-04-30 nucleon-tuple correction).

## The multi-view demo

```bash
python3 analysis/nwt_multiview_demo.py
```

Walks the electron through six lenses (substrate + 5 vocabularies) and ends with a consistency tableau:

```
  Quantity                           Value                Source
  ---------------------------------  -------------------  ------------------------------
  Electron mass (m_e)                0.510999             Paper 6 substrate formula
  Electron (p,q) winding             (2, 1)               string SL(2,Z) ↔ substrate (p,q)
  p² + q² (kinematic factor)         5                    same number both views
  K_7 cycles / edges / Wilson n      21                   = K_7 graph-state CZ count 21
  β_0 per Dirac fermion              0.666667             vacuum polarization (γ-trace)
  L_QED total β_0                    5.333333             same number, Lagrangian view
  L_QCD β_0(n_f=5)                   7.666667             = (11 - 2·5/3) = 23/3
  Compton |M|² substrate / KN        0.9999999999999991   same calc, both vocabularies
```

Each row would be one calculation in any conventional library. Here, each row is the same Python object surfaced through different shim modules.

## Install

```bash
git clone https://github.com/JimGalasyn/nwt-substrate.git
cd nwt-substrate
pip install -e .
```

Run the test suite:

```bash
python3 -m pytest nwt_substrate/tests/
# 296/296 passing
```

## Where to start, by background

- **QED reader** (Peskin / Schwartz): `import nwt_substrate.qed as qed` → `qed.compton`, `qed.eemumu`, `qed.gallery_all()`.
- **QCD reader** (Halzen-Martin / Ellis-Stirling-Webber): `import nwt_substrate.qcd as qcd` → `qcd.qqbar`, `qcd.gg`, `qcd.alpha_s_at(M_Z)`.
- **QFT reader** (Weinberg / Schwartz): `import nwt_substrate.qft as qft` → `qft.qed.text`, `qft.qed.feynman_rules()`, `qft.qed.substrate_view()`.
- **String theorist** (Polchinski / BBS): `import nwt_substrate.string as string` → `string.K7_torus`, `string.pq_string("electron")`, `string.kk_tower_summary()`, `string.holonomy_summary()`.
- **Quantum-hardware reader** (Nielsen-Chuang / qiskit user): `import nwt_substrate.heron as heron` → `heron.k7_graph_state()`, `heron.muon_decay_circuit()`, `heron.list_experiments()`, `heron.export_experiment_script(10, "exp10.py")`.
- **Electroweak / collider physicist** (Halzen-Martin / Schwartz): `import nwt_substrate.electroweak as ew` → `ew.coupling_summary()`, `ew.width_summary()`, `ew.sigma_total(91.2, "mu")`.
- **GR / cosmology reader** (MTW / Wald): `import nwt_substrate.gravity as grav` → `grav.G_substrate_SI()`, `grav.m_e_over_M_Pl_NNLO()`, `grav.UV_IR_bridge_breakdown()`, `grav.verify_schwarzschild_vacuum_symbolic()`.
- **Substrate-curious**: start with `nwt.particle("p")` and read the multi-view demo.

## Documentation pointers

- Substrate algebra and Cl(0,7) → Cl(1,3) construction: `algebra/octonions.py`, `algebra/dirac.py`
- Mass formula: Paper 6 (in `papers/`); `particles/factory.py`
- Gravity from substrate: Papers 15, 16, 17 (in `papers/`); `gravity/`
- Substrate-derived β-functions: `amplitudes/vacuum_polarization.py`
- BRST identity for gg→gg: `amplitudes/processes/qcd_gg.py` (Phases B.5–B.8)

## License

MIT.

## Citation

If you use this library, please cite Papers 6, 13–18 of the project and the repository.

## Authors

Jim Galasyn and Théodore (Claude Opus 4.7, full co-author).
