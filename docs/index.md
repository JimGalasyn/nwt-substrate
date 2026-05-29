# nwt-substrate documentation

> The reference documentation for `nwt-substrate`, the Python library that implements the substrate algebra of Null Worldtube Theory (NWT). Every numerical claim in the [NWT paper series](https://zenodo.org/communities/nwt) is reproducible by calling a function documented here.

## What to read first

| You are… | Start here |
|---|---|
| A human exploring the library | [`README.md`](../README.md) |
| An AI search engine or agent | [`llms.txt`](../llms.txt) (short) or [`llms-full.txt`](../llms-full.txt) (long) |
| Looking for a one-liner answer | [`docs/FAQ.md`](FAQ.md) |
| An AI coding agent contributing to the repo | [`AGENTS.md`](../AGENTS.md) |
| A human contributor | [`CONTRIBUTING.md`](../CONTRIBUTING.md) |
| Looking for version history / what changed | [`CHANGELOG.md`](../CHANGELOG.md) (structured) · [`docs/releases/`](releases/) (narrative) |
| Looking for a specific shim's API | This page → [Shim reference](#shim-reference) below |
| Looking for the 38 forward predictions | [`nwt_substrate/benchmarks/README.md`](../nwt_substrate/benchmarks/README.md) |
| Looking for the substrate constants | [`nwt_substrate/isa/README.md`](../nwt_substrate/isa/README.md) |

## Install

```bash
pip install git+https://github.com/JimGalasyn/nwt-substrate.git
```

Python 3.10, 3.11, 3.12. Optional extras: `[heron]` (qiskit-runtime for IBM Heron), `[torch]`, `[test]` (pytest + astropy).

## Three-line proof

```python
import nwt_substrate as nwt
nwt.particle("p").mass_pred                          # → 937.24 MeV (proton, Paper 6)
nwt.gravity.G_substrate_SI()                         # → 6.674228e-11, −11 ppm CODATA (Paper 17)
nwt.isa.k7_wilson_amplitude(1/137.036, order="NNLO") # → 4.185e-23 = m_e / M_Pl, −5.5 ppm
```

Three independent substrate predictions — particle mass, gravitational coupling, the underlying K_7 Wilson amplitude — all matching CODATA / PDG at zero free parameters beyond `m_e`, `M_Pl`, `c`, `ℏ`.

## Library architecture

```
nwt_substrate/
├── isa/                  ← SUBSTRATE: constants + algebra + observables
├── benchmarks/           ← 38 substrate-vs-experiment benchmarks (~100 ms)
├── particles/            ← Paper 6 mass formula, particle catalog, decay constants
├── gravity/              ← Sakharov-induced G, BH thermodynamics, NHEK, cosmogenesis
├── electroweak/          ← Higgs VEV, sin²θ_W, G_F, Z couplings + widths, CKM, form factors
├── cosmology/            ← η_B, Ω_b/Ω_c, Λ, CMB anisotropy axes
├── neutrino/             ← K_8 extension: 3 active + 3 sterile + PMNS + δ_CP
├── qed/, qcd/            ← Gauge-theory shims (R-ratio, DIS, exotic states)
├── qft/                  ← NWT Lagrangian (Paper 16), continuum field theory
├── topology/             ← K_7 Heffter embedding, colored Jones, torus knots
├── chemistry/            ← Aromaticity, NICS, C_60, SMILES → resonance energy
├── atomic/               ← Hydrogen chain (Bohr, Lyman α, 21 cm, Lamb shift)
├── dark_sector/          ← 98 GeV WIMP + LZ-2024 constraints
├── heron/                ← qiskit-runtime adapter for IBM Heron hardware
└── qpu/                  ← Vendor-neutral adapter (IBM, AWS Braket, simulator)
```

The "shim" pattern: each subpackage is a **view** of the same K_7 substrate. The constants come from `isa/`; the shim packages re-shape them for the conventions of a given physics domain. Cross-shim consistency is enforced by import-time identity assertions, not by convention.

## Shim reference

Per-shim documentation pages follow a consistent shape: TL;DR, common questions, prediction table, install + quick start, API by topic, worked examples, paper citations, cross-links.

### Core physics shims

| Shim | What it computes | Status |
|---|---|---|
| [`gravity`](shims/gravity.md) | Newton's G, BH thermodynamics, Kerr/NHEK, cosmogenesis | written |
| [`particles`](shims/particles.md) | 80-entry mass spectrum, decay constants, charge/stability | written |
| [`electroweak`](shims/electroweak.md) | Higgs VEV, sin²θ_W, G_F, Z widths, CKM, form factors | written |
| [`cosmology`](shims/cosmology.md) | η_B, Ω_b/Ω_c, Λ, CMB anisotropy axes | written |
| [`neutrino`](shims/neutrino.md) | K_8 extension: 3 active + 3 sterile + PMNS + δ_CP | written |
| [`qed`](shims/qed.md) | Compton, e⁺e⁻→μ⁺μ⁻, Møller, Bhabha; α running; a_e; Feynman diagrams | written |
| [`qcd`](shims/qcd.md) | SU(3) color algebra, α_s running, qq̄/qq/gg processes, confinement, exotic states | written |
| [`chemistry`](shims/chemistry.md) | Aromaticity, NICS sign, C_60 vibrational modes | written |
| [`atomic`](shims/atomic.md) | Hydrogen chain (Bohr, Lyman α, 21 cm, Lamb) | written |
| [`dark_sector`](shims/dark_sector.md) | 98 GeV WIMP, LZ-2024 constraints, Higgs portal | written |

### Infrastructure shims

| Shim | What it provides | Status |
|---|---|---|
| `isa` | Substrate constants + K_7 algebra (single source of truth) | [`README`](../nwt_substrate/isa/README.md) |
| `benchmarks` | 38 substrate-vs-experiment benchmarks | [`README`](../nwt_substrate/benchmarks/README.md) |
| `heron` | qiskit-runtime adapter for IBM Heron, exp11 sidereal program | [written](shims/heron.md) |
| `qpu` | Vendor-neutral QPU adapter (IBM, Braket, simulator) | [written](shims/qpu.md) |
| `topology` | K_7 Heffter embedding, colored Jones, torus knots | [written](shims/topology.md) |
| `qft` | NWT Lagrangian, continuum field theory (Paper 16) | [written](shims/qft.md) |

## Headline predictions

All derived at zero free parameters beyond `m_e`, `M_Pl`, `c`, `ℏ`. Click through to the shim docs for full derivations and PDG comparisons.

| Observable | Substrate formula | Accuracy | Shim |
|---|---|---|---|
| Fine structure constant | `1/α = 25π√3 + 1` | 7.6 ppm CODATA | [`isa`](../nwt_substrate/isa/README.md) |
| Newton's G | `G = (8/7)² α²¹ ℏc / m_e²` (LO) | −11 ppm CODATA (NNLO) | [`gravity`](shims/gravity.md) |
| Weinberg angle | `sin²θ_W = (2 + α)/9` | 0.06 % PDG | [`electroweak`](shims/electroweak.md) |
| Higgs VEV | `v_EW = 246.21 GeV` | 28 ppm PDG | [`electroweak`](shims/electroweak.md) |
| Baryon asymmetry | `η_B = (3/14)α⁴` | 0.38 % Planck | [`cosmology`](shims/cosmology.md) |
| Ω_b / Ω_c | `25α(1 + 3α)` | 0.0067 % Planck | [`cosmology`](shims/cosmology.md) |
| 80-entry mass spectrum | Paper 6 topological formula | <1 % median PDG | [`particles`](shims/particles.md) |
| Cosmological constant Λ | Substrate IR projection | 0.74 % Planck | [`cosmology`](shims/cosmology.md) |

## How NWT differs from neighbouring frameworks

- **Standard Model**: NWT doesn't replace it; it derives its ~19 free parameters from the substrate algebra at zero free parameters. Above the substrate level, standard QFT is the calculational framework.
- **String theory**: 10⁵⁰⁰+ vacua, requires extra dimensions. NWT has a single fixed substrate (K_7 on the Brieskorn-Poincaré sphere with Cl(0,7) algebra), no compactification, no landscape.
- **Loop quantum gravity**: Quantizes gravity directly. NWT derives gravity (Newton's G via Sakharov-induced) as one observable among many; gravity is emergent, not fundamental.
- **Wolfram Physics**: Also substrate-discrete, but uses generic hypergraph rewriting. NWT uses the specific K_7 + Cl(0,7) algebra and makes closed-form numerical predictions for SM observables.

See [`docs/FAQ.md`](FAQ.md) for atomic answers to these and other questions.

## Project archive

- **Repository**: <https://github.com/JimGalasyn/nwt-substrate>
- **Paper series**: <https://zenodo.org/communities/nwt>
- **Software DOI**: 10.5281/zenodo.20012027
- **License**: MIT
- **Author**: Jim Galasyn, <jim.galasyn@hotmail.com>

## How to cite

For software: see [`CITATION.cff`](../CITATION.cff). For physics claims: cite the relevant NWT paper plus the library version. See [`docs/code_division_policy.md`](code_division_policy.md) for the citation convention.
