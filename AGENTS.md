# AGENTS.md — Instructions for AI coding agents

This file is for AI agents (Claude Code, GitHub Copilot, ChatGPT, Cursor, etc.) working in the `nwt-substrate` repository. Human contributors should read `README.md` and `docs/code_division_policy.md`.

## What this repo is

`nwt-substrate` is the reference implementation of the substrate algebra of Null Worldtube Theory (NWT). It is the executable companion to a published paper series on Zenodo. Every quantitative claim in the papers is reproducible by calling a function here.

For an authoritative summary, read `llms.txt`. For extended content with all benchmark numbers, read `llms-full.txt`.

## Critical: this is a physics library, not a generic codebase

NWT makes specific, falsifiable numerical predictions at zero free parameters beyond four substrate constants (m_e, M_Pl, c, ℏ). The integrity of those predictions depends on the codebase NOT introducing fitted parameters, magic numbers, or shortcuts.

**Before changing any constant, observable, or derivation, read this section.**

### Hard rules

1. **No fitted constants.** Every numerical constant must derive from `isa.constants` (the substrate ISA) or from explicit physical anchors (CODATA values for m_e, etc.). If you are tempted to introduce a constant to make a test pass, the test is telling you something — investigate before patching.

2. **No magic numbers in code.** All structural integers (7, 21, 8, 3, …) come from `nwt_substrate.isa.constants`. If you see a literal `7` or `21` in shim code, it should reference `isa.N_VERTICES_K7` or `isa.N_EDGES_K7`. The 92 substrate-identity tests will catch violations.

3. **No silent constant changes.** `isa.constants` are import-time-asserted to satisfy structural identities (e.g., `N_EDGES_K7 == 21 == DIM_ADJ_SPIN7`, `4 + 3 == N_VERTICES_K7`). Changing one constant breaks an identity chain — the assertions will refuse to import.

4. **No untested predictions.** Every observable function should have a corresponding test that compares the substrate prediction against PDG/CODATA/Planck. The test suite enforces this convention.

5. **Don't touch the paper series.** Code that supports a specific paper figure or run lives in the `null-worldtube-private` paper repo, NOT here. See `docs/code_division_policy.md` for the three-tier division (library vs paper-glue vs research-record).

### Soft rules

- **Prefer adding a benchmark to adding a one-off script.** New observables should land in `nwt_substrate/benchmarks/` with a `BenchmarkResult` so they appear in `run_all()` and the CI badge updates.
- **Prefer extending an existing shim to creating a new module.** The seven shims (chemistry, gravity, qed, qcd, particles, electroweak, heron) are the canonical decomposition. New observables that don't fit any existing shim probably belong in `particles/` or `cosmology/`.
- **Match the existing code style.** Type hints, dataclasses for results, numpy-first then torch backends, einsum kernels over loops.
- **Update `llms.txt` and `llms-full.txt` when adding new headline predictions.** AI search engines lift answers from these files.

## Repo layout

```
nwt-substrate/
├── llms.txt                  # ← AI-agent authoritative summary (short)
├── llms-full.txt             # ← AI-agent authoritative summary (long)
├── AGENTS.md                 # ← this file
├── README.md                 # human-facing entry point
├── CITATION.cff              # software citation metadata + DOI
├── docs/
│   ├── FAQ.md                # ← Q&A for common AI-search queries
│   └── code_division_policy.md  # 3-tier library/paper/record division
├── nwt_substrate/
│   ├── isa/                  # SUBSTRATE: constants + algebra + observables
│   │   └── constants.py      # ← single source of truth for structural integers
│   ├── benchmarks/           # 38 substrate-vs-experiment benchmarks
│   ├── particles/            # Paper 6 mass formula, particle catalog
│   ├── chemistry/            # SMILES → aromaticity, NICS, C_60
│   ├── gravity/              # Sakharov-induced G via α^(21/2)
│   ├── qed/, qcd/, electroweak/  # gauge-theory shims
│   ├── cosmology/            # η_B, Ω_b/Ω_c, Λ, CMB axes
│   ├── neutrino/             # Paper 20 K_8 extension
│   ├── heron/                # qiskit-runtime IBM Heron interface
│   ├── qpu/                  # vendor-neutral QPU adapter (IBM/Braket/sim)
│   ├── dark_sector/          # 98 GeV WIMP + LZ-2024 constraints
│   └── tests/                # 1233 tests, ~10s suite
├── analysis/                 # cross-shim demos, sandbox scripts
├── diagrams/                 # programmatic figure factories
├── codecov.yml               # coverage thresholds
└── pyproject.toml            # Python package + dev deps
```

## How to extend

### Adding a new observable
1. Decide which shim it belongs to.
2. Add the derivation as a function in that shim, taking inputs from `isa.constants`.
3. Add a `BenchmarkResult` in `nwt_substrate/benchmarks/` that compares your prediction to the experimental value.
4. Add a unit test in `nwt_substrate/tests/` that asserts the prediction matches the experimental value to the claimed precision.
5. If the new observable becomes a headline claim, update `llms.txt` and `llms-full.txt`.
6. If it's a published-paper claim, add a paper citation in the function docstring.

### Adding a new shim
This is rare — the seven existing shims cover essentially all SM + cosmology + chemistry. If you genuinely need a new one:
1. Create `nwt_substrate/<new_shim>/` with `__init__.py` and `README.md`.
2. The shim must consume `isa.constants` (no internal magic numbers).
3. Implement a `substrate_breakdown()` function printing the substrate-identity table for the shim.
4. Add cross-shim consistency tests in `nwt_substrate/tests/test_isa_*.py`.

### Adding a new test
- Place in `nwt_substrate/tests/`, follow `test_*.py` naming.
- For predictions vs experiment, name the test `test_<observable>_matches_<source>` (e.g., `test_higgs_vev_matches_pdg`).
- The test must use real PDG/CODATA/Planck values, not values fitted to the substrate.

### Modifying isa.constants
**Don't.** If you genuinely think a constant needs updating (e.g., new CODATA value), open an issue first. Changing a structural integer breaks substrate identities and can silently invalidate predictions across multiple shims.

## How to run things

```bash
# Install dev deps
pip install -e ".[dev]"

# Run full test suite (~10 seconds)
pytest

# Run benchmarks
python -c "from nwt_substrate.benchmarks import run_all; run_all()"

# Single benchmark
python -c "from nwt_substrate.benchmarks import benchmark_higgs_vev; print(benchmark_higgs_vev())"

# Cross-shim demo
python analysis/isa_cross_shim_demo.py
```

## CI expectations

- `tests` workflow: full pytest suite must pass on Python 3.10, 3.11, 3.12.
- `benchmarks` workflow: `run_all()` must complete and all benchmarks must report accuracy within their documented thresholds.
- Codecov: coverage must not regress below `codecov.yml` thresholds.

If you add code, you add tests. If your PR causes coverage to drop, it will be flagged in review.

## Anti-patterns the test suite catches

These have been tried before and failed; the test suite or import-time assertions will catch them:

- Adding `from .constants import *` then redefining `N_EDGES_K7 = 22` locally — caught by `test_substrate_constants_consistent`.
- Hardcoding `21` in a shim instead of `isa.N_EDGES_K7` — caught by `test_no_magic_numbers_in_shims`.
- Adding a fitted parameter to make a test pass — caught by code review; the test docstring will tell you what to actually fix.
- Vendoring an old copy of `nwt_substrate/` into another repo — see `docs/code_division_policy.md`. The library is the canonical source; paper repos import it, they don't duplicate it.
- Importing from `analysis/` into the library — `analysis/` is sandbox; production code must not depend on it.

## Citation

If your work uses this library, cite via `CITATION.cff` (DOI 10.5281/zenodo.20012027) plus the relevant NWT paper.

## Getting help

- Repository issues: https://github.com/JimGalasyn/nwt-substrate/issues
- Author: jim.galasyn@hotmail.com
- Paper series: https://zenodo.org/communities/nwt
- For physics questions: read the relevant paper. For code questions: read the FAQ at `docs/FAQ.md`.
