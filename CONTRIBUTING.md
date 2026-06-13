# Contributing to nwt-substrate

> Thank you for considering a contribution. This is a physics library, not a generic Python project — the integrity of every numerical prediction depends on the codebase **not** introducing fitted parameters, magic numbers, or shortcuts. This guide explains what that means in practice and how to land a clean PR.

For AI coding agents (Claude Code, GitHub Copilot, Cursor, ChatGPT, etc.), see [`AGENTS.md`](AGENTS.md) — same rules, agent-targeted phrasing.

## TL;DR

1. **Open an issue first** for any non-trivial change so we can agree on scope before you write code.
2. **Install dev deps**: `pip install -e ".[test]"`.
3. **Run the suite**: `pytest` (~10 s, 1352 tests).
4. **Run the benchmarks**: `python -m nwt_substrate.benchmarks` (~100 ms, 38 substrate-vs-experiment checks).
5. **Add tests for new code** — coverage must not regress below the `codecov.yml` patch target (80 %).
6. **Add a CHANGELOG entry** under `[Unreleased]` describing your change.
7. **Don't change `isa.constants`** without a written justification — those integers are load-bearing across every shim.

## Common questions

### How is `nwt-substrate` different from a typical Python library?

It's the executable companion to a published [paper series](https://zenodo.org/communities/nwt). Every numerical claim in the papers is reproducible by calling a function here. That means the library has commitments most libraries don't:

- Every observable has a derivation that traces back to four substrate constants (`m_e`, `M_Pl`, `c`, `ℏ`).
- No "magic" numerical constants — every integer in the code (`7`, `21`, `8`, `3`, `25`, …) is a substrate primitive sourced from `nwt_substrate.isa.constants`.
- No fitted parameters. If you're tempted to tune a number to make a test pass, the test is telling you something — investigate before patching.
- Substrate identities are asserted at import time, so a change that breaks them will refuse to even load.

### What should I work on?

Good candidates:
- **Add a new benchmark** for an observable that doesn't yet appear in `nwt_substrate/benchmarks/`. Pick something with an established experimental value (PDG / CODATA / Planck / FLAG).
- **Improve test coverage** of a partially-tested module — see the [Codecov dashboard](https://codecov.io/gh/JimGalasyn/nwt-substrate).
- **Document a shim** — `docs/shims/<shim>.md`. The template is the four existing pages (`gravity`, `particles`, `electroweak`, `cosmology`).
- **Performance improvements** — replace Python loops with numpy einsum kernels where the substrate algebra allows.
- **Fix a bug** flagged by the test suite or by a paper-author review.

Avoid:
- New free parameters anywhere in the library.
- New shims (the seven existing ones cover essentially all SM + cosmology + chemistry).
- Vendoring `nwt_substrate/` into another repo — see [`docs/code_division_policy.md`](docs/code_division_policy.md).

### How long should I wait for review?

Issues get a first-pass response within ~1 week, PRs within ~2 weeks. If you haven't heard back in that window, ping the issue.

## Development setup

```bash
git clone https://github.com/JimGalasyn/nwt-substrate.git
cd nwt-substrate
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

Optional extras:
- `[heron]` — IBM Heron qiskit-runtime integration
- `[torch]` — torch backend for batched substrate kernels
- `[all]` — everything

Python 3.10, 3.11, 3.12 are supported. CI runs all three.

## Running tests

```bash
# Full suite (~10 s, 1352 tests)
pytest

# With coverage report
pytest --cov=nwt_substrate --cov-branch --cov-report=term

# A single file
pytest nwt_substrate/tests/test_particles_decay_constants.py

# A single test
pytest nwt_substrate/tests/test_particles_decay_constants.py::test_m_tau_substrate_matches_pdg_sub_percent
```

If you hit a `pytest` plugin-autoload conflict locally (some users have an `omegaconf` / `hydra` install that clashes), try:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -p pytest_cov nwt_substrate/tests/
```

## Running benchmarks

```bash
# Full table (~100 ms)
python -m nwt_substrate.benchmarks

# Summary only
python -m nwt_substrate.benchmarks --summary

# JSON output (CI artifact format)
python -m nwt_substrate.benchmarks --json
```

## Hard rules (substrate algebra integrity)

These rules exist because violations silently invalidate predictions across multiple shims. The test suite catches most violations; the rest catch you in review.

### 1. No fitted constants

Every numerical constant must derive from `nwt_substrate.isa.constants` (substrate primitives) or from an explicit physical anchor (CODATA `m_e`, PDG masses, Planck cosmological parameters). **If you are tempted to introduce a constant to make a test pass, the test is telling you something** — find the actual derivation before patching.

Bad:

```python
# We're off by ~3 %. Adding a fudge factor.
SOME_OBSERVABLE = 1.234 * 0.97  # don't
```

Good:

```python
# The substrate gives the factor 25 = H_V_SO7² = q_cinq² explicitly.
from nwt_substrate.isa.constants import H_V_SO7
SOME_OBSERVABLE = H_V_SO7 ** 2 * other_factor
```

### 2. No magic numbers in shim code

All structural integers (`7`, `21`, `8`, `3`, …) come from `nwt_substrate.isa.constants`. If you see a literal `7` or `21` in shim code, it should reference `isa.N_VERTICES_K7` or `isa.N_EDGES_K7`. The 92 substrate-identity tests will catch violations.

Bad:

```python
# 21 K_7 edges
amplitude = alpha ** 21  # magic number
```

Good:

```python
from nwt_substrate.isa.constants import N_EDGES_K7
amplitude = alpha ** N_EDGES_K7
```

### 3. No silent constant changes

`isa.constants` are import-time-asserted to satisfy structural identities (e.g., `N_EDGES_K7 == 21 == DIM_ADJ_SPIN7`, `4 + 3 == N_VERTICES_K7`). Changing one constant breaks an identity chain — the import will refuse. If you genuinely need to update a constant, open an issue first; we will most likely tell you "no, the constant is forced by the substrate," but if there's an actual CODATA update, we'll discuss the right migration.

### 4. No untested predictions

Every observable function should have a corresponding test that compares the substrate prediction against PDG / CODATA / Planck / FLAG. The test name should be `test_<observable>_matches_<source>`, e.g., `test_higgs_vev_matches_pdg`. The test must use a real experimental value, not a value fitted to the substrate.

### 5. Don't touch the paper series

Code that supports a specific paper figure or experimental run lives in the `null-worldtube-private` paper repo, **not** here. See [`docs/code_division_policy.md`](docs/code_division_policy.md) for the three-tier division (library / paper-glue / research-record). The library is the canonical source; paper repos import it, they don't duplicate it.

## Soft rules (style + workflow)

- **Prefer adding a benchmark to adding a one-off script.** New observables should land in `nwt_substrate/benchmarks/` with a `BenchmarkResult` so they appear in `run_all()` and the CI badge updates.
- **Prefer extending an existing shim to creating a new module.** The seven existing shims (chemistry, gravity, qed, qcd, particles, electroweak, heron) are the canonical decomposition.
- **Type hints everywhere.** Use `from __future__ import annotations` at the top of new files for lazy evaluation.
- **Dataclasses for results.** Don't return tuples or dicts when a `@dataclass(frozen=True)` would do.
- **numpy-first.** Use einsum kernels over explicit loops where the substrate algebra is naturally tensorial. Add a torch backend if you need GPU.
- **Match the existing code style.** Read a neighboring file before writing a new one.
- **Update `llms.txt` and `llms-full.txt` when adding new headline predictions.** AI search engines lift answers from these files.
- **Update `CHANGELOG.md`** under `[Unreleased]` for any user-visible change. (Maintainers: cutting a tagged release follows [`docs/RELEASING.md`](docs/RELEASING.md).)

## How to extend the library

### Adding a new observable

1. Decide which shim it belongs to.
2. Add the derivation as a function in that shim, taking inputs from `isa.constants` and explicit physical anchors (CODATA / PDG).
3. Add a `BenchmarkResult` in `nwt_substrate/benchmarks/compute_speed.py` that compares your prediction to the experimental value. Wire it into `run_all()` and export it from `benchmarks/__init__.py`.
4. Add a unit test in `nwt_substrate/tests/` that asserts the prediction matches the experimental value to the claimed precision.
5. If the new observable becomes a headline claim, update `llms.txt`, `llms-full.txt`, `docs/FAQ.md`, and the relevant `docs/shims/<shim>.md` page.
6. If it's a published-paper claim, add a paper citation in the function docstring.
7. Add a `CHANGELOG.md` entry under `[Unreleased]`.

### Adding a new shim

This is rare — the seven existing shims cover essentially all SM + cosmology + chemistry. If you genuinely need a new one:

1. Open an issue first.
2. Create `nwt_substrate/<new_shim>/` with `__init__.py` and `README.md`.
3. The shim must consume `isa.constants` (no internal magic numbers).
4. Implement a `substrate_breakdown()` function printing the substrate-identity table for the shim.
5. Add cross-shim consistency tests in `nwt_substrate/tests/test_isa_*.py`.
6. Write `docs/shims/<new_shim>.md` following the existing template.

### Modifying `isa.constants`

**Don't.** Changing a structural integer breaks substrate identities and can silently invalidate predictions across multiple shims. If you genuinely think a constant needs updating (e.g., new CODATA value), open an issue first — most "needs updating" candidates turn out to be forced by the substrate.

## PR workflow

1. **Open an issue** for non-trivial changes (anything beyond a typo, doc fix, or bug fix). Describe what you want to change and why.
2. **Fork + branch.** Branch names: `feature/<short>` or `fix/<short>`.
3. **Write tests first** where possible. Tests document the substrate identity you're enforcing.
4. **Commit messages.** Imperative mood, short summary line (≤ 72 chars), longer body if needed. Examples from the repo:
   - `particles.decay_constants: consolidate light + heavy + vector + B_c sectors`
   - `ci: pass CODECOV_TOKEN explicitly to codecov-action`
   - `tests: production-code coverage pass (gravity/nhek, heron/sidereal, qed+qcd/diagram)`
5. **Run the suite locally** before pushing: `pytest && python -m nwt_substrate.benchmarks`.
6. **Push and open a PR.** The PR template will ask for a summary, the relevant issue, and a checklist (tests added, CHANGELOG updated, llms.txt updated if relevant).
7. **CI must be green** before review:
   - `tests` workflow on Python 3.10, 3.11, 3.12
   - `benchmarks` workflow runs `run_all()` and verifies all benchmarks within their documented thresholds
   - Codecov patch coverage ≥ 80 %
8. **Iterate on review feedback.** We are friendly but firm about the substrate-integrity rules.

## Anti-patterns the test suite catches

These have been tried before and failed; the test suite or import-time assertions catch them in CI.

- Adding `from .constants import *` then redefining `N_EDGES_K7 = 22` locally — caught by `test_substrate_constants_consistent`.
- Hardcoding `21` in a shim instead of `isa.N_EDGES_K7` — caught by `test_no_magic_numbers_in_shims`.
- Adding a fitted parameter to make a test pass — caught by code review; the test docstring will tell you what to actually fix.
- Vendoring an old copy of `nwt_substrate/` into another repo — see [`docs/code_division_policy.md`](docs/code_division_policy.md).
- Importing from `analysis/` into the library — `analysis/` is sandbox; production code must not depend on it.

## Citation

If your contribution leads to a publication, cite both the library (`CITATION.cff`, DOI [10.5281/zenodo.20012027](https://doi.org/10.5281/zenodo.20012027)) and the relevant NWT paper(s). See [`docs/code_division_policy.md`](docs/code_division_policy.md) for the citation convention.

## Code of conduct

Be kind. Disagree respectfully. Substrate-integrity arguments come with citations, not invective. Reviews critique code, not authors.

## Getting help

- **Repository issues**: <https://github.com/JimGalasyn/nwt-substrate/issues>
- **Paper series**: <https://zenodo.org/communities/nwt>
- **Author**: Jim Galasyn, <jim.galasyn@hotmail.com>
- **FAQ**: [`docs/FAQ.md`](docs/FAQ.md) — atomic Q&A summaries
- **For physics questions**: read the relevant paper. For code questions: read the FAQ or open an issue.

## License

By contributing, you agree that your contributions will be licensed under the project's [MIT License](LICENSE).
