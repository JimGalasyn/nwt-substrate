# Code-division policy — library vs paper vs record

**Status: RATIFIED (VV draft + NWT sign-off, 2026-05-25).** Fold into `CONTRIBUTING`.
A cross-session governance proposal for how code is divided between
`nwt-substrate` (this library), the `null-worldtube` paper repos, and the
research record.

## The problem

Code has been pushed to the paper repos to support the papers. As the program
matures, "the supporting code" splits into three genuinely different things,
and conflating them causes drift (ad-hoc constant redefinitions), bit-rot, and
split-brain (e.g. a stale vendored `nwt_substrate/` copy diverging from the
canonical repo). We want one durable, citable source of truth for the *physics*,
without bloating the library with one-off paper glue.

## Three tiers

1. **Library — `nwt-substrate`.** The *canonical physics*: constants
   (`isa.constants`), algebra, and the derived-observable functions
   (`cosmology.eta_B`, the K₈ code, bridge geometry, …). Tested, versioned,
   CI'd, Zenodo-DOI'd. **Every quantitative claim a paper makes cites the
   library, by version/DOI.**
2. **Paper-glue — `null-worldtube` `analysis/`.** Paper-specific
   *orchestration*: figure generation, data fetch/reduction (HEALPix, Gaia),
   the run that produced a specific figure, a one-off falsifier Monte Carlo, and
   **hardware experiment drivers** (the Steane/Heron/Braket submit–poll–decode
   harnesses, which import the `qpu` library). These **import the physics from
   the library** and add only plotting/data/glue/orchestration. Cited only for
   figure or run reproduction.
3. **Research record.** Eliminated routes, diagnostics, superseded attempts
   (e.g. the α_KL radius route-eliminations). Preserved as the record; **never
   cited as a paper's method, never in the library.**

## The discriminator

> **Reusable + canonical + mature?** → library.
> **Paper-specific orchestration?** → paper repo (thin, imports the library).
> **Exploration / superseded?** → record only.

### Litmus test (operational)

> *If you deleted the library, how much physics would the paper script still
> contain?* Ideally ~none — just orchestration. If a paper script re-derives α,
> re-implements an observable, or rebuilds a code, that physics is a **migration
> candidate**: promote it to the library, reduce the script to a thin caller.

## What does NOT belong in the library (deliberate carve-outs)

- **One-off, non-reusable calcs** — physics, but not *canonical/reusable*; a
  bespoke single-paper Monte Carlo is paper-glue. Don't saddle the library with
  a tests+DOI maintenance contract on things that don't warrant it.
- **In-progress / speculative results** — keep in the paper/exploration tier
  until they stabilize. (Live example: the α_KL disk radius, 4 derivation routes
  eliminated, is *out* of the library on purpose.) Promote only mature,
  load-bearing results.

## Process

- Cite the library **by version + Zenodo DOI** for quantitative claims; cite
  paper-glue scripts only for figure reproduction.
- Paper `analysis/` = **import-from-library + figures/data**; no re-defined
  constants (import from `isa.constants`).
- **No `nwt_substrate/` package code in the paper repos — ever.** Every library
  change lands in the canonical `nwt-substrate` repo; paper repos only consume
  it. (The 2026-05-25 split-brain's root cause was not just vendoring but
  *different sessions committing library code to different repos*.)
- **Treat the library as a pinned dependency, never vendor it** —
  `pip install nwt-substrate==X` or a submodule pinned to a tag. (The 2026-05-25
  split-brain — a vendored copy at 3b6a178 vs canonical dc84602, with an
  egg-info resolving to the stale copy — is the failure mode this prevents.)
- **Transition state (2026-05-25):** the library is `0.1.x.dev` and un-DOI'd, and
  the paper `analysis/` scripts currently reach it via a `sys.path.insert(...)`
  hack. Before "cite by version + DOI" and "pinned dependency" are real, two
  concrete steps: (i) cut a tagged `nwt-substrate` release + Zenodo DOI; (ii)
  migrate the analysis scripts off `sys.path.insert` onto the installed package.
- **Triage the existing `analysis/` backlog** against the discriminator:
  promote mature physics (started with `cosmology.eta_B`; Ω_b/Ω_c and Λ_cc
  next), leave the glue, archive the eliminated routes.

## First application: the 21a / 21b / 22a bundle

Use the bundle as the first test of the policy:
- Every quantitative claim in 21a/21b/22a cites `nwt-substrate` (version + DOI).
- The cosmology observables 22a cites (η_B, Ω_b/Ω_c, Λ_cc, inflation sector,
  v_EW-from-M_Pl, K₈ tower) live in `cosmology.*`, importing `isa.constants`.
- 22a's figure/data scripts stay in the paper repo, importing the library.
- Anything still re-deriving physics in a paper script before release is a
  promote-or-justify item.
- **Release prerequisite:** the bundle cannot ship a "cite by version + DOI"
  reference until `nwt-substrate` has a tagged release + Zenodo DOI — cutting one
  is now a gating item for the 21a/21b/22a release.
