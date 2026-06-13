# Releasing `nwt-substrate`

The maintainer runbook for cutting a tagged, Zenodo-archived release. (Contributor
workflow — how to land a change — lives in [`CONTRIBUTING.md`](../CONTRIBUTING.md);
this is the orthogonal "how to ship a version" process.)

Versions follow [SemVer](https://semver.org/): **minor** bump for new sectors /
observables / modules, **patch** for fixes and corrections that preserve the
public surface, value-preserving throughout (no fitted constants — see
`AGENTS.md`).

## Key facts (read once)

- **The version is git-tag-driven** (`setuptools-scm`). There is **no version
  string to bump in code** — `nwt_substrate.__version__` comes from the tag.
  The *only* manually-edited version field is `CITATION.cff`.
- **`main` is protected** by the "Copilot review for default branch" ruleset
  (requires passing **CodeQL** + tests, a Copilot review, and code-quality).
  Every change — including the release prep and the DOI backfill — goes through
  a **PR**; you cannot push to `main` directly, and `--admin` will not bypass a
  pending CodeQL result. If CodeQL is ever missing, merges block until
  `.github/workflows/codeql.yml` exists.
- **Zenodo mints the DOI on a published GitHub _Release_, not on a bare tag.**
  Pushing the tag alone does nothing; you must create the Release.
- Concept DOI **`10.5281/zenodo.20012027`** resolves to the latest version and
  never changes; each version gets its own version DOI, backfilled after release.
- **PyPI publishes from the same Release event** via trusted publishing (OIDC, no
  stored token — `.github/workflows/publish-pypi.yml`). The tag drives the version
  (`setuptools-scm`), so there's nothing to bump. A one-time *pending publisher*
  registration on PyPI is needed before the first release (see step 4).

## Steps

### 1. Decide it's release-worthy & pick the version

```bash
git log "$(git describe --tags --abbrev=0)"..main --oneline   # what's unreleased
```

New modules/observables → minor (e.g. `v0.5.0`); fixes/corrections → patch. Make
sure `main` CI is green first.

### 2. Prep PR — docs + metadata (no code)

Branch `release/vX.Y.Z`. Edit:

- **`CHANGELOG.md`** — move the `[Unreleased]` content into a new
  `## [X.Y.Z] - YYYY-MM-DD` section with a one-paragraph narrative intro; leave a
  fresh empty `## [Unreleased]` on top. Use a DOI placeholder
  (`_backfilled after the GitHub release_`) — it isn't minted yet.
- **`docs/releases/vX.Y.Z.md`** — narrative release notes (copy the structure of
  the previous one: title, DOI line, Install/upgrade, "The story", one section
  per headline item, "What's next", "See also").
- **`CITATION.cff`** — bump `version:` and `date-released:`. Leave `identifiers:`
  for the backfill.
- **`README.md`** — the four things that reliably go stale:
  1. the "current release **vX.Y.Z**, DATE" line near the top;
  2. the **test count + timing** (two places: the headline bullet and the
     `## Tests + coverage` snippet) — match the latest `main` CI
     (`N passed … in Ts`);
  3. the **"What's implemented"** list — add a bullet per new sector, tagged
     `(vX.Y new)`;
  4. the `## Status` paragraph and the `## Changelog and releases` "latest is…"
     pointer.

Get the current test count from CI:

```bash
RID=$(gh run list --branch main --workflow tests --limit 1 --json databaseId -q '.[0].databaseId')
gh run view "$RID" --log | grep -oE "[0-9]+ passed.*" | tail -1
```

Open the PR, let CI go green, merge (squash).

### 3. Tag the release commit (annotated)

```bash
git checkout main && git pull --ff-only
git tag -a vX.Y.Z <merge-commit> -m "nwt-substrate vX.Y.Z — <one-line theme>

<short body: the headline sectors/fixes>"
git push origin vX.Y.Z
```

Tag subject convention: `nwt-substrate vX.Y.Z — <theme>` (matches prior tags).

### 4. Publish the GitHub Release (this triggers Zenodo **and PyPI**)

```bash
gh release create vX.Y.Z \
  --title "vX.Y.Z — <theme>" \
  --notes-file docs/releases/vX.Y.Z.md
```

Publishing the Release fires two workflows off the one event: Zenodo mints the
version DOI, and `.github/workflows/publish-pypi.yml` builds + publishes the
sdist/wheel to PyPI via **trusted publishing** (OIDC — no token). Confirm both:

```bash
gh run list --workflow publish-pypi.yml --limit 1   # build+publish green?
curl -s https://pypi.org/pypi/nwt-substrate/json -o /dev/null -w "PyPI HTTP %{http_code}\n"
```

> **One-time PyPI setup (before the FIRST release only).** The project doesn't
> exist on PyPI yet, so register a **pending publisher**: pypi.org → Account →
> Publishing → *Add a pending publisher* with Project `nwt-substrate`, Owner
> `JimGalasyn`, Repository `nwt-substrate`, Workflow `publish-pypi.yml`,
> Environment `pypi`. Also create a GitHub Environment named `pypi`
> (repo Settings → Environments). After the first successful publish PyPI
> converts it to a normal trusted publisher; later releases need nothing.

### 5. Backfill the Zenodo version DOI

Zenodo mints it within a few minutes of the Release publishing. Find it via the
concept record (sanity-check that the returned v(N-1) DOI matches the last
release):

```bash
# all versions under the concept, newest first:
curl -s "https://zenodo.org/api/records/?q=conceptrecid:20012027&all_versions=true&sort=mostrecent&size=5" \
  | python3 -c "import json,sys; [print(r['metadata']['version'], r['doi']) for r in json.load(sys.stdin)['hits']['hits']]"
```

Then, on a `chore/backfill-vX.Y.Z-doi` branch, replace the placeholders with
`Version DOI [10.5281/zenodo.NNNNN](https://doi.org/10.5281/zenodo.NNNNN)` in:

- `CHANGELOG.md` (the `[X.Y.Z]` header line),
- `docs/releases/vX.Y.Z.md` (the header DOI line),
- `CITATION.cff` (prepend a new `identifiers:` entry, most-recent-first).

Open the PR, CI green, merge. Then refresh the published Release body too:

```bash
gh release edit vX.Y.Z --notes-file docs/releases/vX.Y.Z.md
```

### 6. Downstream (optional)

- **`jax-solitons`** pins this package as its cross-engine oracle by commit SHA in
  the `oracle` extra. After a release that touches `solitons.faddeev` (or its
  deps), repoint that pin from the bare SHA to `@vX.Y.Z` — a citable tag — and
  bump it intentionally so the equivalence gate re-validates against the new oracle.
- If the release changes a headline prediction, update `llms.txt` / `llms-full.txt`.

## Quick checklist

- [ ] `main` green; version chosen (semver)
- [ ] Prep PR: CHANGELOG + `docs/releases/vX.Y.Z.md` + CITATION (version/date) + README (4 spots) — merged
- [ ] Annotated tag pushed (`nwt-substrate vX.Y.Z — theme`)
- [ ] (first release only) PyPI pending publisher + `pypi` GitHub Environment registered
- [ ] GitHub Release published (triggers Zenodo **and** PyPI)
- [ ] PyPI publish workflow green; package resolves on pypi.org
- [ ] DOI backfilled (CHANGELOG + release notes + CITATION) via PR; Release body refreshed
- [ ] Downstream pins / `llms.txt` updated if applicable
