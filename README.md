# nwt-substrate

A substrate-algebraic computation library for Null Worldtube Theory (NWT).

`nwt-substrate` is the reference implementation of the substrate algebra
described in the NWT paper series: a Cl(0,7) octonion Clifford algebra
with K_7 graph state on the Heegaard torus of the Brieskorn-Poincaré
sphere S^3 / 2I, supporting particle / scattering / decay /
gravitational-coupling computations from a single internally consistent
codebase.

## Install

```bash
pip install nwt-substrate          # not yet on PyPI; for now:
pip install git+https://github.com/JimGalasyn/nwt-substrate.git
```

## Quick start

Particle masses from substrate quantum numbers:

```python
>>> import nwt_substrate as nwt
>>> p = nwt.particle("p")
>>> p.mass_pred
937.24...                             # MeV, Paper 6 mass formula
>>> p.J, p.Q, p.B
(0.5, 1, 1)
```

Connected-sum composition law for molecular bound states:

```python
>>> p, n = nwt.particle("p"), nwt.particle("n")
>>> d = nwt.compose(p, n, op="#", name="d", m_obs=1875.61)
>>> d.mass_pred                       # ~1874.48 MeV
>>> d.mass_residual                   # ~ -0.06 % vs PDG
```

Gravitational coupling from substrate alone:

```python
>>> from nwt_substrate.gravity import G_substrate_SI
>>> G_substrate_SI()                  # 6.674228e-11 m^3 kg^-1 s^-2
                                       # -11 ppm of CODATA, inside ±22 ppm
                                       # experimental error bar
```

Substrate diagrams:

```python
>>> fig = nwt.diagrams.figure_paper18_unified()
>>> fig.savefig("paper18_fig1.pdf")
```

## What's implemented

- **Particles** -- Paper 6 mass formula, charge via extended GMN, the
  full SM hadronic + leptonic + exotic catalog.
- **Compositions** -- knot connected-sum (#) for molecular bound states
  (deuteron, X(3872), Pc family), Hopf-link with Λ_QCD = 313 MeV per
  crossing for nuclear / strongly-bound exotic regimes.
- **Walk-phase scattering** -- substrate-algebraic Compton (matches
  Klein-Nishina to 1e-9), Møller / Bhabha, V-A muon decay matching
  Sargent rate, neutron decay with g_A = 1.27.
- **Gauge-theory shims** -- nwt.qed, nwt.qcd (incl. gg→gg), nwt.electroweak
  (Z resonance + chiral couplings), nwt.qft (Lagrangian view),
  nwt.string (string-theoretic view), nwt.gravity (Sakharov-induced G).
- **Heron experiments** -- qiskit-runtime interface and an experiment
  registry for IBM Heron processors.  Supports Experiments 4 / 5 / 9
  / 10 / 11 from the paper series.
- **Diagrams** -- programmatic figure factories for the canonical
  substrate visualisations (torus knots, K_7 traversals, Heegaard-torus
  unification).

## Tests

```bash
pytest nwt_substrate/tests/ -q
# 341 passed in ~5s
```

## Citation

If you use this library in a publication, please cite both:

- The relevant NWT paper(s) -- typically one of
  [Paper 14--19](https://zenodo.org/communities/nwt) for the result
  you're using.
- The library Zenodo record (auto-archived per release):

```bibtex
@software{nwt_substrate,
  author       = {Galasyn, Jim and others},
  title        = {{nwt-substrate}: a substrate-algebraic computation
                  library for Null Worldtube Theory},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.PLACEHOLDER}
}
```

A `CITATION.cff` is included in this repo for tools that auto-resolve
software citations.

## Papers

The library implements the computations described in:

- **Paper 6** -- topological mass formula (1.06 % median residual on the
  24-particle compendium).
- **Paper 14** -- α^(21/2) heptafoil amplitude.
- **Paper 15** -- Wilson amplitude on K_7 graph state.
- **Paper 16** -- NWT three-field Lagrangian (BPS critical coupling).
- **Paper 17** -- m_e / m_Pl closed form: G to -11 ppm CODATA (inside the ±22 ppm experimental band).
- **Paper 18** -- Sakharov-induced Einstein gravity from substrate
  matter sector.  *Includes the canonical "Heegaard torus, two
  sectors" figure rendered by `nwt.diagrams.figure_paper18_unified()`.*
- **Paper 19** -- substrate monism via library demonstration.

The Zenodo community for the full series is at
https://zenodo.org/communities/nwt (collected DOIs).

## Status

Pre-1.0 software.  API surface is stable across the modules listed
above (particles, compositions, walk_phase, gauge shims, gravity,
diagrams) but minor breaking changes may still occur; we aim for
semver discipline post-1.0.

The main private development monorepo, where new analyses and paper
drafts live before promotion, is `null-worldtube-private` (not
public).  Polished analyses and paper-supporting computations are
promoted to this repo; exploratory work stays private.

## Contributing

Issues and pull requests welcome.  Please run the test suite
(`pytest nwt_substrate/tests/`) before submitting, and include a
short note describing the physics motivation for any new feature
(this is a physics library; please don't add tooling that has no
substrate-algebraic content).

## License

MIT.  See [LICENSE](./LICENSE).
