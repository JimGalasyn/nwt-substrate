# nwt_substrate.chemistry

> Substrate-algebraic chemistry: **O(1) / O(N) group-theoretic predictions** for observables that DFT/CCSD compute at O(N³)–O(N⁷). **Tier A** is closed-form (aromaticity class, C_60 orbit/vibrational counts, McKay coordinations, Wade SEP counts, NICS sign, irrep labels — the substrate algebra gives the answer outright). **Tier B** is skeleton + 1 calibration (aromatic resonance energies via the Hopf-pair rule — 5/5 linear acenes within 6 % RMS = 2.91 % from one benzene calibration). Includes a minimal **SMILES front-end**, π Hückel **bond orders**, **pericyclic** Woodward-Hoffmann selection rules, **polyhedral/Wade** borane rules, **molecular knots**, **transition-metal** electron counts, **NMR** ring-current signs, and a substrate-vs-DFT timing **benchmark** (speedups **1e6–4e11** demonstrated).

[← Back to index](../index.md) · Source: [`nwt_substrate/chemistry/`](../../nwt_substrate/chemistry/) · Papers: [6](https://zenodo.org/records/15376291) (carrier-knot `n_q` table), [13](https://zenodo.org/records/19635239) (Spin(7) rep-class ladder), [community series](https://zenodo.org/communities/nwt)

## Common questions

### What chemistry can a particle-physics substrate predict?

The same group-theory invariants that fix the NWT particle catalog (binary polyhedral groups `2T/2O/2I`, the K_7 graph, the Spin(7) representation-class ladder) also fix a surprising amount of structural chemistry — *without any electronic-structure calculation*. C_60 **is** the binary-icosahedral `2I` canonical vertex orbit, so all its face/edge/mode counts are pure group theory. Aromaticity is a Hopf-pair parity. Stable transition-metal electron counts `{18, 16, 14, 12, 32}` are `N_EDGES_K7 − k` with `k` walking the Spin(7) ladder `{3, 5, 7, 9}`. None of these need DFT.

### Tier A vs Tier B?

**Tier A — closed form (substrate O(1), DFT O(Nᵏ)).** The algebra gives the answer outright: magic numbers, irrep labels, selection rules, orbit sizes, allowed coordinations, NICS *sign*, Wade SEP counts. Speedup 10⁹–10¹¹.

**Tier B — substrate skeleton + 1 calibration (substrate O(N), DFT O(Nᵏ)).** The substrate fixes the *shape*; one experimental number sets the magnitude; everything else is algebraic. The flagship is aromatic resonance energy: count Hopf pairs, multiply by a single benzene-calibrated `kcal/mol` per pair, add a K_7-toroidal correction for hub ring-sets. 5/5 linear acenes land within 6 % (RMS 2.91 %).

### How does the Hopf-pair aromaticity rule work?

Each π pair = 2 electrons = one Hopf link in the carrier-knot picture. The Hopf-link parity classifies the ring:

```
odd Hopf parity   (4n+2 π electrons) → AROMATIC      (Hückel)
even Hopf parity  (4n   π electrons) → ANTI-AROMATIC
single ring + 1 Möbius twist          → MÖBIUS aromatic
```

`aromaticity_class("benzene")` → `parity='odd'`, `n_hopf_links=3` → aromatic. The *same* parity rule re-appears for pericyclic transition states (`pericyclic.selection_rule`) and for NICS ring-current signs (`nmr.nics_sign_from_hopf_parity`) — one substrate rule, three chemistry contexts.

### Can I feed it SMILES?

Yes. `parse_smiles` builds an atom/bond graph (organic subset, lowercase = aromatic, ring-closure digits, branches, charges), `find_aromatic_ring_systems` extracts conjugated cycles via `networkx` cycle-basis, and `smiles_to_aromaticity` / `smiles_resonance_energy[_batch]` apply the Hopf-pair rule to *arbitrary* molecules, not just the curated table. `ring_system_to_so7_adjacency` maps a ring system onto a 7×7 K_7-style adjacency for the substrate accounting.

### How big are the DFT speedups?

Substrate cost is measured by a real timer; DFT cost is estimated from published O(N³)–O(N⁷) scaling. From `benchmark.report_suite()`: aromaticity classification **2.4e9×**, C_60 vibrational modes **2.7e9×**, NICS sign rule **6.7e9×**, transition-metal electron count **4.6e9×**, molecular-knot accessibility **4.1e11×** (DFT-NICS estimate vs a sub-microsecond table lookup). The cheapest combinatorial counts (deltahedra) are only ~1e5× because the DFT baseline there is itself trivial.

## Prediction table

| Observable | Substrate route | Substrate value | Accuracy |
|---|---|---|---|
| Aromaticity classification | Hopf-pair parity (4n+2 / 4n) | benzene aromatic, cyclobutadiene anti-aromatic | 20/20 Hückel/Möbius matches |
| Aromatic RE (linear acenes) | Tier B: Hopf pairs × 1 benzene calibration | benzene 36, naphthalene 60, anthracene 84, coronene 144 kcal/mol | 5/5 within 6 %, RMS 2.91 % |
| C_60 combinatorics | `2I` orbit `\|A_5\|`, `\|2I\|/\|D_n\|` | 60 vertices, 12 pentagons, 20 hexagons, 90 edges | exact |
| C_60 vibrational modes | `I_h` irrep decomposition | 174 total, 4 IR (T_1u), 10 Raman (2 A_g + 8 H_g), 32 silent | exact vs IR/Raman spectroscopy |
| C_60 anion magic states | `2I` irrep filling | n=3 (A_3 SC), n=6 (hexafulleride), n=9/12 predicted | empirical match on n=3, n=6 |
| McKay coordinations | binary-polyhedral orbit sizes | allowed `{1,2,3,4,6,12}`; 5 = fluxional | 26/26 closed-shell match |
| Benzene π bond order | cyclic Hückel, all-bonding | 2/3 = `2/RANK_SO7` (= 0.6667) | exact rational |
| Wade SEP / borane class | `n_seps − n_vertices` offset | closo (6,7), nido (5,7), arachno (4,7), hypho (4,8) | 8/8 closo set exact |
| Closo 3D-aromaticity ladder | Spin(7) rep-class `(n, n+1)` | B_5–B_8 = (5,6)(6,7)(7,8)(8,9) | 5/8 double-canonical, p_rand ≈ 1.5e-4 |
| Transition-metal stable counts | `N_EDGES_K7 − k`, `k ∈ {3,5,7,9}` | 18, 16, 14, 12 (+ 32 f-block = 35 − 3) | reproduces empirical series |
| Molecular-knot accessibility | carrier `n_q` table | 5_1 sweet-spot (n_q=5), 6_1 gap (n_q=6) | P1 + P2 confirmed |
| NICS sign | Hopf-pair parity | benzene diatropic, cyclobutadiene paratropic | 14/14 reference set |

Each row is asserted in the chemistry test suite: [`test_chemistry_smoke.py`](../../nwt_substrate/tests/test_chemistry_smoke.py), [`test_chemistry_bond_orders.py`](../../nwt_substrate/tests/test_chemistry_bond_orders.py), [`test_chemistry_smiles.py`](../../nwt_substrate/tests/test_chemistry_smiles.py), [`test_chemistry_pericyclic.py`](../../nwt_substrate/tests/test_chemistry_pericyclic.py), [`test_chemistry_polyhedral.py`](../../nwt_substrate/tests/test_chemistry_polyhedral.py), [`test_chemistry_molecular_knots.py`](../../nwt_substrate/tests/test_chemistry_molecular_knots.py), [`test_chemistry_transition_metal.py`](../../nwt_substrate/tests/test_chemistry_transition_metal.py), [`test_chemistry_nmr.py`](../../nwt_substrate/tests/test_chemistry_nmr.py).

## Quick start

```python
import nwt_substrate.chemistry as chem

# Tier A — closed-form classifications
chem.aromaticity_class("benzene")
# → AromaticityResult(classification='aromatic', n_pi_pairs=3,
#                     n_hopf_links=3, parity='odd', ...)

chem.fullerene_orbit_counts(60)
# → FullereneOrbit(n_vertices=60, n_pentagons=12, n_hexagons=20,
#                  n_edges=90, point_group='I_h', is_2I_orbit=True, ...)

chem.c60_vibrational_summary()
# → VibrationalSummary(total_modes=174, ir_active_modes=4,
#                      raman_active_modes=10, silent_modes=32, ...)

chem.allowed_coordinations()          # → (1, 2, 3, 4, 6, 12)
chem.is_admissible_coordination(5)    # → False  (fluxional / Berry pseudorotation)
chem.bond_order("benzene")            # → 1.5

# Tier B — substrate skeleton + 1 calibration
chem.aromatic_resonance_energy("naphthalene", calibration_kcal=36.0)
# → ResonanceEnergy(kcal_per_mol=180.0, n_hopf_pairs=5, ...)

# SMILES front-end
chem.smiles_to_aromaticity("c1ccccc1").classification     # → 'aromatic'
chem.smiles_resonance_energy("c1ccc2ccccc2c1")            # → 60.0 (naphthalene)

# Substrate-vs-DFT benchmark
chem.benchmark.compare_to_dft("aromaticity", n_molecules=20).speedup
# → ~2.5e9
```

## API by topic

### Aromaticity (`chemistry.aromaticity`)

| Function / object | Returns |
|---|---|
| `aromaticity_class(name)` | `AromaticityResult` — classification, `n_pi_pairs`, `n_hopf_links`, `parity` |
| `bond_order(name, bond_type="C-C")` | π bond order (benzene → 1.5) |
| `aromatic_resonance_energy(name, calibration_kcal=12.0)` | `ResonanceEnergy` — Tier-B RE in kcal/mol |
| `AromaticityResult`, `ResonanceEnergy` | result dataclasses |

### SMILES (`chemistry.smiles`)

| Function / object | Returns |
|---|---|
| `parse_smiles(smiles)` | `Molecule` (atom/bond graph) |
| `smiles_to_aromaticity(smiles, mobius_twists=0)` | `SmilesAromaticityResult` |
| `smiles_resonance_energy(smiles, calibration_kcal=12.0, ...)` | RE in kcal/mol |
| `smiles_resonance_energy_batch(smiles_list, ...)` | array / details of REs |
| `ring_system_to_so7_adjacency(ring_system)` | `(7×7 ndarray, n)` K_7-style adjacency |
| `Molecule`, `Atom`, `Bond`, `AromaticRingSystem`, `SmilesAromaticityResult` | graph + result dataclasses |

### Fullerenes (`chemistry.fullerenes`)

| Function / object | Returns |
|---|---|
| `fullerene_orbit_counts(n_vertices)` | `FullereneOrbit` — vertices/pentagons/hexagons/edges from `2I` orbit |
| `c60_anion_magic_states()` | list of `AnionState` (n=3 SC, n=6 hexafulleride, n=9/12 predicted) |
| `c60_combinatorial_summary()` | dict of all C_60 combinatorial counts |
| `FullereneOrbit`, `AnionState` | result dataclasses |

### Vibrational (`chemistry.vibrational`)

| Function / object | Returns |
|---|---|
| `c60_vibrational_summary()` | `VibrationalSummary` — 174 modes, 4 IR, 10 Raman, 32 silent + `I_h` decomposition |
| `point_group_mode_count(n_atoms, point_group="I_h")` | dict: `total_modes`, `ir_active`, `raman_active`, `silent`, `decomposition` |
| `VibrationalSummary` | result dataclass |

### McKay coordinations (`chemistry.mckay`)

| Function / object | Returns |
|---|---|
| `allowed_coordinations()` | `(1, 2, 3, 4, 6, 12)` — McKay-admissible only |
| `is_admissible_coordination(n)` | bool (5 → False, fluxional) |
| `mckay_orbit_size(group_name)` | orbit size (`"2T"`→4, `"2O"`→6, `"2I"`→12) |
| `check_coordination(n)` | `McKayCheck` (full diagnostic) |

### Bond orders (`chemistry.bond_orders`)

| Function / object | Returns |
|---|---|
| `huckel_bond_orders(adjacency, n_electrons=None)` | `HuckelBondOrderResult` — Coulson π bond orders |
| `cyclic_pi_bond_order(n_atoms, n_electrons=None)` | scalar (benzene → 2/3) |
| `smiles_pi_bond_orders(smiles)` | dict `{(i,j): P_ij}` from SMILES |
| `cc_bond_length_from_pi_order(p_pi)` | C–C length in Å (linear interp) |
| `PI_BOND_ORDER_BENZENE_CLASS` | `2/3 = 2/RANK_SO7` (≈ 0.6667) |

### Pericyclic — Woodward-Hoffmann (`chemistry.pericyclic`)

| Function / object | Returns |
|---|---|
| `selection_rule(n_electrons)` | `PericyclicSelectionRule` (parity, thermal/photo topology) |
| `reaction_selection_rule(label)` | rule by name (e.g. `"[4+2]_cycloaddition_Diels_Alder"`) |
| `electrocyclic_rotation_mode(n_electrons, thermal=True)` | `"conrotatory"` / `"disrotatory"` |
| `PERICYCLIC_TS_ELECTRON_COUNT` | TS electron-count catalog |

### Polyhedral / Wade-Mingos PSEPT (`chemistry.polyhedral`)

| Function / object | Returns |
|---|---|
| `wade_classification(n_vertices, n_seps)` | `WadeClass` (CLOSO / NIDO / ARACHNO / HYPHO) |
| `closo_borane_sep_count(n_vertices)` | SEP count (`n+1`) |
| `deltahedron_edge_count(n)` / `_face_count(n)` / `_euler_chi(n)` | `3n−6`, `2n−4`, `2` |
| `closo_borane_substrate_canonical(n_vertices)` | `ClosoCanonicalResult` — Spin(7) ladder labels |
| `CLOSO_POLYHEDRA` | dict (n = 5..12) of canonical closo polyhedra |

### Molecular knots (`chemistry.molecular_knots`)

| Function / object | Returns |
|---|---|
| `KNOT_REFERENCE` | dict of `MolecularKnotEntry` (`3_1`, `4_1`, `5_1`, `6_1`, …) |
| `carrier_class_of(crossing_number)` | substrate `n_q` carrier class |
| `substrate_predicted_tier(crossing_number)` | `AccessibilityTier` (SWEET_SPOT at n_q=5, GAP at n_q=6) |
| `accessibility_for_knot(alexander_briggs)` | full `MolecularKnotEntry` (synthesis + host-guest data) |
| `AccessibilityTier`, `HostGuestData`, `MolecularKnotEntry` | dataclasses / enum |

### Transition metal (`chemistry.transition_metal`)

| Function / object | Returns |
|---|---|
| `electron_count_class(n_electrons)` | `ElectronCountClass` (18e/16e/14e/12e/32e/other) |
| `ladder_k_for_count(n_electrons)` | Spin(7) ladder `k` (18 → 3) |
| `substrate_canonical_form(n_electrons)` | string identity (18 → `"N_EDGES_K7 − RANK_SO7 = 21 − 3"`) |
| `is_substrate_predicted_stable(n_electrons)` | bool |
| `transition_metal_entry(formula)` | `OrganometallicEntry` (e.g. `"Fe(CO)5"`) |
| `TRANSITION_METAL_REFERENCE` | reference dict (20 organometallics) |

### NMR ring-current sign (`chemistry.nmr`)

| Function / object | Returns |
|---|---|
| `nics_sign_from_hopf_parity(n_pi_electrons)` | `NICSSign` (DIATROPIC / PARATROPIC) |
| `nics_reference(name)` | `NICSReference` (incl. distinctive-integer hits) |
| `NICS_REFERENCE` | reference dict (14 molecules) |
| `NICSSign`, `StructurallyDistinctiveHit`, `NICSReference` | enums / dataclass |

> The NMR module exposes only the **sign** rule + two narrow distinctive-integer magnitude hits (coronene outer ≈ −18, benzene = −8 = −`DIM_OCTONION`); broad NICS *magnitude* prediction failed the rational-density audit and is deliberately not exposed.

### Benchmark (`chemistry.benchmark`)

| Function | Returns |
|---|---|
| `compare_to_dft(observable, n_molecules=1)` | `BenchmarkResult` |
| `timing_report(observable)` | `BenchmarkResult` |
| `run_full_benchmark_suite()` | list of 11 `BenchmarkResult` |
| `report_suite(results=None)` | pretty-printed comparison table |
| `benchmark_aromaticity`, `benchmark_c60_vibrational`, `benchmark_pah_resonance_energies`, `benchmark_wade_classification`, `benchmark_mckay_admissibility`, `benchmark_nics_sign_rule`, `benchmark_transition_metal_electron_count`, `benchmark_molecular_knot_accessibility`, … | per-observable benchmarks |

## Worked examples

### Benzene aromaticity + Tier-B resonance energy

```python
import nwt_substrate.chemistry as chem

a = chem.aromaticity_class("benzene")
print(a.classification, a.n_hopf_links, a.parity)   # aromatic 3 odd

# Calibrate the Hopf-pair scale on benzene's 36 kcal/mol, then predict
re = chem.aromatic_resonance_energy("naphthalene", calibration_kcal=36.0)
print(re.kcal_per_mol, re.n_hopf_pairs)             # 180.0 5

# With the default 12 kcal/mol scale (benchmark convention):
chem.aromatic_resonance_energy("coronene")          # 144.0 kcal/mol base + K_7 hub correction
```

### C_60 orbit + vibrational summary

```python
import nwt_substrate.chemistry as chem

orb = chem.fullerene_orbit_counts(60)
print(orb.n_vertices, orb.n_pentagons, orb.n_hexagons, orb.n_edges)
# 60 12 20 90   (pure |2I|/|D_n| group theory; is_2I_orbit=True)

v = chem.c60_vibrational_summary()
print(v.total_modes, v.ir_active_modes, v.raman_active_modes, v.silent_modes)
# 174 4 10 32   (matches empirical C_60 IR/Raman spectroscopy exactly)

for s in chem.c60_anion_magic_states():
    print(s.charge, s.closed_shell, s.empirical_match)
# 3  False  A_3 C_60 SC: K_3 (19 K), Rb_3 (28 K), Cs_3 (38 K)
# 6  True   A_6 C_60 hexafullerides (Bausch+ 1991)
# 9  False  UNTESTED — predicted SC under extreme reduction
# 12 True   UNTESTED — predicted under extreme electrochemistry
```

### SMILES round-trip

```python
import nwt_substrate.chemistry as chem

mol = chem.parse_smiles("c1ccccc1")              # benzene
print(len(mol.atoms), len(mol.bonds))            # 6 6

res = chem.smiles_to_aromaticity("c1ccccc1")
print(res.classification, res.total_hopf_count)  # aromatic 3

# Batch resonance energies (default 12 kcal/mol per Hopf pair)
chem.smiles_resonance_energy_batch(["c1ccccc1", "c1ccc2ccccc2c1"])
# → array([36., 60.])   benzene, naphthalene
```

### DFT-speedup benchmark

```python
from nwt_substrate.chemistry import benchmark as bm

r = bm.compare_to_dft("aromaticity", n_molecules=20)
print(r.observable, f"{r.speedup:.2e}", r.substrate_accuracy)
# aromaticity_classification 2.48e+09 '20/20 Hückel/Möbius matches in validation set'

print(bm.report_suite())   # full 11-observable substrate-vs-DFT table
#   molecular_knot_accessibility   speedup 4.10e+11  (max)
#   nics_sign_rule                 speedup 6.69e+09
#   c60_vibrational_modes          speedup 2.67e+09
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `60` | `\|A_5\|` = `2I` canonical vertex orbit (C_60 vertices) | binary icosahedral `2I` (E_8 McKay) |
| `12`, `20`, `90` | `\|2I\|/\|D_5\|`, `\|2I\|/\|D_3\|`, `60·6/2` (pentagons/hexagons/edges) | `2I` orbit / point group `D_n` |
| `{1, 2, 3, 4, 6, 12}` | McKay-admissible coordinations (`Z_n`/`BD_n`/`2T`/`2O`/`2I`) | McKay ADE correspondence |
| `2/3` | benzene π bond order `= 2/RANK_SO7` | `isa.RANK_SO7 = 3` |
| `{3, 5, 7, 9}` | Spin(7) rep-class ladder `(rank, dual Coxeter, dim_V, pos roots)` | `isa.RANK_SO7, H_V_SO7, N_VERTICES_K7, N_POS_ROOTS_SO7` |
| `18, 16, 14, 12` | stable TM counts `= N_EDGES_K7 − k`, `k ∈ {3,5,7,9}` | `isa.N_EDGES_K7 = 21` |
| `32` | f-block count `= K_7 triangles − RANK_SO7 = 35 − 3` | K_7 35 triangles |
| `n_q = 2` (Hopf link) | meson / aromatic-π-pair carrier | particle carrier table (`particles.md`) |
| `n_q = 5` (cinquefoil) | molecular-knot synthesis sweet-spot = proton-stable carrier | particle carrier table |
| `−8` | benzene NICS `= −DIM_OCTONION`; coronene-outer `≈ −18` (K_7 hub) | `isa.DIM_OCTONION = 8`; `K7_TR_M2_THRESHOLD = −24` hub indicator |
| `56.0` | K_7-toroidal RE correction (W_6 wheel signature) | `isa.K7_STABILIZATION_KCAL = 56.0` |

## Papers

The chemistry shim's group-theory inputs are shared with the particle sector — there is no standalone published chemistry paper in the source yet (an in-repo draft, `papers/paper23_chemistry.tex`, is in progress as of 2026-05). The carrier-knot `n_q` table and the Spin(7) representation-class ladder it reuses are documented in:

- **Paper 6** — carrier-knot `n_q ∈ {0..6}` table (shared with molecular-knot accessibility and aromaticity)
- **Paper 13** — Standard Model capstone / Spin(7) representation-class ladder (Wade boranes, transition-metal counts)

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`particles`](particles.md) — the carrier-knot `n_q` table (Hopf link = n_q 2, cinquefoil = n_q 5) that the aromaticity and molecular-knot modules reuse
- [`topology`](topology.md) — Hopf links, trefoils, and the knot-theory substrate
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`RANK_SO7`, `N_EDGES_K7`, `DIM_OCTONION`, Spin(7) ladder, `K7_STABILIZATION_KCAL`)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_chemistry`, `benchmark_aromatic_resonance_energies`, `benchmark_c60_vibrational_modes`, `benchmark_nmr_chemical_shifts`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
