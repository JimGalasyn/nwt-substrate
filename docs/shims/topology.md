# nwt_substrate.topology

> The topological substrate primitives of NWT: torus-knot family / Seifert-genus / crossing-number / Hopf-charge classification (the **carrier knots** behind the particle catalog), the **K_7 Heffter genus-1 torus embedding** (`V=7, E=21, F=14`, `χ=0`), and the SU(2)_k colored-Jones / cabled-state machinery (quantum integers, torus channels, Rosso-Jones cabling-space dimension `n_q^q`, colored Jones + q-trace) at `DEFAULT_LEVEL = 5` (= SU(2)_5, the so(7)-forced Chern-Simons level with `k+2 = |V(K_7)| = 7`).

[← Back to index](../index.md) · Source: [`nwt_substrate/topology/`](../../nwt_substrate/topology/) · Papers: [6](https://zenodo.org/records/15376291) (carrier-knot mass / Hopf charge), [7](https://zenodo.org/communities/nwt) (`n_q^q` confinement = Rosso-Jones cabling), [21a](https://zenodo.org/communities/nwt) (cabled-state synthesis)

## Common questions

### What are torus knots in NWT?

Every NWT particle is hosted by a **carrier knot** — a torus knot `T(p, q)` whose crossing number `n_q ∈ [0, 6]` selects one of `7 = N_VERTICES_K7` carrier topologies (one per K_7 vertex). This module supplies the classical invariants of that knot:

| `n_q` | Carrier | Physics |
|---|---|---|
| 0 | unknot | leptons |
| 1 | unknot (extended) | charged-lepton tower |
| 2 | Hopf link | mesons (π, K, …) |
| 3 | trefoil `T(2,3)` | baryons (proton, neutron, …) |
| 4 | figure-eight | tetraquarks |
| 5 | cinquefoil `T(2,5)` | nucleons + pentaquarks |
| 6 | `6_1` knot | hexaquarks / dibaryons |

The carrier table here is the topological backbone of the [`particles`](particles.md) `(p, q, m, n_q)` mass formula. `knot_family(2, 3)` returns `"trefoil"`, `crossing_number(2, 3) = 3`, `seifert_genus(2, 3) = 1`.

### What is the K_7 Heffter embedding?

`K_7` (the complete graph on 7 vertices) admits Heffter's 1891 triangular embedding on the **torus**: `V = 7`, `E = 21`, `F = 14`, Euler characteristic `χ = V − E + F = 0`, hence **genus 1**, with **all 14 faces triangular**. The rotation system at vertex `i` has cyclic neighbor order `(i+1, i+3, i+2, i+6, i+4, i+5) mod 7`. This is the substrate's background geometry: the 21 edges carry the gravity-sector Wilson amplitude `α^(21/2)`, and the 14 faces are where `(p, q)` torus knots host vortex solitons. `is_genus_one_embedding(heffter_rotation())` returns `(True, {...})`.

### What is the colored Jones polynomial used for?

The classical invariants (genus, crossing number) are supplemented by the **quantum** SU(2)_k Witten-Reshetikhin-Turaev invariant, computed mechanically from the explicit `U_q(sl2)` R-matrix — no closed-form recall. A torus knot `T(strands, power)` is the closure of the `strands`-strand braid `(σ_1 … σ_{strands−1})^power`, each strand colored by spin `s`. The **Rosso-Jones cabling-space dimension** is `n_q^q` — which is *exactly* the NWT mass-formula confinement factor (Paper 7). The `colored_jones` Markov trace `Σ_J [2J+1] c_J` and its q-trace cross-check `colored_jones_qtrace` agree to `< 1e-7` and reproduce the ordinary Jones polynomial ratios at the root of unity `ζ_7`.

### What level k does it use and why?

`DEFAULT_LEVEL = 5`, i.e. **SU(2)_5**. This is the so(7)-forced Chern-Simons level: `k = h^∨(so7)`, `k+1 = h(so7)`, and crucially `k + 2 = 7 = |V(K_7)|`. So `q = exp(iπ/(k+2)) = exp(iπ/7)` is a 7th root of unity, and the quantum integer is `[n] = sin(nπ/7) / sin(π/7)`. The level-5 truncation wraps: `[5] = [2]`, `[6] = [1]`, `[7] = 0` — leaving 6 integrable anyons `w ∈ {0..5}`.

### What SU(2)_5 modular data emerges?

The 6 quantum dimensions are `[1, 1.802, 2.247, 2.247, 1.802, 1]`, giving total quantum dimension `D = √(Σ d²) = 4.3118` and chiral central charge `c = 15/7 = 2.142857` (from the Gauss sum). This is the SU(2)_5 modular tensor category that the K_7 substrate realizes; see `benchmark_modular_data`.

## Prediction table

| Structural fact | Substrate value | Meaning |
|---|---|---|
| K_7 Heffter embedding | `V=7, E=21, F=14, χ=0, genus=1` | Heffter 1891 triangular torus embedding |
| K_7 faces | 14 triangles | all `len(face) == 3` (vortex-soliton hosts) |
| K_7 edges | 21 | gravity Wilson amplitude `α^(21/2)` |
| trefoil `T(2,3)` | `crossing=3, genus=1, family="trefoil"` | baryon carrier (`n_q=3`) |
| cinquefoil `T(2,5)` | `crossing=5, genus=2` | nucleon/pentaquark carrier (`n_q=5`) |
| Hopf link `T(2,2)` | `family="Hopf"` | meson carrier (`n_q=2`) |
| Hopf charge `Q_H` | `p·m` (e.g. `hopf_charge(2,1)=2`) | L3 Hopf charge (Paper 6/7) |
| quantum integers `[1..6]` at k=5 | `[1, 1.8019, 2.2470, 2.2470, 1.8019, 1]` | `[n]=sin(nπ/7)/sin(π/7)`; `[7]=0` |
| `|colored_jones|` trefoil `T(2,3)` s=½ | `3.17771` | unreduced Markov trace at ζ_7 |
| `|colored_jones|` Hopf `T(2,2)` s=½ | `2.24698` (`= [3]`) | meson-carrier quantum invariant |
| `|colored_jones|` cinquefoil `T(2,5)` s=½ | `2.65498` | nucleon-carrier quantum invariant |
| `|colored_jones|` `T(3,5)` s=½ | `3.17771` | lands on trefoil value (famous ζ_7 coincidence) |
| cabling-space dim `n_q^q` | `T(3,4)→81`, `T(2,3)→8`, `T(5,3)→125` | Rosso-Jones rep dimension = Paper 7 `n_q^q` |
| SU(2)_5 total quantum dim | `D = 4.3118` | `√(Σ d²)`, `D² = 18.59` |
| SU(2)_5 chiral central charge | `c = 15/7 = 2.142857` | from Gauss sum (`benchmark_modular_data`) |

Each row asserted in [`nwt_substrate/tests/test_torus_knots.py`](../../nwt_substrate/tests/test_torus_knots.py), [`test_K7.py`](../../nwt_substrate/tests/test_K7.py), and [`test_colored_jones.py`](../../nwt_substrate/tests/test_colored_jones.py).

## Quick start

```python
import nwt_substrate.topology as t

# --- Torus-knot classification (the carrier knots) ---
t.knot_family(2, 3)        # → "trefoil"   (baryon carrier, n_q=3)
t.crossing_number(2, 3)    # → 3
t.seifert_genus(2, 3)      # → 1
t.is_torus_knot(2, 3)      # → True   (gcd=1, single component)
t.is_torus_knot(2, 4)      # → False  (gcd=2, multi-component link)
t.hopf_charge(2, 1)        # → 2      (Q_H = p·m)
t.knot_family(2, 2)        # → "Hopf"
t.knot_family(2, 5)        # → "cinquefoil"

# --- K_7 Heffter genus-1 torus embedding ---
rot = t.heffter_rotation()           # {0: [1,3,2,6,4,5], 1: [...], ...}
passes, info = t.is_genus_one_embedding(rot)
passes                               # → True
info["V"], info["E"], info["F"]      # → (7, 21, 14)
info["chi"], info["genus"]           # → (0, 1)
len(t.trace_K7_faces(rot))           # → 14   (all triangular)

# --- Colored Jones / SU(2)_5 (DEFAULT_LEVEL = 5) ---
t.DEFAULT_LEVEL                      # → 5
t.quantum_integer(2).real            # → 1.801938  (= 2 cos(π/7))
t.quantum_integer(7)                 # → ~0+0j     (k=5 truncation: [7]=0)
abs(t.colored_jones(0.5, 2, 3))      # → 3.17771   (trefoil at ζ_7)
abs(t.colored_jones(0.5, 2, 2))      # → 2.24698   (Hopf link = [3])
t.cabling_space_dimension(3, 4)      # → 81        (Rosso-Jones n_q^q)
```

## API by topic

### Torus knots ([`torus_knots.py`](../../nwt_substrate/topology/torus_knots.py))

| Function | Returns |
|---|---|
| `knot_family(p, q)` | Named carrier family: `"unknot"`, `"Hopf"`, `"trefoil"`, `"cinquefoil"`, `"heptafoil"`, `"septafoil"`, or generic `"T(p,q)"` |
| `seifert_genus(p, q)` | Seifert genus `g = (p−1)(q−1)/2` for `gcd(p,q)=1`; `0` for multi-component links |
| `crossing_number(p, q)` | `min(p(q−1), q(p−1))` for `p,q ≥ 2`; `0` for the unknot |
| `is_torus_knot(p, q)` | `True` iff `gcd(p, q) = 1` (single-component knot, not a link) |
| `hopf_charge(p, m)` | L3 Hopf charge `Q_H = p·m` |

### K_7 embedding ([`K7.py`](../../nwt_substrate/topology/K7.py))

| Function | Returns |
|---|---|
| `heffter_rotation()` | Dict `{v: [n1..n6]}` — Heffter's cyclic rotation system, `(v+b) mod 7` for `b ∈ [1,3,2,6,4,5]` |
| `trace_K7_faces(rotation)` | List of faces (tuples of vertices) traced from a rotation system |
| `is_genus_one_embedding(rotation)` | `(passes, info)` — checks `χ=0` (torus) + all triangular faces; `info` carries `V, E, F, chi, face_sizes, all_triangles, genus` |

### Colored Jones / SU(2)_k ([`colored_jones.py`](../../nwt_substrate/topology/colored_jones.py))

| Function | Returns |
|---|---|
| `DEFAULT_LEVEL` | `5` — the so(7)-forced Chern-Simons level (SU(2)_5, `k+2=7`) |
| `quantum_integer(n, level=5)` | `[n] = (q^n − q^−n)/(q − q^−1) = sin(nπ/(k+2))/sin(π/(k+2))` |
| `torus_channels(spin, strands)` | Dict `{J: multiplicity}` of spins `J` in `V_s^{⊗strands}` |
| `cabling_space_dimension(n_q, q)` | `n_q^q` — Rosso-Jones cabling-space dim = Paper 7 confinement factor |
| `cabled_state(spin, strands, power, level=5, kmax=None)` | Wilson-line-basis state `{J: c_J(power)}` for the closure of `(σ_1…σ_{strands−1})^power` |
| `colored_jones(spin, strands, power, level=5)` | Unreduced Markov trace `Σ_J [2J+1]·c_J(power)` of `T(strands, power)` |
| `colored_jones_qtrace(spin, strands, power, level=5)` | Same invariant via direct quantum trace `Tr(ρ(β)^power · μ^{⊗strands})` — cross-check of `colored_jones` |

## Worked examples

### Classify the proton's trefoil carrier

```python
import nwt_substrate.topology as t

# The proton (n_q=3) is carried by the trefoil knot T(2,3).
p, q = 2, 3
print(t.knot_family(p, q))     # trefoil
print(t.crossing_number(p, q)) # 3
print(t.seifert_genus(p, q))   # 1
print(t.is_torus_knot(p, q))   # True   (gcd(2,3)=1)
print(t.hopf_charge(p, 1))     # 2      (Q_H = p·m)
```

### Verify the K_7 Heffter V/E/F genus-1 embedding

```python
import nwt_substrate.topology as t

rot = t.heffter_rotation()
faces = t.trace_K7_faces(rot)
passes, info = t.is_genus_one_embedding(rot)

print(passes)                                   # True
print(info["V"], info["E"], info["F"])          # 7 21 14
print(info["chi"], info["genus"])               # 0 1
print(len(faces), all(len(f) == 3 for f in faces))  # 14 True
print(faces[:3])  # [(0, 1, 3), (0, 3, 2), (0, 2, 6)]
```

### Evaluate colored Jones of the trefoil at k=5

```python
import nwt_substrate.topology as t

# Trefoil T(2,3) = T(strands=2, power=3), colored by spin 1/2
cj  = t.colored_jones(0.5, 2, 3)        # state-sum Markov trace
cjq = t.colored_jones_qtrace(0.5, 2, 3) # direct quantum-trace cross-check
print(abs(cj))                          # 3.177709
print(abs(cj - cjq) < 1e-7)             # True  (the two routes agree)

# Hopf link T(2,2) — the meson carrier — equals [3] at k=5
print(abs(t.colored_jones(0.5, 2, 2)))  # 2.246980  (= quantum_integer(3))

# n_q^q confinement factor = cabling-space dimension
print(t.cabling_space_dimension(2, 3))  # 8
print(t.cabling_space_dimension(3, 4))  # 81
```

## Substrate constants used here

| Symbol | Substrate identity | Source |
|---|---|---|
| `7` | `\|V(K_7)\| = N_VERTICES_K7`; also `= k+2` (SU(2)_5 order) and `N_CARRIER_TYPES` (one carrier per vertex) | `isa.N_VERTICES_K7` |
| `21` | `\|E(K_7)\| = C(7,2)` (gravity Wilson amplitude `α^(21/2)`) | `isa.N_EDGES_K7` |
| `14` | `F(K_7)` Heffter triangular faces (`V−E+F = 0` ⟹ genus 1) | derived (`is_genus_one_embedding`) |
| `n_q^q` | Rosso-Jones cabling-space dimension = Paper 7 mass-formula confinement factor | `cabling_space_dimension` |
| `k = 5` | SU(2)_5 Chern-Simons level `= h^∨(so7)` (so(7)-forced) | `topology.DEFAULT_LEVEL` |
| `k + 2 = 7` | root-of-unity order `q = exp(iπ/7) = ζ_7`; `[7] = 0` | `isa.N_VERTICES_K7` |
| `n_q ∈ [0, 6]` | crossing number = carrier index; `MAX_CROSSING_NUMBER = N_VERTICES_K7 − 1` | `isa.MAX_CROSSING_NUMBER` |

## Papers

- **Paper 6** — Carrier-knot mass formula; L3 Hopf charge `Q_H = p·m`, `(p, q, m, n_q)` torus-knot indexing ([Zenodo 15376291](https://zenodo.org/records/15376291)).
- **Paper 7** — The `n_q^q` confinement factor identified as the Rosso-Jones cabling-space dimension (the dimension of the rep the Rosso-Jones operator acts on).
- **Paper 21a** — Cabled-state synthesis: colored-Jones Markov trace reproducing Jones-polynomial ratios at the root of unity, including multiplicity-`>1` cablings.

[Full series on Zenodo](https://zenodo.org/communities/nwt).

## See also

- [`particles`](particles.md) — the `(p, q, m, n_q)` carrier-knot mass spectrum that this module's classical invariants underpin
- [`chemistry`](chemistry.md) — aromaticity / NICS / C_60
- [`heron`](heron.md) — IBM Heron adapter; the Hopf-link colored-Jones S-matrix reconstruction runs on hardware
- [`isa`](../../nwt_substrate/isa/README.md) — substrate constants (`N_VERTICES_K7 = 7`, `N_EDGES_K7 = 21`, `N_CARRIER_TYPES`, `MAX_CROSSING_NUMBER`)
- [`benchmarks`](../../nwt_substrate/benchmarks/README.md) — `benchmark_k7_face_structure`, `benchmark_modular_data`
- [`docs/FAQ.md`](../FAQ.md) — atomic Q&A summaries
