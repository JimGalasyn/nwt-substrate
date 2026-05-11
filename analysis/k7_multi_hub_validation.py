"""K_7 toroidal ring-set detection — validation on larger PAHs.

Tests:
  - coronene (baseline, 1 K_7 set)
  - dicoronylene (2 coronenes σ-linked → 2 separate systems, each 1 K_7)
  - HBC attempts (hexa-peri-hexabenzocoronene)
  - ovalene-like (elongated, expected 0 K_7 hubs)

Key analytical result (derived before running):
  In HBC topology, ring-adjacency degrees are:
    center: 6 (all 6 petals)        ← only hub
    petals: 5 (center + 2 petals + 2 outer benzos)
    outer benzos: 2 (2 adjacent petals)
  So algorithmic prediction: 1 K_7 hub, not 2.

  This contradicts the earlier informal "HBC = 2 K_7 sets" claim
  which was a calibration guess to fit +128 kcal/mol empirical
  gap.  Per the current model, HBC's extra stabilization beyond
  +56 kcal/mol must come from a different mechanism.
"""
import nwt_substrate.chemistry as chem


CORONENE = "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"

# Dicoronylene: two coronenes σ-linked.  Ring labels 1-6 are free
# after the first coronene closes, so re-using them in the second
# coronene is legal SMILES.
DICORONYLENE = f"{CORONENE}-{CORONENE}"

# Hexa-peri-hexabenzocoronene (HBC, "superbenzene", C_42H_18, 13 rings).
#
# Building HBC SMILES from coronene + 6 peri-fused benzo groups.
# Each outer benzo shares 3 atoms (one full edge of one petal + 1
# vertex from each of 2 adjacent petals).
#
# This is a topology we attempt to construct; parser may not produce
# 13 perfect 6-rings (some SSSR decompositions yield different shapes).
HBC_ATTEMPT = (
    # Inner coronene-like part, but with extended rim:
    "c1cc2ccc3ccc4ccc5ccc6ccc7ccc8ccc9ccc%10ccc%11ccc%12ccc1"
    "c1c2c3c4c5c6c7c8c9c%10c%11c%121"
)

# Elongated 5x2 PAH — circumanthracene-like, no central hub expected.
# 10-ring elongated structure.
OVALENE_ATTEMPT = (
    "c1ccc2ccc3ccc4ccc5ccc6c1c2c3c4c5c6"
)


def describe(name: str, smiles: str):
    print(f"\n=== {name} ===")
    print(f"SMILES: {smiles[:80]}{'...' if len(smiles) > 80 else ''}")
    try:
        r = chem.smiles_to_aromaticity(smiles)
    except Exception as e:
        print(f"  PARSE ERROR: {e}")
        return None
    print(f"  atoms: {r.n_atoms}  aromatic atoms: {r.n_aromatic_atoms}")
    print(f"  ring systems: {len(r.aromatic_ring_systems)}")
    print(f"  classification: {r.classification}")
    print(f"  total π pairs: {r.total_pi_pairs}")
    total_k7 = sum(s.k7_ring_sets for s in r.aromatic_ring_systems)
    for i, s in enumerate(r.aromatic_ring_systems):
        ring_sizes = sorted([len(rr) for rr in s.rings])
        print(f"  system {i}: {len(s.rings)} rings  sizes={ring_sizes}  "
              f"π_pairs={s.n_pi_pairs}  Clar={s.clar_sextets}  "
              f"K_7={s.k7_ring_sets}  six_rings={s.n_six_rings}")
    print(f"  total K_7 sets: {total_k7}")
    re = chem.smiles_resonance_energy(smiles, calibration_kcal=12.0)
    print(f"  predicted RE: {re} kcal/mol "
          f"(base = {r.total_pi_pairs * 12} + K_7 = {total_k7 * 56})")
    return r


def main():
    print("K_7 multi-hub validation on larger PAHs\n" + "=" * 60)

    r_coronene = describe("coronene (baseline)", CORONENE)
    r_dicor = describe("dicoronylene (2 coronenes σ-linked)", DICORONYLENE)
    r_hbc = describe("HBC attempt (hexa-peri-hexabenzocoronene)", HBC_ATTEMPT)
    r_oval = describe("ovalene-like (elongated, expected 0 K_7)", OVALENE_ATTEMPT)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  coronene K_7 sets:     "
          f"{sum(s.k7_ring_sets for s in r_coronene.aromatic_ring_systems) if r_coronene else 'PARSE FAIL'}")
    if r_dicor:
        print(f"  dicoronylene K_7 sets: "
              f"{sum(s.k7_ring_sets for s in r_dicor.aromatic_ring_systems)} "
              f"(expected 2: one per system)")
    if r_hbc:
        print(f"  HBC K_7 sets:          "
              f"{sum(s.k7_ring_sets for s in r_hbc.aromatic_ring_systems)} "
              f"(by topology: only center has degree ≥ 6)")


if __name__ == "__main__":
    main()
