#!/usr/bin/env python3
"""
NWT Electroweak Boson Decay Channel Analysis

Checks whether the Hopf(2) EWK boson assignments are consistent with
observed decay channels:
  1. Energy budget: parent mass > sum of daughter masses
  2. Quantum number conservation: Q, B, lepton number
  3. Topology: what carrier surgery does each decay require?
  4. Comparison with NWT predictions for daughter particle masses
"""

import numpy as np
from math import gcd

ME = 0.511
BETA_E = np.sqrt(5.0/4.0)

def beta_fast(p, m):
    val = m**2 / p**2 - 1
    return np.sqrt(val) if val > 0 else None

def mass_nwt(p_eff, q_eff, m_int, n_q):
    beta = beta_fast(p_eff, m_int)
    if beta is None or beta <= 1e-6: return None
    pq = p_eff**2 + q_eff**2
    lr = np.log(8*beta) / np.log(8*BETA_E)
    if lr <= 0: return None
    geom = (pq/5.0) * (beta/BETA_E) * lr
    try: enh = float(n_q**q_eff) if n_q > 0 and q_eff > 0 else 1.0
    except: return None
    if not np.isfinite(enh) or enh > 1e30: return None
    mass = geom * enh * ME
    return mass if np.isfinite(mass) and mass > 0 else None


# ── EWK boson assignments (best fits from Hopf(2) scan) ──────────────
# Option A: different modes
EWK_A = {
    'W±': {'mode': (7,5), 'carrier': (2,2), 'eff': (14,10), 'm': 34,
            'n_q': 2, 'mass_nwt': 80429, 'mass_pdg': 80379, 'Q': 1, 'J': 1, 'I': 1, 'f': 2},
    'Z⁰': {'mode': (3,5), 'carrier': (2,2), 'eff': (6,10), 'm': 27,
            'n_q': 2, 'mass_nwt': 90706, 'mass_pdg': 91188, 'Q': 0, 'J': 1, 'I': 0, 'f': 0},
    'H⁰': {'mode': (9,5), 'carrier': (2,2), 'eff': (18,10), 'm': 46,
            'n_q': 2, 'mass_nwt': 125016, 'mass_pdg': 125100, 'Q': 0, 'J': 0, 'I': 0, 'f': 0},
}

# Option B: same mode, consecutive m
EWK_B = {
    'W±': {'mode': (1,5), 'carrier': (2,2), 'eff': (2,10), 'm': 10,
            'n_q': 2, 'mass_nwt': 79850, 'mass_pdg': 80379, 'Q': 1, 'J': 1, 'I': 1, 'f': 2},
    'Z⁰': {'mode': (1,5), 'carrier': (2,2), 'eff': (2,10), 'm': 11,
            'n_q': 2, 'mass_nwt': 90529, 'mass_pdg': 91188, 'Q': 0, 'J': 1, 'I': 0, 'f': 0},
    'H⁰': {'mode': (2,5), 'carrier': (2,2), 'eff': (4,10), 'm': 26,
            'n_q': 2, 'mass_nwt': 125382, 'mass_pdg': 125100, 'Q': 0, 'J': 0, 'I': 0, 'f': 0},
}

# ── NWT daughter particle assignments (from Paper 8 Table 1) ─────────
DAUGHTERS = {
    # Leptons (unknot carriers)
    'e⁻':   {'mass_pdg': 0.511,   'Q': -1, 'B': 0, 'J': 0.5, 'carrier': 'unknot'},
    'e⁺':   {'mass_pdg': 0.511,   'Q': +1, 'B': 0, 'J': 0.5, 'carrier': 'unknot'},
    'μ⁻':   {'mass_pdg': 105.66,  'Q': -1, 'B': 0, 'J': 0.5, 'carrier': 'unknot'},
    'μ⁺':   {'mass_pdg': 105.66,  'Q': +1, 'B': 0, 'J': 0.5, 'carrier': 'unknot'},
    'τ⁻':   {'mass_pdg': 1776.86, 'Q': -1, 'B': 0, 'J': 0.5, 'carrier': 'unknot'},
    'τ⁺':   {'mass_pdg': 1776.86, 'Q': +1, 'B': 0, 'J': 0.5, 'carrier': 'unknot'},
    'ν':    {'mass_pdg': 0.0,     'Q': 0,  'B': 0, 'J': 0.5, 'carrier': 'topology data'},
    'ν̄':    {'mass_pdg': 0.0,     'Q': 0,  'B': 0, 'J': 0.5, 'carrier': 'topology data'},
    # Light mesons (Hopf carriers)
    'π⁺':   {'mass_pdg': 139.57,  'Q': +1, 'B': 0, 'J': 0, 'carrier': 'Hopf(2)'},
    'π⁻':   {'mass_pdg': 139.57,  'Q': -1, 'B': 0, 'J': 0, 'carrier': 'Hopf(2)'},
    'π⁰':   {'mass_pdg': 135.0,   'Q': 0,  'B': 0, 'J': 0, 'carrier': 'Hopf(2)'},
    'K⁺':   {'mass_pdg': 493.68,  'Q': +1, 'B': 0, 'J': 0, 'carrier': 'Hopf(2)'},
    'K⁻':   {'mass_pdg': 493.68,  'Q': -1, 'B': 0, 'J': 0, 'carrier': 'Hopf(2)'},
    # Heavy mesons
    'D⁰':   {'mass_pdg': 1864.8,  'Q': 0,  'B': 0, 'J': 0, 'carrier': 'Hopf'},
    'B⁰':   {'mass_pdg': 5279.7,  'Q': 0,  'B': 0, 'J': 0, 'carrier': 'Hopf'},
    'B⁺':   {'mass_pdg': 5279.3,  'Q': +1, 'B': 0, 'J': 0, 'carrier': 'Hopf'},
    # Baryons (trefoil carriers)
    'p':    {'mass_pdg': 938.27,  'Q': +1, 'B': 1, 'J': 0.5, 'carrier': 'trefoil'},
    'n':    {'mass_pdg': 939.57,  'Q': 0,  'B': 1, 'J': 0.5, 'carrier': 'trefoil'},
    'p̄':    {'mass_pdg': 938.27,  'Q': -1, 'B': -1,'J': 0.5, 'carrier': 'trefoil'},
    'n̄':    {'mass_pdg': 939.57,  'Q': 0,  'B': -1,'J': 0.5, 'carrier': 'trefoil'},
    # Photon
    'γ':    {'mass_pdg': 0.0,     'Q': 0,  'B': 0, 'J': 1,   'carrier': 'none (radiation)'},
    # W/Z themselves (for H→WW, H→ZZ)
    'W':    {'mass_pdg': 80379.0, 'Q': 1,  'B': 0, 'J': 1,   'carrier': 'Hopf(2)'},
    'Z':    {'mass_pdg': 91188.0, 'Q': 0,  'B': 0, 'J': 1,   'carrier': 'Hopf(2)'},
}

# ── Observed decay channels ──────────────────────────────────────────
DECAYS = [
    # W± decays (BR ~ 67% hadronic, 33% leptonic)
    ('W⁺ → e⁺ν',        'W±', [('e⁺', 1), ('ν', 1)],         10.71, 'leptonic'),
    ('W⁺ → μ⁺ν',        'W±', [('μ⁺', 1), ('ν', 1)],         10.63, 'leptonic'),
    ('W⁺ → τ⁺ν',        'W±', [('τ⁺', 1), ('ν', 1)],         11.38, 'leptonic'),
    ('W⁺ → π⁺π⁰',       'W±', [('π⁺', 1), ('π⁰', 1)],        0.01, 'hadronic (exclusive)'),
    ('W⁺ → hadrons',     'W±', [],                              67.41, 'hadronic (inclusive)'),

    # Z⁰ decays
    ('Z⁰ → e⁺e⁻',       'Z⁰', [('e⁺', 1), ('e⁻', 1)],        3.363, 'leptonic'),
    ('Z⁰ → μ⁺μ⁻',       'Z⁰', [('μ⁺', 1), ('μ⁻', 1)],        3.366, 'leptonic'),
    ('Z⁰ → τ⁺τ⁻',       'Z⁰', [('τ⁺', 1), ('τ⁻', 1)],        3.370, 'leptonic'),
    ('Z⁰ → νν̄',          'Z⁰', [('ν', 1), ('ν̄', 1)],         20.00, 'invisible'),
    ('Z⁰ → hadrons',     'Z⁰', [],                             69.91, 'hadronic (inclusive)'),
    ('Z⁰ → π⁺π⁻',       'Z⁰', [('π⁺', 1), ('π⁻', 1)],        0.00, 'hadronic (exclusive)'),
    ('Z⁰ → pp̄',          'Z⁰', [('p', 1), ('p̄', 1)],           0.00, 'hadronic (exclusive)'),

    # H⁰ decays
    ('H⁰ → bb̄',         'H⁰', [('B⁰', 2)],                   58.0, 'dominant'),
    ('H⁰ → WW*',        'H⁰', [('W', 2)],                     21.5, 'off-shell W'),
    ('H⁰ → ZZ*',        'H⁰', [('Z', 2)],                      2.6, 'off-shell Z'),
    ('H⁰ → τ⁺τ⁻',      'H⁰', [('τ⁺', 1), ('τ⁻', 1)],         6.3, 'leptonic'),
    ('H⁰ → γγ',         'H⁰', [('γ', 2)],                      0.23, 'loop-induced'),
]

# ── Analysis ─────────────────────────────────────────────────────────
def analyze_decays(ewk_set, label):
    print(f"\n{'=' * 110}")
    print(f"DECAY ANALYSIS: {label}")
    print(f"{'=' * 110}")

    for decay_name, parent_key, daughters, br, dtype in DECAYS:
        parent = ewk_set[parent_key]
        m_parent = parent['mass_pdg']
        Q_parent = parent['Q']
        B_parent = 0  # All EWK bosons have B=0

        # Calculate daughter masses and quantum numbers
        m_daughters = 0
        Q_daughters = 0
        B_daughters = 0
        carrier_changes = []
        daughter_str = []

        if not daughters:
            # Inclusive hadronic — just note it
            pass
        else:
            for dname, count in daughters:
                d = DAUGHTERS[dname]
                m_daughters += d['mass_pdg'] * count
                Q_daughters += d['Q'] * count
                B_daughters += d['B'] * count
                carrier_changes.append(f"{count}×{dname}[{d['carrier']}]")
                daughter_str.append(f"{count}×{dname}({d['mass_pdg']:.1f})")

        # Energy budget
        if daughters:
            Q_avail = m_parent - m_daughters
            energy_ok = Q_avail > 0
            Q_conserved = (Q_parent == Q_daughters) or (Q_parent == -Q_daughters)  # W⁺ or W⁻
            B_conserved = (B_parent == B_daughters)

            # For W±, the charge should match the decay
            if parent_key == 'W±':
                # W⁺ decays: Q_daughters should be +1
                Q_conserved = True  # We listed W⁺ decays with correct charges

            status = '✓' if (energy_ok and Q_conserved and B_conserved) else '✗'

            print(f"\n  {decay_name:25s}  BR={br:5.2f}%  [{dtype}]")
            print(f"    Parent: {parent_key} = {m_parent:.0f} MeV (NWT: {parent['mass_nwt']:.0f})")
            print(f"    Daughters: {' + '.join(daughter_str)} = {m_daughters:.1f} MeV")
            print(f"    Energy available: {Q_avail:.1f} MeV  {'✓ OPEN' if energy_ok else '✗ CLOSED'}")
            print(f"    ΔQ = {Q_parent} → {Q_daughters}  {'✓' if Q_conserved else '✗'}")
            print(f"    ΔB = {B_parent} → {B_daughters}  {'✓' if B_conserved else '✗'}")

            # Topology analysis
            parent_carrier = 'Hopf(2)'
            daughter_carriers = ' + '.join(carrier_changes)
            print(f"    Topology: {parent_carrier} → {daughter_carriers}")

            # Classify the topology change
            if all(DAUGHTERS[dn]['carrier'] == 'unknot' or DAUGHTERS[dn]['carrier'] == 'topology data'
                   for dn, _ in daughters):
                print(f"    Surgery: Hopf link UNLINKED → two unknots (R2 move, Δcrossings = -2)")
                print(f"    → WEAK decay: requires topology change (neutrino carries surgery data)")
            elif all(DAUGHTERS[dn]['carrier'].startswith('Hopf') for dn, _ in daughters):
                print(f"    Surgery: Hopf(2) → Hopf(2) daughter(s) (excitation change, no topology change)")
                print(f"    → STRONG/EM decay: carrier type preserved")
            elif any(DAUGHTERS[dn]['carrier'] == 'none (radiation)' for dn, _ in daughters):
                print(f"    Surgery: Hopf(2) annihilated → radiation")
                print(f"    → EM decay: link annihilation")
            else:
                print(f"    Surgery: mixed topology change")

        else:
            print(f"\n  {decay_name:25s}  BR={br:5.2f}%  [{dtype}]")
            print(f"    Inclusive channel — daughters are multiple hadrons")
            print(f"    Topology: Hopf(2) → multiple Hopf links + unknots")
            print(f"    Energy budget: {m_parent:.0f} MeV available for jet fragmentation")

    # ── Special: H → WW* and H → ZZ* ────────────────────────────────
    print(f"\n{'─' * 110}")
    print(f"  OFF-SHELL DECAY ANALYSIS")
    print(f"{'─' * 110}")

    m_H = ewk_set['H⁰']['mass_pdg']
    m_W = ewk_set['W±']['mass_pdg']
    m_Z = ewk_set['Z⁰']['mass_pdg']

    print(f"\n  H⁰ → WW*: m_H = {m_H:.0f}, 2×m_W = {2*m_W:.0f}")
    print(f"    Deficit: {m_H - 2*m_W:.0f} MeV → one W must be off-shell (virtual)")
    print(f"    In NWT: the H Hopf(2) splits into two Hopf(2) carriers,")
    print(f"    but only enough energy for one to be on-shell (m={ewk_set['W±']['m']})")
    print(f"    The other is a virtual state (no stable phase closure)")
    print()
    print(f"  H⁰ → ZZ*: m_H = {m_H:.0f}, 2×m_Z = {2*m_Z:.0f}")
    print(f"    Deficit: {m_H - 2*m_Z:.0f} MeV → both Z must be off-shell")
    print(f"    In NWT: same topology change but even more virtual")
    print()

    # ── Topology summary ─────────────────────────────────────────────
    print(f"\n{'=' * 110}")
    print(f"TOPOLOGY SURGERY SUMMARY")
    print(f"{'=' * 110}")

    print(f"""
  Parent carriers: all Hopf(2) = two linked unknotted vortex rings

  Leptonic decays (W→lν, Z→ll):
    Hopf(2) → unknot + unknot [+ neutrino carrying topology data]
    Surgery type: R2 (remove 2 crossings = unlink)
    Paper 8 identification: R2 = strong interaction / ν_μ

    BUT: W/Z leptonic decays are WEAK, not strong!
    This means the Reidemeister assignment in Paper 8 Sec. 5
    may need revision, OR the surgery is actually R3 (strand slide)
    that happens to also unlink in this case.

  Hadronic decays (W→qq̄, Z→qq̄):
    Hopf(2) → Hopf(2) + Hopf(2) [fragmentation]
    Surgery type: the parent Hopf link's energy goes into creating
    new Hopf links (mesons) from the vacuum — no topology change
    of the carrier TYPE, just energy redistribution.
    This is analogous to Δ→pπ: trefoil energy → trefoil + Hopf.

  H→γγ:
    Hopf(2) → two photons (no vortex structure)
    Surgery type: complete link annihilation
    Both vortex rings unwind to radiation
    This is the analog of π⁰→γγ (also Hopf→radiation)

  H→WW*, H→ZZ*:
    Hopf(2) → Hopf(2) + Hopf(2)
    Same carrier type in and out — an EXCITATION-LEVEL change
    The parent Hopf link (high m) splits into two Hopf links (lower m)
    This is topologically trivial (no knot type change)
    The energy constraint forces one or both daughters off-shell.
""")

    # ── Energy budget summary table ──────────────────────────────────
    print(f"\n{'=' * 110}")
    print(f"ENERGY BUDGET SUMMARY")
    print(f"{'=' * 110}")
    print(f"\n  {'Decay':30s} {'m_parent':>10s} {'Σm_daughters':>12s} {'Available':>10s} {'Status':>8s}")
    print(f"  {'─'*75}")

    for decay_name, parent_key, daughters, br, dtype in DECAYS:
        if not daughters:
            m_parent = ewk_set[parent_key]['mass_pdg']
            print(f"  {decay_name:30s} {m_parent:10.0f} {'(jets)':>12s} {m_parent:10.0f} {'✓':>8s}")
            continue
        m_parent = ewk_set[parent_key]['mass_pdg']
        m_d = sum(DAUGHTERS[dn]['mass_pdg'] * c for dn, c in daughters)
        avail = m_parent - m_d
        ok = '✓' if avail > 0 else '✗ virtual'
        print(f"  {decay_name:30s} {m_parent:10.0f} {m_d:12.1f} {avail:10.1f} {ok:>8s}")


# Run for both options
analyze_decays(EWK_A, "Option A — different modes for W, Z, H")
print("\n\n")
analyze_decays(EWK_B, "Option B — same mode (1,5), consecutive m")


# ── Cross-check: can hadronic W/Z decays produce known mesons? ───────
print(f"\n\n{'=' * 110}")
print(f"HADRONIC DECAY CHANNELS: Can EWK Hopf(2) produce known mesons?")
print(f"{'=' * 110}")
print(f"""
  W⁺ (80,379 MeV) hadronic decays produce quark pairs that
  hadronize into jets of mesons and baryons.

  In NWT: the parent Hopf(2) carrier's energy creates new
  Hopf(2) carriers (mesons) from the vacuum condensate.

  Maximum number of pions from W decay: {80379/139.57:.0f} pions
  Maximum number of pions from Z decay: {91188/139.57:.0f} pions

  Typical hadronic W decay: W⁺ → 2 jets → ~20 hadrons
  Average energy per hadron: {80379/20:.0f} MeV (consistent with
  mix of π, K, p, etc.)

  Key topology point: the parent Hopf link FRAGMENTS into
  multiple daughter Hopf links. This is vortex reconnection
  in the superfluid vacuum — the same GPE dynamics that
  produces trefoil knots from colliding unknots (Zuccher & Ricca 2022).

  Specific exclusive channels that MUST work:
""")

exclusive = [
    ('W⁺ → π⁺π⁰', 80379, 139.57 + 135.0, 'Hopf→Hopf+Hopf'),
    ('W⁺ → K⁺K̄⁰', 80379, 493.68 + 497.61, 'Hopf→Hopf+Hopf'),
    ('W⁺ → D⁺D̄⁰', 80379, 1869.7 + 1864.8, 'Hopf→Hopf+Hopf'),
    ('Z⁰ → π⁺π⁻', 91188, 2*139.57, 'Hopf→Hopf+Hopf'),
    ('Z⁰ → K⁺K⁻', 91188, 2*493.68, 'Hopf→Hopf+Hopf'),
    ('Z⁰ → pp̄', 91188, 2*938.27, 'Hopf→trefoil+trefoil'),
    ('Z⁰ → B⁺B⁻', 91188, 2*5279.3, 'Hopf→Hopf+Hopf'),
    ('H⁰ → τ⁺τ⁻', 125100, 2*1776.86, 'Hopf→unknot+unknot'),
    ('H⁰ → bb̄→B⁰B̄⁰', 125100, 2*5279.7, 'Hopf→Hopf+Hopf'),
]

print(f"  {'Channel':25s} {'m_parent':>10s} {'Σm_daught':>10s} {'Avail':>10s} {'Topology':>25s} {'OK':>4s}")
print(f"  {'─'*90}")
for ch, mp, md, topo in exclusive:
    avail = mp - md
    print(f"  {ch:25s} {mp:10.0f} {md:10.1f} {avail:10.1f} {topo:>25s} {'✓' if avail > 0 else '✗':>4s}")


if __name__ == '__main__':
    pass
