# Benchmark-suite triage — 2026-07-12

**Full failure-mode classification of all 38 `benchmark_*` functions in
`nwt_substrate/benchmarks/compute_speed.py`, under the standard set by the
pre-registration kill surface (`benchmarks/surface.py`) and the
derivation-audit gauntlet.**

Provenance of this document: classification performed by an audit pass over
the benchmark file and every underlying module (session 2026-07-12, the same
session that built `surface.py` and `neighbouring_value.py`).  Nine rows were
independently spot-verified line-by-line (`alpha_derivation`, `mass_spectrum`,
`gravitational_constant`, `higgs_vev`, `electron_anomaly`,
`c60_vibrational_modes`, `qcd_constants`, `nmr_chemical_shifts`,
`composite_particles`) — all matched.  The remaining rows carry the auditing
agent's file:line evidence and are Auditor-checkable at this commit's SHA.

## Classes

- **TAUTOLOGY** — the accuracy line compares the prediction against itself or
  the same formula/table that generated it.  Cannot fail; carries zero
  information.
- **STANDARD-RELABELED** — the result is textbook physics/chemistry (group
  theory, 3N−6 counts, shell/Hückel combinatorics, LO QED formulas) that
  holds with no K₇/substrate structure at all (gauntlet mode 3: the structure
  is inert).
- **MENU-FITTED** — a closed form assembled from the structural-integer menu
  and matched to a measured target; discounted by the measured look-elsewhere
  volume (`neighbouring_value.py`: ~83% of random targets fit inside G's bar
  on the full menu).
- **MEASURED-INPUT** — the computation takes a measured value as input
  (ALPHA_QED defaults, PDG masses, hardcoded α_s) and its "prediction"
  therefore smuggles the target or a sibling measured quantity.
- **METROLOGY** (secondary only) — an impressive %/ppm accuracy quoted against
  a target measured far more precisely; dead-as-exact in experimental σ.
- **QUALITATIVE/SCORE** (secondary only) — N/N or sign-rule scores.
- **GENUINE-FORWARD** — no measured witness exists yet; untestable today but
  uncontaminated as a test.

## Classification (source order)

| benchmark | primary class | secondary | evidence (file:line, one clause) |
|---|---|---|---|
| benchmark_alpha_derivation | MENU-FITTED | METROLOGY | compute_speed.py:56,63 — 1/α=25π√3+1 (isa/constants.py:100) scored "7.6 ppm vs CODATA", i.e. dead at ~50,000σ in experimental σ |
| benchmark_mass_spectrum | MENU-FITTED | MEASURED-INPUT | particles/compendium.py:6-8 — per-particle integers (p,q,m,n_q) chosen per particle, nucleons post-hoc corrected (1,4,\*,3)→(1,3,\*,5) to drop residual 1.06%→0.76%; m_e anchor (mass.py:22) |
| benchmark_modular_data | STANDARD-RELABELED | TAUTOLOGY | compute_speed.py:120-134 — textbook SU(2)_k quantum-group S/T/spin formulas; accuracy "exact (closed-form)" has no external witness, only internal consistency c=15/7 |
| benchmark_ckm_cabibbo | MENU-FITTED | METROLOGY | compute_speed.py:153 — λ²=7·ALPHA_SUBSTRATE, integer 7 from menu, matched to PDG λ=0.2253 (actual 0.22601 = 0.32%, not "~0.1%") |
| benchmark_k7_face_structure | STANDARD-RELABELED | QUALITATIVE/SCORE | topology/K7.py:2-11 — re-verifies Heffter 1891 toroidal embedding of K_7, pure known combinatorics, no physics claim tested |
| benchmark_wimp_tower | GENUINE-FORWARD | MENU-FITTED, MEASURED-INPUT | compute_speed.py:211 — masses = α^(N_e/2)·M_PL_GEV with per-species free integer N_e; 98 GeV/38 keV rungs unmeasured; M_PL_GEV is the measured Planck mass (isa/constants.py:117) |
| benchmark_lambda_cc | MENU-FITTED | — | cosmology/lambda_cc.py:57-58 — ρ_Λ=(m_e/M_Pl)⁴·α¹⁶·6, exponent 16=2·DIM_OCTONION and factor h_Cox=6 assembled from menu to hit 1.20e-123 |
| benchmark_omega_b_c | MENU-FITTED | METROLOGY | cosmology/omega_b_c.py:60-62 — 25α+75α² (25=h_v², 75=3·25 from menu) vs Planck ratio; "67 ppm" quoted against a target only known to ~% |
| benchmark_eta_B | MENU-FITTED | — | cosmology/eta_B.py:61-64 — η_B=3α⁴/14 with 3=RANK_SO7, 14=2·7 (each admitting "two readings"), matched to Planck 6.10e-10 |
| benchmark_full_ckm | MENU-FITTED | — | electroweak/substrate_ckm.py:90-101,126 — λ²=7α, A²=9/14, \|ρ̄+iη̄\|²=21α, δ=π/3 "conjecture"; every entry a small-rational×α^n form matched to PDG |
| benchmark_higgs_vev | MENU-FITTED | MEASURED-INPUT, METROLOGY | substrate_gf.py:82-84 — v=25·m_e/α²·(1+25α/(4√3)) with PDG m_e input (line 50); "27.7 ppm" vs PDG value known ~10× better with α already 7.6 ppm wrong |
| benchmark_fermi_constant | MENU-FITTED | MEASURED-INPUT, METROLOGY | substrate_gf.py:110-111 — G_F=α⁴/(625√2 m_e²(1+25α/4√3)²), same menu form + PDG m_e; +55 ppm vs G_F known to 0.5 ppm |
| benchmark_z_boson_width | STANDARD-RELABELED | MEASURED-INPUT | electroweak/decays.py:6,34 + constants.py:25,46 — standard SM tree Γ_Z with measured G_F, M_Z, effective sin²θ_W=0.23121, α_s=0.1179 all PDG inputs; substrate inert |
| benchmark_higgs_mass_vs_98gev | MENU-FITTED | MEASURED-INPUT | compute_speed.py:431-434 — λ_H=18α ("DNA integer 18") into m_h=√(2λ_H)·V_HIGGS_GEV with measured v=246.22 (electroweak/constants.py:57) |
| benchmark_muon_lifetime | MENU-FITTED | STANDARD-RELABELED | compute_speed.py:463-469 — standard tree τ=192π³/(G_F²m_μ⁵) fed entirely by menu-fitted inputs (substrate G_F + Paper-6 m_μ); honestly reported ~10% miss |
| benchmark_neutrino_sector | MEASURED-INPUT | GENUINE-FORWARD, MENU-FITTED | neutrino/__init__.py:74-80,166-167 — m_2,m_3 built from measured Δm²_21/Δm²_31 (NuFIT); m_1 and steriles are α^(N_e/2)·M_Pl forwards with no measured witness for the "0.04%" headline |
| benchmark_pmns_angles | MENU-FITTED | STANDARD-RELABELED | neutrino/__init__.py:248-252 — θ_13=asin(√(3α)) menu form scored 0.7%; θ_12, θ_23 are the textbook tri-bimaximal values (atan(1/√2), π/4) |
| benchmark_decay_constants | MENU-FITTED | MEASURED-INPUT | decay_constants.py:31-46,148-167 — f_X = PDG m_X × √(7α) or /5^(1/4); heavy states get per-state integer N∈{10,11,24,25} chosen from the menu, PDG masses as input |
| benchmark_vector_meson_decay | MENU-FITTED | MEASURED-INPUT | decay_constants.py:38-46,311-322 — each vector assigned its own rational C∈{1,7/8,7/6,10/7,15/2,40,3,7/2,25/4,8,16} against its measured f; textbook look-elsewhere case |
| benchmark_atomic_hydrogen | STANDARD-RELABELED | MEASURED-INPUT, METROLOGY | atomic/hydrogen.py:88-97,124-137 — Bohr/Rydberg/Fermi-hyperfine textbook formulas with PDG m_e, M_p, g_p (lines 57-67); only α is substrate, and its 7.6 ppm offset makes the "±15 ppm" claims dead in σ |
| benchmark_electron_anomaly | TAUTOLOGY | — | compute_speed.py:689-690 — err computed vs SCHWINGER_REFERENCE = ALPHA_SUBSTRATE/(2π), the identical formula electron_a_e_one_loop returns (hydrogen.py:163); 0.0000% by construction |
| benchmark_qcd_constants | MEASURED-INPUT | TAUTOLOGY | qcd/constants.py:29 — alpha_s = 0.1179 is the PDG 2022 literal, then compute_speed.py:723 scores it against 0.1179 (0.00%); the "K_7 Wilson loop" note is false; Λ_QCD=0.087 vs quoted 0.210 silently 59% off |
| benchmark_sin2_theta_W | MENU-FITTED | MEASURED-INPUT | electroweak/constants.py:36 — (2+α)/9 menu form; compared at compute_speed.py:750 to 1−M_W²/M_Z² built from measured M_W, M_Z |
| benchmark_black_hole_thermodynamics | STANDARD-RELABELED | MEASURED-INPUT | gravity/black_holes.py:23-51 — textbook Hawking/Schwarzschild formulas, accuracy "exact (closed form)"; G_sub from G_substrate_SI() whose chain defaults to measured ALPHA_QED (coupling.py:89,66) |
| benchmark_qed_compton_scattering | STANDARD-RELABELED | METROLOGY, MEASURED-INPUT | amplitudes/cross_sections.py:136-144 — Thomson σ=(8π/3)(α/m_e)², pure textbook with PDG m_e; ppm framing vs σ_T known far better |
| benchmark_qed_eemumu | STANDARD-RELABELED | — | amplitudes/cross_sections.py:161 — σ=4πα²/(3s), the standard QED Born formula; substrate contributes only α, and no Z-pole physics at √s=200 GeV |
| benchmark_muon_decay_rate | MENU-FITTED | STANDARD-RELABELED, MEASURED-INPUT, METROLOGY | compute_speed.py:880-889 — standard Kinoshita-Sirlin f(x)·(1+δ_QED) with PDG m_μ, m_e; the 0.011% headline is inherited from the menu-fitted substrate G_F (substrate_gf.py:111) vs τ_μ measured to 1 ppm |
| benchmark_cosmogenesis | MENU-FITTED | GENUINE-FORWARD | gravity/cosmogenesis.py:158-170 — a\*=√(1−(2/κ_p)²) with κ_p=3(1+√α)·κ_Macken assembled to land at Thorne's 0.998; f_J, κ_parent have no measured witness |
| benchmark_nmr_chemical_shifts | TAUTOLOGY | STANDARD-RELABELED, QUALITATIVE/SCORE | compute_speed.py:963-964 compares nics_sign_from_hopf_parity(n_π) to ref.predicted_sign, which nmr.py:171 generated with the same function on the same n_π — N/N guaranteed; underlying rule is Hückel 4n+2 relabeled "Hopf parity" |
| benchmark_c60_vibrational_modes | STANDARD-RELABELED | QUALITATIVE/SCORE | chemistry/vibrational.py:49-55 — hard-coded standard I_h irrep table; 174=3·60−6 and 4 T_1u/10 Raman are textbook group theory, no K_7 content |
| benchmark_composite_particles | STANDARD-RELABELED | MENU-FITTED | compositions/connected_sum.py:11,107-112 — deuteron mass_pred = m_p+m_n+0 (binding=0); additivity is trivially true to 0.12% for any weakly bound state, constituents themselves menu-fitted |
| benchmark_exotic_states | MENU-FITTED | MEASURED-INPUT | qcd/exotic_states.py:126-157 — m_X²=(4m_π⁰)²·N with per-state N from 7 to 389 including N=17,18,19 flagged "(open)"; a ~0.5%-spaced ladder can match any mass; m_π⁰ formula takes PDG m_e |
| benchmark_bhabha_scattering | STANDARD-RELABELED | — | qed/process.py:296-322 — standard LO s+t-channel Bhabha at substrate α; accuracy field is "LO QED formula", no external target scored at all |
| benchmark_moller_scattering | STANDARD-RELABELED | — | qed/process.py:246-262 — standard LO t+u-channel Møller, same as above; substrate contributes only α |
| benchmark_aromatic_resonance_energies | STANDARD-RELABELED | TAUTOLOGY, MENU-FITTED | compute_speed.py:1137 — the "experimental" coronene RE=144 is annotated "K_7-toroidal target (substrate)", i.e. the substrate's own number (12 pairs×12, aromaticity.py:47), so that entry is self-scored; rest is benzene-calibrated additive RE (standard), 56 kcal correction "calibrated from coronene/HBC" (aromaticity.py:173-174) |
| benchmark_substrate_dimensions | TAUTOLOGY | QUALITATIVE/SCORE | compute_speed.py:1172-1187 — "10/10 verified" checks constants against the definitions the same isa/constants.py file asserts at import (e.g. DIM_ADJ_SPIN7==21==N_EDGES_K7); nothing external can fail |
| benchmark_chemistry | STANDARD-RELABELED | QUALITATIVE/SCORE, TAUTOLOGY | chemistry/benchmark.py:116-145,539-552 — aromaticity classification is Hückel/Möbius counting, "14/14" NICS sign is the tautological rule above, C_60 counts are exact standard combinatorics |
| benchmark_gravitational_constant | MEASURED-INPUT | MENU-FITTED, METROLOGY | compute_speed.py:1254-1255 calls G_substrate_NNLO_natural()/G_substrate_SI() with no args, whose defaults are alpha=ALPHA_QED — the measured CODATA α (gravity/coupling.py:66,89,129) — plus measured m_e; the "-11 ppm" is a menu bracket (8/7)α^(21/2)(1+α/7+(21/8)α²) tuned inside CODATA's ±22 ppm band |

## Summary

- **Primary-class counts (38 benchmarks): MENU-FITTED 18, STANDARD-RELABELED
  13, TAUTOLOGY 3, MEASURED-INPUT 3, GENUINE-FORWARD 1.**
- Every ppm-class claim (α, v_EW, G_F, G, hydrogen, Thomson) is
  METROLOGY-contaminated: with ALPHA_SUBSTRATE excluded at ~50,000σ, all
  α-propagated ppm matches are excluded at thousands of σ in experimental σ.
- Three benchmarks are outright self-comparisons (electron_anomaly,
  nmr_chemical_shifts, substrate_dimensions), plus a partial fourth
  (aromatic_resonance_energies' coronene reference is the substrate's own
  number annotated as the "experimental" target).
- Notable false provenance note: benchmark_qcd_constants hardcodes the PDG
  α_s = 0.1179 literal, scores it against itself, and its notes attribute it
  to a "K_7 Wilson loop"; its Λ_QCD is silently 59% off the value its own
  comment quotes.
- Lower-confidence primaries (reasonable alternative in parentheses):
  muon_decay_rate MENU-FITTED (STANDARD-RELABELED), qcd_constants
  MEASURED-INPUT (TAUTOLOGY), neutrino_sector MEASURED-INPUT
  (GENUINE-FORWARD), composite_particles STANDARD-RELABELED (MENU-FITTED).
- **The only entry untainted as a test is benchmark_wimp_tower's unmeasured
  DM rungs** (with the sterile-ν/f_J sub-outputs as secondaries) — and even
  those carry free per-rung integers N_e that will absorb future data unless
  frozen BEFORE any measurement exists.  Freezing them is the single
  perishable action item this triage produces.

## Consequence

The suite README's central claim — "38 substantive physical observables …
genuinely predictive — not fitted, not numerology" — is quantitatively
inverted by its own contents: 18 fitted, 13 standard results relabeled, 3
self-scored, 3 measured-input, 1 forward-untested.  No benchmark demonstrates
discriminating substrate content under the appropriate gauntlet test.  This
document is evidence for the Auditor brief and a blocking input to any
de-editorializing pass on the README and benchmark notes.
