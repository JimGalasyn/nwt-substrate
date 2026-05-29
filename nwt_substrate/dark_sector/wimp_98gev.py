"""
NWT 98 GeV WIMP via L_NWT Higgs portal — direct-detection + LHC predictions.

THE WIMP
========

The substrate's K_8 mass tower (VV's enumeration: 20 N_e rungs spanning
~14 meV to ~M_Pl) has a rung at

    m_χ = α^8 · M_Pl ≈ 98.23 GeV    (N_e = 16, the "Higgs-sector" exponent)

This sits in the WIMP-band of dark-matter direct-detection experiments.
Other rungs of the same K_8 tower predict additional species at 38 keV
(warm DM), 1.15 TeV, and beyond.

L_NWT HIGGS PORTAL
==================

The minimal substrate-consistent coupling of the K_8 DM scalar χ to the
K_7 SM sector goes via the L_NWT complex scalar ψ (Paper 16, L_1):

    L_portal  =  - g_Hχχ · |ψ|² · χ²

After electroweak symmetry breaking |ψ|² → v² + 2 v · h + h²,

    L_portal  =  - g_eff · h · χ²   ,   g_eff  ≡  g_Hχχ · v

with v = 246 GeV the Higgs VEV.

KEY OBSERVABLES (all derive from g_eff at tree level)
=====================================================

1. Direct-detection (LZ-G3 / XENONnT):
        σ_SI  =  (g_eff² · f_N²) / (π · m_h⁴) · μ_χN²

   where m_h = 125 GeV, f_N ≈ 0.30 (Higgs-nucleon form factor; lattice +
   chiral pert. theory), μ_χN = m_χ·m_N/(m_χ+m_N) is the reduced mass.

2. LHC production via Higgs portal (pp → h* → χχ + jet):
        σ_LHC  ~  σ_h^SM(at √s_eff) · BR(h → χχ)
   with BR(h→χχ) determined by g_eff if 2 m_χ < m_h (kinematically
   forbidden for our 98 GeV WIMP since 2·98 > 125, so production goes
   through OFF-SHELL Higgs — much smaller).

3. Relic abundance via geometric annihilation:
        σ_ann · v  ~  π · (1/m_χ)²
   geometrically saturated; tunes to right Ω_DM h² (WIMP miracle).

THE OPEN INPUT
==============

g_Hχχ is the substrate-derived dimensionless Higgs-portal coupling.
Natural candidates from substrate algebra (σ_SI values below are the
`sigma_si_higgs_portal` output for each tier — NOT round-number lore;
see the parity check in `test_dark_sector_wimp.py`):
  - O(1):    naive Yukawa scale       → σ_SI ≈ 2.4e-33 cm² (EXCLUDED, ratio ~6e14)
  - α ≈ 1/137: EM-suppressed          → σ_SI ≈ 1.3e-37 cm² (EXCLUDED, ratio ~3e10)
  - α²:      double-EM suppression    → σ_SI ≈ 6.8e-42 cm² (EXCLUDED, ratio ~1.6e6) [default]
  - α^4 = α^(N_e/4): K_8-mass-scaled  → σ_SI ≈ 1.9e-50 cm² (ALLOWED, ~4.6e-3 of limit)

The current LZ (2024) limit at m_χ = 100 GeV is σ_SI < 4.2 × 10⁻⁴⁸ cm².
Only the α^4 (or smaller) coupling survives this bound: σ_SI scales as
g_Hχχ², so each factor of α drops σ_SI by α² ≈ 5e-5.

This module computes σ_SI as a function of g_Hχχ; the substrate
program's task is to fix g_Hχχ from K_7/K_8 portal algebra.

REFERENCES
==========
- VV's K_8 mass-tower enumeration (substrate-dna-integers-18-and-25 memory)
- Paper 16 L_NWT (paper16_nwt_lagrangian.tex), eq. master L_NWT
- LZ-2024 results: arXiv:2410.17036
- XENONnT: arXiv:2303.14729
"""

from __future__ import annotations
import math
from dataclasses import dataclass

from nwt_substrate.isa.constants import ALPHA_SUBSTRATE, M_PL_GEV

# Standard Model parameters
V_HIGGS_GEV: float = 246.22                  # Higgs VEV
M_HIGGS_GEV: float = 125.10                  # observed Higgs mass
M_NUCLEON_GEV: float = 0.939                 # average proton/neutron
F_N_HIGGS: float = 0.30                      # Higgs-nucleon form factor
                                              # (lattice + chiral pert; FLAG-like)


# Unit conversion: 1 GeV⁻² in cm²
# ℏc = 0.1973 GeV·fm = 0.1973e-13 GeV·cm
# 1 GeV⁻¹ = ℏc/GeV cm = 0.1973e-13 cm
# 1 GeV⁻² = (0.1973e-13)² cm² = 3.894e-28 cm²
GEV_INV2_TO_CM2: float = (0.19732698e-13) ** 2     # = 3.894e-28 cm²


@dataclass(frozen=True)
class WIMP_98GeV:
    """98 GeV WIMP from the K_8 mass tower (N_e = 16 rung).

    Mass derived from m_χ = α^N · M_Pl with N = N_e/2 = 8.
    """
    # Substrate-fixed inputs
    mass_gev: float = ALPHA_SUBSTRATE ** 8 * M_PL_GEV
    n_e: int = 16
    k8_sector: str = "DM-tower (carrier knot in K_8 not K_7)"

    # Higgs-portal coupling — substrate-derived but currently parametric
    g_hxx: float = ALPHA_SUBSTRATE ** 2          # default: α² K_7/K_8 portal suppression

    @property
    def g_eff(self) -> float:
        """g_eff = g_Hχχ · v_Higgs.  Coefficient of h·χ² in expanded L_portal."""
        return self.g_hxx * V_HIGGS_GEV

    @property
    def reduced_mass_with_nucleon(self) -> float:
        return self.mass_gev * M_NUCLEON_GEV / (self.mass_gev + M_NUCLEON_GEV)


def higgs_portal_lagrangian_text() -> str:
    """Return the L_NWT Higgs-portal Lagrangian as a string."""
    return r"""
    L_portal  =  - g_Hχχ · |ψ|² · χ²

    After EWSB:  |ψ|² → v² + 2 v h + h²

    L_portal  →  - g_eff · h · χ²    with  g_eff = g_Hχχ · v_Higgs
    """


# ---------------------------------------------------------------------------
# Direct-detection cross section
# ---------------------------------------------------------------------------

def sigma_si_higgs_portal(wimp: WIMP_98GeV) -> float:
    """Spin-independent χ-nucleon cross section from Higgs-portal coupling.

        σ_SI = (g_eff² · f_N²) / (π · m_h⁴) · μ_χN²

    Returns σ_SI in cm².  At tree-level, no loop corrections.
    """
    mu = wimp.reduced_mass_with_nucleon
    sigma_gev2 = (wimp.g_eff ** 2 * F_N_HIGGS ** 2 * mu ** 2
                  / (math.pi * M_HIGGS_GEV ** 4))
    return sigma_gev2 * GEV_INV2_TO_CM2


def compare_with_lz_limit(sigma_si_cm2: float,
                          lz_limit_at_100gev_cm2: float = 4.2e-48) -> dict:
    """Compare a predicted σ_SI to the LZ-2024 limit at m_χ ≈ 100 GeV.

    Returns dict with predicted, limit, ratio, and verdict ('allowed' or
    'excluded')."""
    return {
        "sigma_si_cm2": sigma_si_cm2,
        "lz_limit_cm2": lz_limit_at_100gev_cm2,
        "ratio_pred_over_limit": sigma_si_cm2 / lz_limit_at_100gev_cm2,
        "verdict": "EXCLUDED" if sigma_si_cm2 > lz_limit_at_100gev_cm2 else "ALLOWED",
    }


# ---------------------------------------------------------------------------
# LHC production
# ---------------------------------------------------------------------------

def lhc_production_cross_section(wimp: WIMP_98GeV,
                                 root_s_gev: float = 13600.0) -> float:
    """Rough estimate of pp → χχ via off-shell Higgs at LHC √s=13.6 TeV.

    For our 98 GeV WIMP, 2 m_χ ≈ 196 GeV > m_h = 125 GeV, so h → χχ is
    kinematically FORBIDDEN on-shell.  Production goes through OFF-shell
    h* with effective mass ≥ 2 m_χ:

        σ(pp → h* → χχ) ~ (g_eff²/8π) · ∫ dM²
            σ_h^SM(M)  ·  β(M)/((M² - m_h²)² + (m_h Γ_h)²)

    where β(M) = √(1 - 4 m_χ²/M²) is the χχ phase-space.

    Estimate via narrow-effective-width approximation: scale of σ_LHC
    ≈ (g_eff·v)² · (1 / (4 m_χ² - m_h²)²) · σ_h^SM × phase space.

    Returns σ_LHC in fb.  HONEST: a real prediction needs MadGraph with
    the substrate Lagrangian as a UFO model.
    """
    sigma_h_sm_lhc_fb = 55_000.0                   # σ(pp→h) at √s=13.6 TeV (fb)
    # Off-shell propagator at effective mass M ≈ 2 m_χ
    delta_m2 = (2 * wimp.mass_gev) ** 2 - M_HIGGS_GEV ** 2
    propagator_suppression = M_HIGGS_GEV ** 4 / max(delta_m2 ** 2, 1e-6)
    # Coupling strength
    coupling_factor = wimp.g_eff ** 2 / (8 * math.pi)
    # PDF + phase-space factor (rough O(1/m_χ²))
    phase_space = 1.0 / wimp.mass_gev ** 2
    return sigma_h_sm_lhc_fb * propagator_suppression * coupling_factor * phase_space


# ---------------------------------------------------------------------------
# main: print structured summary
# ---------------------------------------------------------------------------

def predict_all(wimp: WIMP_98GeV | None = None,
                show_table: bool = True) -> dict:
    """Return predictions + comparison for a given WIMP configuration."""
    if wimp is None:
        wimp = WIMP_98GeV()

    sigma_si = sigma_si_higgs_portal(wimp)
    lz_cmp = compare_with_lz_limit(sigma_si)
    sigma_lhc_fb = lhc_production_cross_section(wimp)

    predictions = {
        "wimp": {
            "mass_gev": wimp.mass_gev,
            "g_hxx": wimp.g_hxx,
            "g_eff": wimp.g_eff,
            "n_e": wimp.n_e,
            "k8_sector": wimp.k8_sector,
        },
        "direct_detection": {
            "sigma_si_cm2": sigma_si,
            "lz_comparison": lz_cmp,
        },
        "lhc": {
            "sigma_offshell_fb_at_13_6tev": sigma_lhc_fb,
        },
    }

    if show_table:
        print(f"  m_χ                = {wimp.mass_gev:>9.3f} GeV    "
              f"(α^8 · M_Pl, K_8 N_e={wimp.n_e})")
        print(f"  g_Hχχ              = {wimp.g_hxx:>9.3e}      "
              f"(substrate input)")
        print(f"  g_eff = g·v        = {wimp.g_eff:>9.3e}")
        print(f"  σ_SI (Higgs portal) = {sigma_si:>9.3e} cm²")
        print(f"  LZ-2024 @ 100 GeV  < {lz_cmp['lz_limit_cm2']:>9.3e} cm²    "
              f"{lz_cmp['verdict']}")
        print(f"  ratio (pred/limit) = {lz_cmp['ratio_pred_over_limit']:>9.3e}")
        print(f"  σ_LHC (off-shell)  = {sigma_lhc_fb:>9.3e} fb @ √s=13.6 TeV")

    return predictions
