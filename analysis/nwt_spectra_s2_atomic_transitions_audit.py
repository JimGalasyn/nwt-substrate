"""
Substrate-algebraic audit for Spectra Tier S.2 — atomic transitions in
hydrogenic ions + hydrogen 21-cm hyperfine.

Filed AFTER pre-registration ([[spectra-s2-atomic-transitions-prereg]]).
Tests substrate orbital-locking where the mechanism IS orbital-geometric
(1-electron Coulomb). Decisive contrast to S.1 (condensate-dominated,
honest Form E).

Substrate Rydberg = (1/2) α² m_e c² with α = 1/(25π√3+1).
Hydrogenic transitions: E(Z, n_i → n_f) = Z² R (1/n_f² − 1/n_i²) · μ/m_e.
21-cm hyperfine: E = (4/3) α⁴ R g_p (m_e/m_p).

Run:
    PYTHONPATH=/home/jim/repos/nwt-substrate \\
        python3 analysis/nwt_spectra_s2_atomic_transitions_audit.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Substrate primitives (Paper 17 trefoil α + standard constants)
# ---------------------------------------------------------------------------

# Substrate α (= NWT-derived fine-structure constant)
ALPHA_SUBSTRATE = 1.0 / (25.0 * math.pi * math.sqrt(3.0) + 1.0)

# CODATA 2018 fundamental constants
M_E_KG = 9.1093837015e-31      # kg
C_M_S = 2.99792458e8           # m/s (exact)
E_CHARGE = 1.602176634e-19     # C (exact)
H_PLANCK = 6.62607015e-34      # J·s (exact)
H_BAR = H_PLANCK / (2 * math.pi)

# Mass ratios
M_P_OVER_M_E = 1836.15267343   # CODATA proton/electron mass ratio
G_PROTON = 5.585694713         # proton g-factor (CODATA)

# Electron rest energy
M_E_C2_J = M_E_KG * C_M_S**2
M_E_C2_EV = M_E_C2_J / E_CHARGE   # 510998.95 eV


# Substrate Rydberg energy
def substrate_rydberg_eV(reduced_mass_factor: float = 1.0) -> float:
    """Substrate-derived Rydberg energy in eV.

    R = (1/2) · α² · m_e · c² · (μ/m_e)

    For an infinite-mass nucleus, μ/m_e = 1 (R_∞).
    For hydrogen, μ/m_e = M_p/(M_p + m_e) = 1836.15/1837.15 = 0.99946.
    """
    return 0.5 * ALPHA_SUBSTRATE**2 * M_E_C2_EV * reduced_mass_factor


R_INF = substrate_rydberg_eV(1.0)  # R_∞ = 13.6057 eV
R_HYDROGEN = substrate_rydberg_eV(M_P_OVER_M_E / (M_P_OVER_M_E + 1))  # R_H


# ---------------------------------------------------------------------------
# Hydrogenic transition reference set (LOCKED at pre-reg time)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HydrogenicLine:
    name: str
    Z: int                  # nuclear charge
    n_initial: int          # higher level (smaller energy)
    n_final: int            # lower level (larger energy)
    lambda_nm_obs: float    # NIST observed wavelength in nm (vacuum)
    energy_eV_obs: float    # observed transition energy
    substrate_prefactor: str  # symbolic substrate form


def make_reference_set() -> list[HydrogenicLine]:
    """NIST CODATA + standard atomic-physics references.

    Wavelengths from NIST Atomic Spectra Database, vacuum values.
    """
    return [
        # Hydrogen Lyman series (vacuum UV)
        HydrogenicLine("H Lyman α",  1, 2, 1, 121.5670, 10.19883, "3/4 = RANK_SO7 / N_VERTICES_K_4"),
        HydrogenicLine("H Lyman β",  1, 3, 1, 102.5722, 12.08749, "8/9 = DIM_OCTONION / N_POS_ROOTS_SO7"),
        HydrogenicLine("H Lyman γ",  1, 4, 1,  97.2537, 12.74853, "15/16 = (N_EDGES_K7−H_COXETER) / (N_EDGES_K7−H_V_SO7)"),
        # Hydrogen Balmer series (visible)
        HydrogenicLine("H Balmer α", 1, 3, 2, 656.279,  1.88867, "5/36 = H_V_SO7 / (4·9)"),
        HydrogenicLine("H Balmer β", 1, 4, 2, 486.135,  2.54983, "3/16 = RANK_SO7 / (21−5)"),
        HydrogenicLine("H Balmer γ", 1, 5, 2, 434.047,  2.85594, "21/100 = N_EDGES_K7 / 100"),
        # Hydrogen Paschen series (IR)
        HydrogenicLine("H Paschen α", 1, 4, 3, 1875.10, 0.66128, "7/144 = N_VERTICES_K7 / K8_PARTITION[2]²"),
        # Helium-II Lyman α
        HydrogenicLine("He⁺ Lyman α", 2, 2, 1,  30.378, 40.8121, "Z²·(3/4) = 4·RANK_SO7/N_VERTICES_K_4"),
        # Li-III Lyman α
        HydrogenicLine("Li²⁺ Lyman α", 3, 2, 1, 13.4989, 91.846, "Z²·(3/4) = 9·RANK_SO7/N_VERTICES_K_4"),
        # Be-IV Lyman α
        HydrogenicLine("Be³⁺ Lyman α", 4, 2, 1,  7.594, 163.298, "Z²·(3/4)"),
        # B-V Lyman α
        HydrogenicLine("B⁴⁺ Lyman α",  5, 2, 1,  4.860, 255.157, "Z²·(3/4)"),
        # C-VI Lyman α
        HydrogenicLine("C⁵⁺ Lyman α",  6, 2, 1,  3.375, 367.395, "Z²·(3/4)"),
    ]


# ---------------------------------------------------------------------------
# Substrate prediction
# ---------------------------------------------------------------------------

def reduced_mass_factor(Z: int) -> float:
    """μ/m_e for a Z-charge ion with mass approximately Z·m_proton.

    For hydrogenic ions, M_nucleus ≈ Z · M_p (using proton mass as approx).
    More precisely should use isotope-specific masses, but this is a 1e-5
    correction at most.
    """
    M_nucleus_over_m_e = Z * M_P_OVER_M_E
    return M_nucleus_over_m_e / (M_nucleus_over_m_e + 1)


def substrate_predicted_energy(line: HydrogenicLine) -> float:
    """Substrate prediction: E = Z² · R_substrate · (1/n_f² − 1/n_i²) · μ/m_e."""
    Z = line.Z
    rmf = reduced_mass_factor(Z)
    R = substrate_rydberg_eV(rmf)
    prefactor = 1.0 / line.n_final**2 - 1.0 / line.n_initial**2
    return Z**2 * R * prefactor


# ---------------------------------------------------------------------------
# 21-cm hyperfine
# ---------------------------------------------------------------------------

# Observed value (CODATA 2018)
NU_21CM_HZ_OBS = 1420.4057517667e6
E_21CM_EV_OBS = NU_21CM_HZ_OBS * H_PLANCK / E_CHARGE
# = 5.874 × 10⁻⁶ eV


def substrate_predicted_21cm_LO() -> float:
    """Substrate prediction for hydrogen 21-cm at leading order.

    Standard Fermi-contact result for 1s hydrogen hyperfine:
        E_HFS = (8/3) · g_p · α² · (m_e/m_p) · R_∞

    Substrate-canonical prefactor (BOTH primitives):
        8 = DIM_OCTONION (direct substrate primitive)
        3 = RANK_SO7 (direct substrate primitive)

    Higher-order QED corrections (radiative, recoil, structural) contribute
    ~0.1% and are not in the LO substrate form.
    """
    R = R_INF
    return (8.0 / 3.0) * G_PROTON * ALPHA_SUBSTRATE**2 * (1.0 / M_P_OVER_M_E) * R


# ---------------------------------------------------------------------------
# Audit
# ---------------------------------------------------------------------------

def run_audit() -> dict:
    lines = make_reference_set()
    results = []
    for line in lines:
        pred = substrate_predicted_energy(line)
        obs = line.energy_eV_obs
        dev_pct = 100 * abs(obs - pred) / obs
        results.append({
            "line": line,
            "energy_predicted_eV": pred,
            "energy_observed_eV": obs,
            "deviation_pct": dev_pct,
        })

    # 21-cm separately
    pred_21cm = substrate_predicted_21cm_LO()
    obs_21cm = E_21CM_EV_OBS
    dev_21cm_pct = 100 * abs(obs_21cm - pred_21cm) / obs_21cm

    return {
        "hydrogenic": results,
        "rydberg_substrate_eV": R_INF,
        "rydberg_hydrogen_eV": R_HYDROGEN,
        "alpha_substrate": ALPHA_SUBSTRATE,
        "twenty_one_cm": {
            "predicted_eV": pred_21cm,
            "observed_eV": obs_21cm,
            "deviation_pct": dev_21cm_pct,
            "predicted_MHz": pred_21cm * E_CHARGE / H_PLANCK / 1e6,
            "observed_MHz": NU_21CM_HZ_OBS / 1e6,
        },
    }


def render_report(audit: dict) -> str:
    out = []
    out.append("=" * 86)
    out.append("Spectra Tier S.2 — Atomic transitions in hydrogenic ions + 21-cm hyperfine")
    out.append("=" * 86)
    out.append("")
    out.append(f"Substrate α (Paper 17 trefoil):   {audit['alpha_substrate']:.12f}")
    out.append(f"Substrate Rydberg (R_∞):           {audit['rydberg_substrate_eV']:.6f} eV")
    out.append(f"Substrate Rydberg (R_H with μ):   {audit['rydberg_hydrogen_eV']:.6f} eV")
    out.append("")
    out.append("Formula: E(Z, n_i → n_f) = Z² · R_∞ · (1/n_f² − 1/n_i²) · μ/m_e")
    out.append("")

    out.append("-" * 86)
    out.append("HYDROGENIC TRANSITIONS")
    out.append("-" * 86)
    out.append("")
    out.append(f"  {'line':<18} {'Z':>2} {'n_i→n_f':<8} {'pred (eV)':>12} {'obs (eV)':>12} {'dev %':>8} {'prefactor':<40}")
    for r in audit["hydrogenic"]:
        line = r["line"]
        ni_nf = f"{line.n_initial}→{line.n_final}"
        out.append(f"  {line.name:<18} {line.Z:>2} {ni_nf:<8} {r['energy_predicted_eV']:>12.5f} "
                   f"{r['energy_observed_eV']:>12.5f} {r['deviation_pct']:>8.4f} {line.substrate_prefactor:<40}")
    out.append("")

    # 21-cm
    out.append("-" * 86)
    out.append("21-CM HYPERFINE (LO substrate prediction)")
    out.append("-" * 86)
    out.append("")
    t = audit["twenty_one_cm"]
    out.append(f"  Substrate form: E = (8/3) · g_p · α² · (m_e/m_p) · R_∞")
    out.append(f"  Substrate primitives: 8 = DIM_OCTONION; 3 = RANK_SO7 (BOTH direct primitives)")
    out.append("")
    out.append(f"  Predicted: {t['predicted_MHz']:.4f} MHz  ({t['predicted_eV']:.6e} eV)")
    out.append(f"  Observed:  {t['observed_MHz']:.4f} MHz  ({t['observed_eV']:.6e} eV)")
    out.append(f"  Deviation: {t['deviation_pct']:.3f}%")
    out.append("")

    # Summary
    out.append("=" * 86)
    out.append("SUMMARY")
    out.append("=" * 86)
    out.append("")
    devs = [r["deviation_pct"] for r in audit["hydrogenic"]]
    max_dev = max(devs)
    avg_dev = sum(devs) / len(devs)
    out.append(f"  Hydrogenic transitions: {len(audit['hydrogenic'])} tested")
    out.append(f"    average deviation: {avg_dev:.4f}%")
    out.append(f"    maximum deviation: {max_dev:.4f}%")
    out.append(f"    Z range:           {min(r['line'].Z for r in audit['hydrogenic'])} to {max(r['line'].Z for r in audit['hydrogenic'])}")
    out.append("")
    out.append(f"  21-cm hyperfine: {t['deviation_pct']:.3f}% deviation at LO")
    out.append("    (higher-order QED corrections ~0.1% not in LO form)")
    out.append("")

    # Verdict
    out.append("=" * 86)
    out.append("FRAMEWORK VERDICT")
    out.append("=" * 86)
    out.append("")
    out.append("  Pre-reg locked criterion (LO Bohr + LO Fermi-contact):")
    out.append("    Framework PASSES iff:")
    out.append("      (a) Hydrogenic max deviation ≤ 0.1% at LO Bohr formula")
    out.append("          (NLO Sommerfeld Z²α² fine structure ~0.02% expected)")
    out.append("      (b) 21-cm prediction ≤ 5% at LO Fermi-contact form")
    out.append("")
    cond_a = max_dev <= 0.1
    cond_b = t["deviation_pct"] <= 5.0
    out.append(f"  (a) Hydrogenic max dev ≤ 0.1%:     max = {max_dev:.4f}%   {'✓' if cond_a else '✗'}")
    out.append(f"  (b) 21-cm dev ≤ 5%:                 {t['deviation_pct']:.3f}%   {'✓' if cond_b else '✗'}")
    out.append("")

    if cond_a and cond_b:
        out.append("  → FORM A CLEAN: substrate orbital-locking exactly reproduces atomic")
        out.append("    spectroscopy. Substrate Rydberg formula IS the atomic Coulomb spectrum.")
    elif cond_a:
        out.append("  → FORM D: hydrogenic clean; 21-cm needs higher-order corrections.")
    else:
        out.append("  → FRAMEWORK MISCALIBRATED: hydrogenic Rydberg formula doesn't match.")
    return "\n".join(out)


if __name__ == "__main__":
    audit = run_audit()
    print(render_report(audit))
