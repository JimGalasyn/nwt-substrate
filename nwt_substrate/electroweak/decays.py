"""
Z-boson partial decay widths and total width.

For each kinematically open channel f f̄:

    Γ(Z → f f̄) = (g_Z² / 48 π) M_Z × (g_V^f² + g_A^f²) × N_c × β_f × (1 + 2 m_f²/M_Z²)

with β_f = √(1 - 4 m_f² / M_Z²).  Summing over all SM fermions reproduces
the PDG Z width Γ_Z = 2.495 GeV to ~1% (residual from radiative corrections
and m_b mass-effect).
"""

from __future__ import annotations
import math

from .constants import M_Z, G_Z
from .couplings import coupling, SM_COUPLINGS


# Approximate fermion masses (GeV); only used for the β_f phase-space factor
# of the heavy fermions
FERMION_MASS_GEV = {
    "e": 0.511e-3, "mu": 0.10566, "tau": 1.77686,
    "nu_e": 0.0, "nu_mu": 0.0, "nu_tau": 0.0,
    "u": 2.16e-3, "d": 4.67e-3, "s": 0.0934,
    "c": 1.27, "b": 4.18, "t": 172.7,
}


def partial_width_Z(fermion_name: str) -> float:
    """
    Γ(Z → f f̄) [GeV].  Returns 0.0 if 2 m_f > M_Z (kinematically forbidden,
    e.g. top quark).
    """
    cpl = coupling(fermion_name)
    m_f = FERMION_MASS_GEV[fermion_name]
    if 2 * m_f >= M_Z:
        return 0.0
    beta_f = math.sqrt(1.0 - 4.0 * m_f ** 2 / M_Z ** 2)
    phase = beta_f * (1.0 + 2.0 * m_f ** 2 / M_Z ** 2)
    prefactor = (G_Z ** 2 / (48.0 * math.pi)) * M_Z
    return prefactor * cpl.gV2_plus_gA2 * cpl.n_color * phase


def total_width_Z() -> float:
    """Sum of partial widths over all SM fermions kinematically open at √s = M_Z."""
    return sum(partial_width_Z(name) for name in SM_COUPLINGS)


def branching_ratios_Z() -> dict:
    """B(Z → f f̄) for each open channel."""
    Gamma_total = total_width_Z()
    return {name: partial_width_Z(name) / Gamma_total for name in SM_COUPLINGS}


def width_summary() -> str:
    """Pretty-print Z partial widths and total width."""
    lines = ["Z-boson partial widths (LO, including m_f phase-space β·(1+2m²/s)):"]
    lines.append("")
    Gamma_total = total_width_Z()
    BR = branching_ratios_Z()
    fermion_groups = [
        ("Charged leptons", ["e", "mu", "tau"]),
        ("Neutrinos       ", ["nu_e", "nu_mu", "nu_tau"]),
        ("Up-type quarks  ", ["u", "c", "t"]),
        ("Down-type quarks", ["d", "s", "b"]),
    ]
    for group, members in fermion_groups:
        sub_total = sum(partial_width_Z(m) for m in members)
        lines.append(f"  {group}:  Γ = {sub_total*1000:>7.2f} MeV  "
                     f"(BR = {sub_total / Gamma_total * 100:>5.2f}%)")
        for m in members:
            gamma = partial_width_Z(m)
            if gamma > 0:
                lines.append(f"      Γ(Z → {m:<7s}) = "
                             f"{gamma*1000:>7.2f} MeV   "
                             f"({gamma / Gamma_total * 100:>5.2f}%)")
            else:
                lines.append(f"      Γ(Z → {m:<7s}) =     -      "
                             f"(2 m > M_Z, kinematically forbidden)")
    lines.append("")
    lines.append(f"  TOTAL Γ_Z (predicted) = {Gamma_total:.4f} GeV "
                 f"({Gamma_total*1000:.1f} MeV)")
    lines.append(f"  PDG Γ_Z (measured)    = 2.4952 GeV (2495.2 MeV)")
    lines.append(f"  agreement: {Gamma_total / 2.4952 * 100:.2f}% of PDG")
    return "\n".join(lines)
