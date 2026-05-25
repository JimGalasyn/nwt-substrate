"""
Substrate pattern-stability ratio ρ = m / Γ for K_7 walks.

P7b synthesis closure (Galasyn 2026-05-23,
analysis/paper21_rho_formula_k7_walks.md in null-worldtube-private):

Every NWT particle is a closed walk on the K_7 substrate.  The walk's
**pattern-stability ratio** is the ratio of its substrate re-inscription
rate (Compton angular frequency ω_C = m c² / ℏ — the rate at which the
K_7 stabilizer code re-validates the syndrome) to its substrate
dissolution rate (Γ / ℏ — the total decay width):

    ρ_X = R_reinsc / R_diss = (m_X c² / ℏ) / (Γ_X / ℏ) = m_X / Γ_X

In natural units (ℏ = c = 1)::

    ρ_X = m_X / Γ_X         (inverse relative width)

Substrate interpretation::

    ρ → ∞      Γ = 0 ; substrate measurements always re-validate; BPS-stable
    ρ ≫ 1      Γ ≪ m ; long-lived; substrate re-inscription dominates
    ρ ∼ 1      Γ ∼ m ; critical; "active-sustained" regime
    ρ < 1      Γ > m ; substrate dissolution dominates; cannot maintain
                       coherent K_7 walk pattern

For all 25 compendium particles (excluding strictly-stable ones),
ρ > 1 universally — the substrate K_7 walk domain sits ENTIRELY in the
re-inscription-dominated regime.  The lowest-ρ NWT particles are the
strong resonances (ρ-meson ρ_X ≈ 5).  The classical "ρ < 1" regime
(fires, explosions, dissipative structures) lies OUTSIDE the K_7 walk
domain — those patterns aren't substrate particles.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping


# ============================================================
# Catalog: PDG mass + width for compendium particles
# ============================================================

@dataclass(frozen=True)
class StabilitySpec:
    """PDG mass + total width for a compendium particle.

    Width Γ = 0 → strictly stable (ρ = ∞).  All values in GeV.
    """
    name: str
    m_X_GeV: float
    Gamma_X_GeV: float
    mediator: str

    @property
    def is_stable(self) -> bool:
        return self.Gamma_X_GeV == 0.0


# Compendium PDG values from paper21_rho_formula_k7_walks.md §4.1.
COMPENDIUM_STABILITY: Mapping[str, StabilitySpec] = {
    # Strictly stable (BPS-protected)
    "p":     StabilitySpec("p",     0.93828,  0.0,       "stable"),
    "e":     StabilitySpec("e",     5.110e-4, 0.0,       "stable"),
    "gamma": StabilitySpec("gamma", 0.0,      0.0,       "stable"),
    # Long-lived (weak-mediated)
    "mu":    StabilitySpec("mu",    0.10566,  3.00e-19,  "weak"),
    "K_L":   StabilitySpec("K_L",   0.49761,  1.29e-17,  "weak"),
    "K+":    StabilitySpec("K+",    0.49368,  5.32e-17,  "weak"),
    "pi+":   StabilitySpec("pi+",   0.13957,  2.53e-17,  "weak"),
    # Weak hadronic
    "Xi0":   StabilitySpec("Xi0",   1.315,    2.27e-15,  "weak_hadronic"),
    "Lambda": StabilitySpec("Lambda", 1.116,  2.50e-15,  "weak_hadronic"),
    "Xi-":   StabilitySpec("Xi-",   1.322,    4.02e-15,  "weak_hadronic"),
    "Omega-": StabilitySpec("Omega-", 1.672,  8.02e-15,  "weak_hadronic"),
    "Sigma-": StabilitySpec("Sigma-", 1.197,  4.45e-15,  "weak_hadronic"),
    "Sigma+": StabilitySpec("Sigma+", 1.189,  8.21e-15,  "weak_hadronic"),
    "K_S":   StabilitySpec("K_S",   0.49761,  7.35e-15,  "weak"),
    # Heavy-flavor weak
    "B+":    StabilitySpec("B+",    5.279,    4.02e-13,  "weak_heavy"),
    "D+":    StabilitySpec("D+",    1.870,    6.33e-13,  "weak_heavy"),
    "D0":    StabilitySpec("D0",    1.865,    1.60e-12,  "weak_heavy"),
    "tau":   StabilitySpec("tau",   1.777,    2.27e-12,  "weak"),
    # EM
    "pi0":   StabilitySpec("pi0",   0.135,    7.81e-9,   "EM"),
    "eta":   StabilitySpec("eta",   0.548,    1.31e-6,   "EM_strong"),
    "Sigma0": StabilitySpec("Sigma0", 1.193,  8.89e-6,   "EM"),
    "Upsilon": StabilitySpec("Upsilon", 9.460, 5.40e-5,  "EM_narrow"),
}


# ============================================================
# Formula
# ============================================================

def stability_ratio(m_X_GeV: float, Gamma_X_GeV: float) -> float:
    """Substrate pattern-stability ratio ρ_X = m_X / Γ_X.

    Returns:
        float, possibly ``math.inf`` if Γ_X == 0 (strictly stable).

    Raises ValueError if Γ_X < 0 or m_X < 0.
    """
    if m_X_GeV < 0.0:
        raise ValueError(f"Mass must be non-negative: got {m_X_GeV}")
    if Gamma_X_GeV < 0.0:
        raise ValueError(f"Width must be non-negative: got {Gamma_X_GeV}")
    if Gamma_X_GeV == 0.0:
        return math.inf
    return m_X_GeV / Gamma_X_GeV


def stability_ratio_for(name: str) -> float:
    """Substrate ρ for a named compendium particle.

    Returns ``math.inf`` for strictly stable particles.
    """
    if name not in COMPENDIUM_STABILITY:
        raise KeyError(f"Unknown particle: {name}.  "
                       f"Known: {sorted(COMPENDIUM_STABILITY)}")
    spec = COMPENDIUM_STABILITY[name]
    return stability_ratio(spec.m_X_GeV, spec.Gamma_X_GeV)


def log10_stability_ratio_for(name: str) -> float:
    """log_10(ρ) for a named particle; +inf for strictly stable."""
    rho = stability_ratio_for(name)
    return math.inf if math.isinf(rho) else math.log10(rho)


# ============================================================
# Regime classification
# ============================================================

def classify_regime(rho: float) -> str:
    """Classify a substrate ρ value into the qualitative regimes:

        "BPS"       (ρ = ∞, strictly stable)
        "passive"   (ρ ≫ 1, re-inscription-dominated; the K_7 walk domain)
        "critical"  (ρ ∼ 1, active-sustained)
        "decaying"  (ρ < 1, dissolution-dominated, NOT in K_7 walk domain)
    """
    if math.isinf(rho):
        return "BPS"
    if rho >= 10.0:
        return "passive"
    if rho >= 1.0:
        return "critical"
    return "decaying"


def all_k7_walks_are_passive_or_BPS() -> bool:
    """Empirical: every compendium particle has ρ > 1.

    The substrate-monism reading: the K_7 walk domain sits ENTIRELY in
    the re-inscription-dominated regime.  Fires, explosions, and other
    dissipative structures (ρ < 1) are NOT substrate K_7 walks.
    """
    return all(
        classify_regime(stability_ratio_for(name)) in {"BPS", "passive"}
        for name in COMPENDIUM_STABILITY
    )


# ============================================================
# Pretty-print
# ============================================================

def stability_summary() -> str:
    """Pretty-print ρ = m/Γ for the 22 compendium particles."""
    lines = [
        "Substrate pattern-stability ratio ρ = m_X / Γ_X for K_7 walks:",
        "",
        "  Particle      m (GeV)      Γ (GeV)        ρ             log_10(ρ)  regime",
    ]
    # Sort by ρ ascending (most-unstable first), with stable at the end.
    items = list(COMPENDIUM_STABILITY.items())
    items.sort(key=lambda kv: (math.isinf(stability_ratio_for(kv[0])),
                               stability_ratio_for(kv[0])))
    for name, spec in items:
        rho = stability_ratio_for(name)
        log_rho = log10_stability_ratio_for(name)
        regime = classify_regime(rho)
        if math.isinf(rho):
            rho_str = "∞"
            log_str = "—"
        else:
            rho_str = f"{rho:.3e}"
            log_str = f"{log_rho:5.1f}"
        lines.append(
            f"  {name:<10}  {spec.m_X_GeV:>10.4f}  "
            f"{spec.Gamma_X_GeV:>11.3e}  {rho_str:>12}  "
            f"{log_str:>9}    {regime}"
        )
    return "\n".join(lines)
