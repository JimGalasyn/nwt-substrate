"""
Substrate light-pseudoscalar decay constants f_K, f_η, f_π.

P7b §2-3 closure (Galasyn 2026-05-23,
analysis/paper21_p7b_decay_amplitudes.md in null-worldtube-private).

Closes the light pseudoscalar sector with two formulas covering distinct
substrate regimes:

(a) Cabibbo-scale law (K, η — ordinary light pseudoscalars)::

        f_X² / m_X² = 7 α              (= λ² Cabibbo-Wilson amplitude)
        f_X        = m_X · √(7α)

(b) Fibonacci-anomaly law (π — Goldstone chiral special case)::

        f_π² / m_π² = 1 / √F_5 = 1/√5
        f_π        = m_π / 5^(1/4)

The pion's anomalously-light Goldstone status (Gell-Mann-Oakes-Renner)
means its f_π²/m_π² ratio doesn't carry α-suppression; instead the
Fibonacci-Goldstone walk-length F_5 = 5 sets the cinquefoil-carrier scale.

Precision vs PDG 2024 (substrate α and PDG m_X anchor):

    f_K  = m_K · √(7α)       111.5 MeV   PDG 110.0   +1.4 %
    f_η  = m_η · √(7α)       123.8 MeV   PDG 130.0   −4.7 %
    f_π  = m_π / 5^(1/4)      93.4 MeV   PDG  92.4   +1.1 %

The 7α structural factor recurs in V_us² (Wolfenstein λ²), v_EW NLO
(25α/(4√3)), Sirlin Δq (α/(2π)(25/4 - π²)), heavy-meson f_X
strangeness factor √(7/4), vector-meson binding scale 7α m_τ², and
form-factor cos^N(θ_C) substrate Cabibbo — universal vocabulary.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

from .substrate_gf import ALPHA_SUBSTRATE
from .heavy_meson_decay_constants import (
    pi_zero_mass, f_pi_substrate)


# ============================================================
# Catalog
# ============================================================

@dataclass(frozen=True)
class LightPseudoscalarSpec:
    """Substrate spec for a light pseudoscalar meson.

    Fields:
      name      — meson label (e.g. "K", "eta", "pi")
      m_X_GeV   — PDG meson mass (GeV)
      f_PDG_MeV — PDG decay constant (MeV)
      regime    — "cabibbo" (Tier-K/η) or "fibonacci" (Tier-π)
    """
    name: str
    m_X_GeV: float
    f_PDG_MeV: float
    regime: str


LIGHT_PSEUDOSCALARS: Mapping[str, LightPseudoscalarSpec] = {
    "pi":  LightPseudoscalarSpec("pi",  0.13957, 92.4,  "fibonacci"),
    "K":   LightPseudoscalarSpec("K",   0.49368, 110.0, "cabibbo"),
    "eta": LightPseudoscalarSpec("eta", 0.54786, 130.0, "cabibbo"),
}


# ============================================================
# Substrate formulas
# ============================================================

def cabibbo_scale_fX(m_X_GeV: float,
                     alpha: float = ALPHA_SUBSTRATE) -> float:
    """Substrate decay constant via Cabibbo-scale law: f = m · √(7α) [GeV].

    Applies to ordinary light pseudoscalars (K, η).  The 7α = λ² is the
    Cabibbo-Wilson amplitude on K_7 (Paper 6b: λ² = 7α).
    """
    return m_X_GeV * math.sqrt(7.0 * alpha)


def fibonacci_anomaly_fX(m_X_GeV: float,
                         walk_length: int = 5) -> float:
    """Substrate decay constant via Fibonacci-anomaly law [GeV].

    Applies to the pion (Goldstone chiral special case): the cinquefoil
    walk-length F_5 = 5 sets the substrate Goldstone scale::

        f_π = m_π / F_walk_length^(1/4) = m_π / 5^(1/4)

    Default walk_length = 5 = F_5 = pion T(2,5) cinquefoil-carrier index.
    """
    return m_X_GeV / math.pow(walk_length, 0.25)


def light_meson_fX_for(name: str,
                       alpha: float = ALPHA_SUBSTRATE) -> float:
    """Substrate f_X for a named light pseudoscalar (K, eta, pi) [GeV]."""
    if name not in LIGHT_PSEUDOSCALARS:
        raise KeyError(f"Unknown light pseudoscalar: {name}.  "
                       f"Known: {sorted(LIGHT_PSEUDOSCALARS)}")
    spec = LIGHT_PSEUDOSCALARS[name]
    if spec.regime == "cabibbo":
        return cabibbo_scale_fX(spec.m_X_GeV, alpha)
    elif spec.regime == "fibonacci":
        return fibonacci_anomaly_fX(spec.m_X_GeV, walk_length=5)
    raise ValueError(f"Unknown regime: {spec.regime}")


# ============================================================
# Precision chain
# ============================================================

def precision_chain(alpha: float = ALPHA_SUBSTRATE) -> dict:
    """Per-meson substrate-vs-PDG comparison for light pseudoscalars."""
    out = {}
    for name, spec in LIGHT_PSEUDOSCALARS.items():
        f_sub_GeV = light_meson_fX_for(name, alpha)
        f_sub_MeV = f_sub_GeV * 1e3
        out[name] = {
            'substrate_MeV': f_sub_MeV,
            'pdg_MeV':       spec.f_PDG_MeV,
            'rel_err':       f_sub_MeV / spec.f_PDG_MeV - 1.0,
            'regime':        spec.regime,
        }
    return out


def verify_light_meson_fX(alpha: float = ALPHA_SUBSTRATE,
                          tol: float = 0.05) -> dict:
    """Verify K, η, π substrate f_X match PDG within ``tol`` (default 5 %)."""
    chain = precision_chain(alpha)
    ok = all(abs(chain[k]['rel_err']) <= tol for k in LIGHT_PSEUDOSCALARS)
    chain['pass'] = ok
    return chain


def precision_chain_summary(alpha: float = ALPHA_SUBSTRATE) -> str:
    """Pretty-print substrate light-pseudoscalar decay-constant precision."""
    chain = precision_chain(alpha)
    lines = [
        "Substrate light-pseudoscalar decay constants:",
        "",
        "  Cabibbo-scale law (K, η):     f_X = m_X · √(7α)",
        "  Fibonacci anomaly (π):        f_π = m_π / 5^(1/4)",
        "",
    ]
    for name, spec in LIGHT_PSEUDOSCALARS.items():
        c = chain[name]
        formula = ("m · √(7α)" if spec.regime == "cabibbo"
                   else "m / 5^(1/4)")
        lines.append(
            f"  f_{name:<5} ({spec.regime:>9})  {formula:<13}  "
            f"= {c['substrate_MeV']:5.1f} MeV  "
            f"vs PDG {c['pdg_MeV']:5.1f}  "
            f"→ {c['rel_err'] * 100:+6.2f} %"
        )
    lines.append("")
    lines.append("  Free params:  0  (α + integer walk-length)")
    return "\n".join(lines)
