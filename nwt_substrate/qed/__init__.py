"""
nwt.qed
=======

Textbook-style QED interface over the substrate algebra.

Designed to feel familiar to anyone reading Peskin & Schroeder or Schwartz —
expose the physics primitives (spinors, propagators, vertex factors) by their
canonical letters, expose the standard processes (Compton, e+e- -> mu+mu-,
Møller, Bhabha) with `M_squared_avg / dsigma_dOmega / sigma_total / event`
methods, and provide a Feynman-diagram interface (matplotlib + TikZ-Feynman)
for both pedagogy and paper inclusion.

Underneath, every call routes through the substrate primitives in
`nwt_substrate.amplitudes`, `nwt_substrate.algebra`, and friends.

Quick start::

    import nwt_substrate.qed as qed

    # Constants
    qed.alpha          # 1/137.036
    qed.alpha_at(91)   # running alpha at M_Z

    # A Compton scattering result
    event = qed.compton.event(omega=0.01, theta=1.0)
    print(event)
    # Event('compton', sqrt(s)=...)

    # Total cross-section for e+e- -> mu+mu- at 10 GeV
    sigma_pb = qed.eemumu.sigma_total(E_cm=10.0)   # ~870 pb

    # Render a diagram + emit TikZ-Feynman
    fig = qed.compton.diagrams.s_channel.render()
    tikz_code = qed.compton.diagrams.s_channel.to_tikz()

    # Multi-panel gallery for Paper 18
    fig = qed.gallery_all()
"""

from __future__ import annotations

import matplotlib.pyplot as plt

# ---- Constants ----
from .constants import (
    alpha, e_charge,
    m_e, m_mu, m_tau,
    m_Z, m_W,
    r_e,
)

# ---- Substrate-ISA identities ----
from ..isa import (
    DIM_OCTONION,
    DIM_S_SPIN7,
    DIM_V_SPIN7,
    N_CL07_IMAGINARY,
    N_LORENTZ_FROM_CL07,
    DIM_INTERNAL_SU2,
    N_VERTICES_K7,
    N_EDGES_K7,
)


def substrate_breakdown() -> str:
    """Pretty-print the substrate Cl(0,7)→Cl(1,3) construction.

    Spells out the structural identities that QED inherits from the
    substrate algebra. Every magic number below is sourced from
    `nwt_substrate.isa.constants` — a refactor that violates them
    breaks tests across all shims.
    """
    lines = ["QED from substrate Cl(0,7) octonion algebra:"]
    lines.append("")
    lines.append(f"    Cl(0,7) has {N_CL07_IMAGINARY} imaginary generators "
                 f"(= N_CL07_IMAGINARY = DIM_V_SPIN7 = |V(K_7)|).")
    lines.append(f"    The octonion algebra is {DIM_OCTONION}-dim "
                 f"(= DIM_OCTONION = DIM_S_SPIN7 = Spin(7) spinor rep).")
    lines.append("")
    lines.append(f"    Lorentz Cl(1,3) embedding:")
    lines.append(f"      pick {N_LORENTZ_FROM_CL07} of {N_CL07_IMAGINARY} imaginary directions "
                 f"({N_LORENTZ_FROM_CL07} = N_LORENTZ_FROM_CL07)")
    lines.append(f"      Wick-rotate one → time direction")
    lines.append(f"      → 4 Dirac γ^μ as {DIM_OCTONION}×{DIM_OCTONION} matrices")
    lines.append("")
    lines.append(f"    Remaining {DIM_INTERNAL_SU2} directions "
                 f"({N_CL07_IMAGINARY} - {N_LORENTZ_FROM_CL07} = {DIM_INTERNAL_SU2}) "
                 f"generate internal SU(2)")
    lines.append(f"    that commutes with all γ^μ. "
                 f"(= DIM_INTERNAL_SU2 = RANK_SO7 = 3)")
    lines.append("")
    lines.append(f"    Sanity: 4 (Lorentz) + 3 (internal) = 7 (Cl(0,7) imaginary)")
    lines.append(f"                                       = |V(K_7)| = DIM_V_SPIN7")
    lines.append("")
    lines.append(f"    Cross-shim ID: the {N_EDGES_K7} edges of K_7 (= so(7) generators)")
    lines.append(f"      drive: gravity's α^({N_EDGES_K7}/2) = α^10.5 Wilson exponent")
    lines.append(f"             chemistry's K_7-hub Tr(M²) ≤ -24 indicator")
    lines.append(f"             QED's 8×8 γ^μ embedded in the same Spin(7) "
                 f"adjoint")
    return "\n".join(lines)


# Expose substrate identities as a namespace too (for `qed.substrate.X` access)
class _SubstrateNamespace:
    """Cross-shim K_7 substrate identities visible from the QED shim.

    All values are sourced from `nwt_substrate.isa.constants`. If one
    of these doesn't match `isa.X`, the cross-shim tests fail.
    """
    DIM_OCTONION = DIM_OCTONION
    DIM_S_SPIN7 = DIM_S_SPIN7
    DIM_V_SPIN7 = DIM_V_SPIN7
    N_CL07_IMAGINARY = N_CL07_IMAGINARY
    N_LORENTZ_FROM_CL07 = N_LORENTZ_FROM_CL07
    DIM_INTERNAL_SU2 = DIM_INTERNAL_SU2
    N_EDGES_K7 = N_EDGES_K7
    N_VERTICES_K7 = N_VERTICES_K7
    breakdown = staticmethod(substrate_breakdown)


substrate = _SubstrateNamespace()

# ---- Running coupling ----
from ..amplitudes.running_couplings import alpha_qed as _alpha_qed_running


def alpha_at(mu: float) -> float:
    """1-loop QED running coupling at scale mu (GeV)."""
    return _alpha_qed_running(mu)


# ---- Building blocks (named with their textbook letters) ----
from ..amplitudes import (
    spinors as spinor,
    polarizations as polarization,
    propagators as propagator,
)
from ..amplitudes.vertices import qed_vertex as vertex


# ---- Mandelstam helper ----
def mandelstam(p1, p2, p3, p4):
    """Compute (s, t, u) from four 4-momenta (mostly-minus metric)."""
    def dot(a, b):
        return float(a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3])
    s = dot(p1+p2, p1+p2)
    t = dot(p1-p3, p1-p3)
    u = dot(p1-p4, p1-p4)
    return s, t, u


# ---- Processes (compton, eemumu, moller, bhabha, muon_decay) ----
from .process import (
    compton, eemumu, moller, bhabha, muon_decay,
    Event,
)
from . import diagram as _diag_mod
Diagram = _diag_mod.Diagram


# ---- One-line gallery for paper figures ----
def gallery_all(figsize=(14, 10)):
    """
    Render every canonical QED process diagram in a single multi-panel
    figure.  Returns a matplotlib Figure suitable for inclusion in Paper 18
    or pedagogical materials.
    """
    fig, axes = plt.subplots(3, 2, figsize=figsize)
    compton.diagrams.s_channel.render(ax=axes[0, 0])
    axes[0, 0].set_title("Compton (s-channel)", fontsize=11)
    compton.diagrams.u_channel.render(ax=axes[0, 1])
    axes[0, 1].set_title("Compton (u-channel)", fontsize=11)
    eemumu.diagrams.s_channel.render(ax=axes[1, 0])
    axes[1, 0].set_title(r"$e^{+}e^{-} \to \mu^{+}\mu^{-}$", fontsize=11)
    muon_decay.diagrams.tree.render(ax=axes[1, 1])
    axes[1, 1].set_title(r"$\mu^{-}$ decay (V-A)", fontsize=11)
    moller.diagrams.t_channel.render(ax=axes[2, 0])
    axes[2, 0].set_title("Møller (t-channel)", fontsize=11)
    bhabha.diagrams.s_channel.render(ax=axes[2, 1])
    axes[2, 1].set_title("Bhabha (s-channel)", fontsize=11)
    fig.suptitle("Substrate-algebra QED tree diagrams", fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0.02, 1, 0.99])
    return fig


__all__ = [
    # constants
    "alpha", "e_charge", "m_e", "m_mu", "m_tau", "m_Z", "m_W", "r_e",
    "alpha_at",
    # building blocks
    "spinor", "polarization", "propagator", "vertex",
    "mandelstam",
    # processes
    "compton", "eemumu", "moller", "bhabha", "muon_decay",
    "Event", "Diagram",
    # gallery
    "gallery_all",
    # substrate identities
    "substrate", "substrate_breakdown",
]
