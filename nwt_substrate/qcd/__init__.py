"""
nwt.qcd
=======

Textbook-style QCD interface over the substrate algebra.

Designed to feel familiar to anyone reading Peskin & Schroeder, Schwartz,
or Halzen & Martin — exposes coupling constants, Gell-Mann matrices and
color algebra, and the standard scattering processes (q qbar -> q' qbar',
q q -> q q, gg -> gg) with the same `M_squared_avg / dsigma_dOmega /
sigma_total / event` quartet that the QED shim uses.

Quick start::

    import nwt_substrate.qcd as qcd

    # Constants
    qcd.alpha_s            # 0.1179 at M_Z
    qcd.alpha_s_at(91.2)   # running coupling
    qcd.C_F, qcd.C_A       # 4/3, 3
    qcd.Lambda_QCD         # 87 MeV (1-loop)

    # Color algebra
    T = qcd.T              # 8 SU(3) generators (lambda^a / 2)
    f = qcd.f              # structure constants

    # Processes (analogous to qed.eemumu, qed.compton)
    qcd.qqbar.sigma_total(E_cm=10.0)        # q qbar -> q' qbar' in pb
    qcd.qq.dsigma_dOmega(E_cm=10.0, theta=1.0)
    qcd.gg.textbook_M_squared(E_cm=100.0, theta=1.5)

    # Diagrams (gluon-as-coil rendering, TikZ-Feynman 'gluon' line type)
    fig = qcd.qqbar.diagrams.s_channel.render()
    tikz = qcd.qq.diagrams.t_channel.to_tikz()

    # Multi-panel gallery
    fig = qcd.gallery_all()
"""

from __future__ import annotations

import matplotlib.pyplot as plt

# ---- Constants ----
from .constants import (
    alpha_s, g_s,
    N_c, C_F, C_A, T_R,
    Lambda_QCD_5flavor as Lambda_QCD,
    m_u, m_d, m_s, m_c, m_b, m_t,
    m_u_constituent, m_d_constituent, m_s_constituent,
    m_Z, m_proton, Lambda_chiral,
)

# ---- Substrate-ISA identities ----
from ..isa import (
    N_C_SU3,
    N_GLUONS,
    C_F_SU3,
    C_A_SU3,
    T_R_SU3,
    DIM_OCTONION,
    DIM_S_SPIN7,
    RANK_SO7,
    N_EDGES_K7,
    N_VERTICES_K7,
)


def substrate_breakdown() -> str:
    """Pretty-print the SU(3)-as-substrate decomposition.

    SU(3) color sits inside G_2 = aut(octonions) ⊂ Spin(7), and the
    substrate identities below trace each QCD constant to the same
    K_7 / Spin(7) primitives that drive chemistry's K_7 indicator,
    gravity's α^(21/2), and QED's 8×8 γ^μ.
    """
    lines = ["QCD from substrate SU(3) ⊂ G_2 ⊂ Spin(7):"]
    lines.append("")
    lines.append(f"    N_c = N_C_SU3 = RANK_SO7 = {N_C_SU3}")
    lines.append(f"      (= 3 colors; same as so(7) rank; same as the")
    lines.append(f"       DIM_INTERNAL_SU2 directions left over in Cl(0,7)")
    lines.append(f"       after the {4}-of-{N_VERTICES_K7} Lorentz partition)")
    lines.append("")
    lines.append(f"    N_gluons = N_GLUONS = DIM_OCTONION = "
                 f"DIM_S_SPIN7 = {N_GLUONS}")
    lines.append(f"      (= N_c² - 1 = 8 gluons in the SU(3) adjoint)")
    lines.append(f"      (SUBSTRATE IDENTITY: same 8 that makes Dirac γ^μ")
    lines.append(f"       into {DIM_OCTONION}×{DIM_OCTONION} matrices)")
    lines.append("")
    lines.append(f"    C_F = DIM_OCTONION / (2 × N_C_SU3) = "
                 f"{DIM_OCTONION}/{2*N_C_SU3} = {C_F_SU3}")
    lines.append(f"      (= (N_c² - 1)/(2 N_c) = 4/3 fundamental Casimir)")
    lines.append("")
    lines.append(f"    C_A = N_C_SU3 = {C_A_SU3}")
    lines.append(f"      (= 3 adjoint Casimir = N_c)")
    lines.append("")
    lines.append(f"    T_R = 1/2")
    lines.append("")
    lines.append(f"    Cross-shim ID: the {DIM_OCTONION}-dim octonion algebra")
    lines.append(f"      gives:")
    lines.append(f"        QED:   {DIM_OCTONION}×{DIM_OCTONION} Dirac γ^μ matrices")
    lines.append(f"        QCD:   {N_GLUONS} gluons in SU(3) adjoint")
    from ..isa.constants import DIM_V_SPIN7
    lines.append(f"        Gravity:  the {DIM_S_SPIN7}/{DIM_V_SPIN7} = "
                 f"DIM_S/DIM_V prefactor in m_e/M_Pl")
    lines.append(f"      and the {N_EDGES_K7} edges of K_7 (= so(7) generators)")
    lines.append(f"      give:")
    lines.append(f"        Gravity:  α^({N_EDGES_K7}/2) K_7 Wilson exponent")
    lines.append(f"        Chemistry: K_7-hub indicator (W_6 wheel Tr(M²)=-24)")
    return "\n".join(lines)


class _SubstrateNamespace:
    """Cross-shim K_7 substrate identities visible from the QCD shim.

    All values are sourced from `nwt_substrate.isa.constants`.
    """
    N_C_SU3 = N_C_SU3
    N_GLUONS = N_GLUONS
    C_F_SU3 = C_F_SU3
    C_A_SU3 = C_A_SU3
    T_R_SU3 = T_R_SU3
    DIM_OCTONION = DIM_OCTONION
    RANK_SO7 = RANK_SO7
    N_EDGES_K7 = N_EDGES_K7
    breakdown = staticmethod(substrate_breakdown)


substrate = _SubstrateNamespace()

# ---- Running coupling ----
from ..amplitudes.running_couplings import alpha_s as _alpha_s_running


def alpha_s_at(mu: float, alpha_s_at_mz: float = 0.1179) -> float:
    """1-loop QCD running coupling at scale mu (GeV), anchored at M_Z."""
    return _alpha_s_running(mu, alpha_s_at_mz)


# ---- Color algebra primitives ----
from ..algebra.su3 import (
    gell_mann_matrices,
    su3_generators as _su3_gens,
    structure_constants as _f,
    d_constants as _d,
    fundamental_casimir,
    adjoint_casimir,
)


# Conventional aliases that QCD practitioners would expect:
def gell_mann():
    """8 Hermitian Gell-Mann matrices lambda^a (a=1..8)."""
    return gell_mann_matrices()


def T():
    """The 8 SU(3) fundamental generators T^a = lambda^a / 2."""
    return _su3_gens()


def f():
    """Structure constants f^abc as a (8,8,8) totally antisymmetric array."""
    return _f()


def d():
    """Symmetric tensor d^abc as a (8,8,8) totally symmetric array."""
    return _d()


# ---- Process objects ----
from .process import (
    qqbar, qq, gg,
    QCDEvent,
)
from . import diagram as _diag_mod
from ..qed.diagram import Diagram


# ---- Substrate confinement theorem (closed-walk on K_7) ----
# P5b closure (Galasyn 2026-05-23): color confinement is a STRUCTURAL
# CONSEQUENCE of the K_7 graph + SU(3) ⊂ G_2 ⊂ Spin(7) chain, not a
# dynamical mechanism.  Every closed K_7 walk is in the SU(3)-singlet
# sector under the canonical v̂ = e_4 partition (3 + 3̄ + 1).
from .confinement import (
    POLAR_VERTEX,
    TRIPLET_VERTICES,
    ANTITRIPLET_VERTICES,
    VERTEX_INDICES,
    color_charge,
    Walk,
    net_color,
    is_singlet_walk,
    enumerate_closed_walks,
    verify_confinement,
    verify_confinement_range,
    GLUEBALL_CANDIDATE_WALK,
    MESON_CANDIDATE_WALK,
    BARYON_CANDIDATE_WALK,
)


# ---- Substrate universal mass formula for hidden-gluon bound states ----
# P7b §7.8 closure (Galasyn 2026-05-23): glueballs, tetraquarks,
# pentaquarks, exotic Z_b all live on m² = (4 m_π⁰)² · N substrate ladder.
from .exotic_states import (
    pi_zero_mass,
    gluon_pair_scale,
    universal_mass,
    universal_N,
    BOUND_STATES_CATALOG,
    precision_chain as exotic_precision_chain,
    verify_universal_mass_formula,
    precision_chain_summary as exotic_precision_chain_summary,
)


# ---- Gallery ----
def gallery_all(figsize=(12, 8)):
    """
    Render every canonical QCD process diagram in a single multi-panel
    figure.  Returns a matplotlib Figure.
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    qqbar.diagrams.s_channel.render(ax=axes[0, 0])
    axes[0, 0].set_title(r"$q\bar{q} \to q'\bar{q}'$ (s-channel)", fontsize=11)
    qq.diagrams.t_channel.render(ax=axes[0, 1])
    axes[0, 1].set_title(r"$qq \to qq$ (t-channel)", fontsize=11)
    gg.diagrams.s_channel.render(ax=axes[1, 0])
    axes[1, 0].set_title(r"$gg \to gg$ (s-channel, two 3-gluon vertices)",
                          fontsize=11)
    gg.diagrams.four_gluon.render(ax=axes[1, 1])
    axes[1, 1].set_title(r"$gg \to gg$ (4-gluon contact)", fontsize=11)
    fig.suptitle("Substrate-algebra QCD tree diagrams", fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0.02, 1, 0.99])
    return fig


__all__ = [
    # constants
    "alpha_s", "g_s",
    "N_c", "C_F", "C_A", "T_R",
    "Lambda_QCD",
    "m_u", "m_d", "m_s", "m_c", "m_b", "m_t",
    "m_u_constituent", "m_d_constituent", "m_s_constituent",
    "m_Z", "m_proton", "Lambda_chiral",
    "alpha_s_at",
    # color algebra
    "gell_mann", "T", "f", "d",
    "fundamental_casimir", "adjoint_casimir",
    # processes
    "qqbar", "qq", "gg",
    "QCDEvent", "Diagram",
    # gallery
    "gallery_all",
    # substrate universal mass formula (exotic states)
    "pi_zero_mass", "gluon_pair_scale",
    "universal_mass", "universal_N",
    "BOUND_STATES_CATALOG",
    "exotic_precision_chain", "verify_universal_mass_formula",
    "exotic_precision_chain_summary",
    # substrate confinement theorem (P5b closure)
    "POLAR_VERTEX", "TRIPLET_VERTICES", "ANTITRIPLET_VERTICES",
    "VERTEX_INDICES", "color_charge",
    "Walk", "net_color", "is_singlet_walk",
    "enumerate_closed_walks", "verify_confinement",
    "verify_confinement_range",
    "GLUEBALL_CANDIDATE_WALK", "MESON_CANDIDATE_WALK",
    "BARYON_CANDIDATE_WALK",
    # substrate identities
    "substrate", "substrate_breakdown",
]
