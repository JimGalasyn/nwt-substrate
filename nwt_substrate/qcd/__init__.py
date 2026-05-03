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
]
