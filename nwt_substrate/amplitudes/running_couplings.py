"""
Running couplings via 1-loop renormalization-group equations.

QED 1-loop beta function (with thresholds):

    d alpha / d ln(mu) = (alpha^2 / (3 pi)) * b_QED(mu)

    b_QED(mu) = sum over active fermions of  N_c(f) * Q_f^2

For SM matter all active:
    leptons (3 x charge 1, no color)         contribute  3
    up-type quarks (3 x Q=2/3, N_c=3)        contribute  3 * 3 * 4/9 = 4
    down-type quarks (3 x Q=-1/3, N_c=3)     contribute  3 * 3 * 1/9 = 1
    total b_QED = 8

QCD 1-loop beta function:

    d alpha_s / d ln(mu) = -(alpha_s^2 / (2 pi)) * (11 - 2 n_f / 3)

    where n_f = number of active quark flavors with m_q < mu.

Solution (1-loop, single threshold sector):

    alpha_s(mu) = alpha_s(mu_0) / (1 + (alpha_s(mu_0)/(2 pi)) * b_0 * ln(mu/mu_0))

Standard PDG values to verify against:

    alpha_QED(M_Z) = 1/127.952     (full PDG, includes hadronic + higher-loop)
    alpha_s(M_Z)   = 0.1179 +- 0.001

We reproduce these to leading-log accuracy with constituent thresholds.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


# CODATA 2018
ALPHA_THOMSON = 1.0 / 137.035999084     # alpha_QED at q^2 -> 0

# Standard mass references (GeV)
M_E = 0.510998928e-3
M_MU = 0.10566
M_TAU = 1.77686
M_Z = 91.1876
M_W = 80.379
M_HIGGS = 125.10
M_TOP = 172.7


# Constituent / effective quark masses for QED running
# Using current m_u, m_d shows the "naive" running; for hadronic corrections
# one uses dispersion relations.  Constituent ~0.3 GeV gives a closer match
# at leading-log level.
QUARK_MASSES = {
    "u": 2.16e-3,        # current mass
    "d": 4.67e-3,
    "s": 0.0934,
    "c": 1.27,
    "b": 4.18,
    "t": M_TOP,
}

QUARK_MASSES_CONSTITUENT = {
    "u": 0.336,
    "d": 0.336,
    "s": 0.486,
    "c": 1.55,
    "b": 4.73,
    "t": M_TOP,
}


@dataclass
class FermionThreshold:
    name: str
    mass: float            # threshold mass (GeV)
    qed_b: float           # contribution to b_QED = N_c * Q^2
    n_f_qcd: int           # contribution to QCD n_f (1 for quark, 0 for lepton)


def standard_thresholds(use_constituent: bool = False) -> list[FermionThreshold]:
    """SM fermion thresholds in increasing mass order."""
    qm = QUARK_MASSES_CONSTITUENT if use_constituent else QUARK_MASSES
    thresholds = [
        FermionThreshold("e",   M_E,    1.0,      0),
        FermionThreshold("u",   qm["u"], 4.0/3.0, 1),
        FermionThreshold("d",   qm["d"], 1.0/3.0, 1),
        FermionThreshold("s",   qm["s"], 1.0/3.0, 1),
        FermionThreshold("mu",  M_MU,   1.0,      0),
        FermionThreshold("c",   qm["c"], 4.0/3.0, 1),
        FermionThreshold("tau", M_TAU, 1.0,       0),
        FermionThreshold("b",   qm["b"], 1.0/3.0, 1),
        FermionThreshold("t",   qm["t"], 4.0/3.0, 1),
    ]
    return sorted(thresholds, key=lambda x: x.mass)


# ---------------------------------------------------------------------------
# QED running coupling
# ---------------------------------------------------------------------------

def alpha_qed(mu: float,
              alpha_at_mu0: float = ALPHA_THOMSON,
              mu0: float = M_E,
              thresholds: Iterable[FermionThreshold] | None = None) -> float:
    """
    1-loop QED running coupling at scale mu.

    Starts from alpha_at_mu0 at the lightest threshold (electron) and
    integrates upward through each active fermion mass threshold.

    1/alpha(mu) = 1/alpha_0 - (1/(3 pi)) sum_{intervals} b_QED * ln(mu_2/mu_1)
    """
    if thresholds is None:
        thresholds = standard_thresholds()
    thresholds = sorted(thresholds, key=lambda x: x.mass)

    if mu <= mu0:
        return alpha_at_mu0

    one_over_alpha = 1.0 / alpha_at_mu0
    b_running = 0.0
    prev_scale = mu0

    for th in thresholds:
        if th.mass <= mu0:
            # Already counted as part of alpha_at_mu0 starting point;
            # but if mu0 > th.mass, threshold is below starting scale -> skip
            b_running += th.qed_b
            continue
        if th.mass >= mu:
            break
        # Run from prev_scale up to th.mass with current b_running
        if th.mass > prev_scale:
            one_over_alpha -= (b_running / (3.0 * np.pi)) * np.log(th.mass / prev_scale)
        b_running += th.qed_b
        prev_scale = th.mass

    # Final segment from prev_scale to mu
    if mu > prev_scale and b_running > 0:
        one_over_alpha -= (b_running / (3.0 * np.pi)) * np.log(mu / prev_scale)

    return 1.0 / one_over_alpha


# ---------------------------------------------------------------------------
# QCD running coupling
# ---------------------------------------------------------------------------

def alpha_s(mu: float,
            alpha_s_at_mz: float = 0.1179,
            thresholds: Iterable[FermionThreshold] | None = None) -> float:
    """
    1-loop QCD running coupling at scale mu, anchored at M_Z.

      alpha_s(mu) = alpha_s(mu_0) / [1 + (alpha_s(mu_0)/(2pi)) * b_0 * ln(mu/mu_0)]
      b_0 = 11 - 2 n_f / 3

    Uses standard quark thresholds; n_f changes when crossing m_b, m_c, etc.
    """
    if thresholds is None:
        thresholds = standard_thresholds()
    thresholds = sorted(thresholds, key=lambda x: x.mass)

    # Quark-only thresholds (QCD only)
    quark_th = [t for t in thresholds if t.n_f_qcd == 1]

    def b0(n_f: int) -> float:
        return 11.0 - 2.0 * n_f / 3.0

    # Determine n_f at M_Z (count quarks below M_Z)
    n_f_mz = sum(1 for t in quark_th if t.mass < M_Z)

    # Run from M_Z either up or down
    alpha = alpha_s_at_mz
    current = M_Z
    n_f = n_f_mz

    if mu == M_Z:
        return alpha

    if mu > M_Z:
        # Run upward through any thresholds above M_Z (just top usually)
        upper_thresholds = sorted([t for t in quark_th if M_Z < t.mass < mu],
                                  key=lambda x: x.mass)
        for t in upper_thresholds:
            if t.mass > current:
                alpha = alpha / (1.0 + (alpha / (2.0 * np.pi)) * b0(n_f)
                                  * np.log(t.mass / current))
                current = t.mass
                n_f += 1
        if mu > current:
            alpha = alpha / (1.0 + (alpha / (2.0 * np.pi)) * b0(n_f)
                              * np.log(mu / current))
        return alpha
    else:
        # Run downward
        lower_thresholds = sorted([t for t in quark_th if mu < t.mass < M_Z],
                                  key=lambda x: -x.mass)
        for t in lower_thresholds:
            if t.mass < current:
                alpha = alpha / (1.0 + (alpha / (2.0 * np.pi)) * b0(n_f)
                                  * np.log(t.mass / current))
                current = t.mass
                n_f -= 1
        if mu < current:
            alpha = alpha / (1.0 + (alpha / (2.0 * np.pi)) * b0(n_f)
                              * np.log(mu / current))
        return alpha


# ---------------------------------------------------------------------------
# Reference table
# ---------------------------------------------------------------------------

def pdg_alpha_qed_at_mz() -> float:
    """PDG value 1/alpha_QED(M_Z) = 127.952 -> alpha = 0.0078126."""
    return 1.0 / 127.952


def pdg_alpha_s_at_mz() -> float:
    """PDG value alpha_s(M_Z) = 0.1179."""
    return 0.1179


def lambda_qcd(alpha_s_at_mz: float = 0.1179, n_f: int = 5) -> float:
    """
    Estimate Lambda_QCD from alpha_s(M_Z).
    1-loop:  Lambda = M_Z * exp(-2 pi / (b_0 alpha_s(M_Z)))
    """
    b_0 = 11.0 - 2.0 * n_f / 3.0
    return M_Z * np.exp(-2.0 * np.pi / (b_0 * alpha_s_at_mz))
