"""
QCD process objects for the nwt.qcd shim.

Each process exposes the uniform quartet:
    .M_squared_avg(...)       spin- and color-averaged |M|^2
    .dsigma_dOmega(...)        differential cross-section in pb/sr
    .sigma_total(...)          total cross-section in pb (with optional cutoff)
    .event(...)                full kinematic + amplitude object

Plus  .diagrams.<channel>  with named Diagram objects for each channel.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Optional

import numpy as np

from ..algebra.octonions import make_octonion_table
from ..algebra.dirac import make_lorentz_gammas
from ..amplitudes import processes as _proc
from ..amplitudes import cross_sections as _xs
from ..qed.diagram import Diagram as _Diagram
from . import diagram as _diag_mod


_T_OCT = make_octonion_table()
_GAMMAS = make_lorentz_gammas(_T_OCT)
GEV_INV_SQ_TO_PB = 3.894e8


@dataclass
class QCDEvent:
    """A QCD scattering event (analog of the QED Event)."""
    process: str
    s: float
    t: float
    u: float
    M_sq_avg: float
    dsigma_dOmega_pb: float
    sigma_total_pb: Optional[float] = None
    extra: dict = None

    def __repr__(self):
        return (f"QCDEvent({self.process!r}, sqrt(s)={np.sqrt(self.s):.3f} GeV, "
                f"|M|^2={self.M_sq_avg:.3e}, "
                f"dsigma/dOmega={self.dsigma_dOmega_pb:.3e} pb/sr)")


def _mandelstam(p1, p2, p3, p4):
    def dot(a, b):
        return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
    return (float(dot(p1+p2, p1+p2)),
            float(dot(p1-p3, p1-p3)),
            float(dot(p1-p4, p1-p4)))


# ---------------------------------------------------------------------------
# q qbar -> q' qbar' (different flavor, s-channel only)
# ---------------------------------------------------------------------------

class _QQbarProcess:
    name = "qqbar"

    def M_squared_avg(self, E_cm: float, theta: float,
                      m_q: float = 2.16e-3, m_qprime: float = 4.67e-3,
                      alpha_s: float = 0.1179) -> float:
        M, *_ = _proc.qcd_qqbar.M_squared_avg(
            E_cm, theta, m_q, m_qprime, _GAMMAS, alpha_s)
        return M

    def dsigma_dOmega(self, E_cm: float, theta: float,
                      m_q: float = 2.16e-3, m_qprime: float = 4.67e-3,
                      alpha_s: float = 0.1179, units: str = "pb") -> float:
        M_sq = self.M_squared_avg(E_cm, theta, m_q, m_qprime, alpha_s)
        s = E_cm ** 2
        ds = _xs.dsigma_dOmega_cm(M_sq, s, (m_q, m_q), (m_qprime, m_qprime))
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total(self, E_cm: float,
                    m_q: float = 2.16e-3, m_qprime: float = 4.67e-3,
                    alpha_s: float = 0.1179, units: str = "pb") -> float:
        s = E_cm ** 2
        sigma = _xs.sigma_2to2_cm(
            lambda th: self.M_squared_avg(E_cm, th, m_q, m_qprime, alpha_s),
            s, (m_q, m_q), (m_qprime, m_qprime),
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    def event(self, E_cm: float, theta: float,
              m_q: float = 2.16e-3, m_qprime: float = 4.67e-3,
              alpha_s: float = 0.1179) -> QCDEvent:
        M_sq, p1, p2, p3, p4 = _proc.qcd_qqbar.M_squared_avg(
            E_cm, theta, m_q, m_qprime, _GAMMAS, alpha_s)
        s, t, u = _mandelstam(p1, p2, p3, p4)
        return QCDEvent(
            process=self.name, s=s, t=t, u=u,
            M_sq_avg=M_sq,
            dsigma_dOmega_pb=self.dsigma_dOmega(E_cm, theta, m_q, m_qprime, alpha_s),
            sigma_total_pb=self.sigma_total(E_cm, m_q, m_qprime, alpha_s),
            extra={"E_cm": E_cm, "theta": theta, "alpha_s": alpha_s},
        )

    def __init__(self):
        self.diagrams = SimpleNamespace(
            s_channel=_Diagram(
                process_name="qqbar", channel="s_channel",
                expression=_diag_mod.EXPRESSIONS["qqbar"],
                feynman_rules={"propagator": "gluon",
                               "vertex": "QCD quark-gluon",
                               "color_factor_avg": "2/9"},
                _render_fn=_diag_mod.render_qqbar_s,
                _tikz_template=_diag_mod.TIKZ_QQBAR_S,
            ),
        )


# ---------------------------------------------------------------------------
# q q -> q q (same flavor, t + u channels with substrate doubling)
# ---------------------------------------------------------------------------

class _QQProcess:
    name = "qq"

    def M_squared_avg(self, E_cm: float, theta: float,
                      m_q: float = 2.16e-3, alpha_s: float = 0.1179) -> float:
        M, *_ = _proc.qcd_qq.M_squared_avg(E_cm, theta, m_q, _GAMMAS, alpha_s)
        return M

    def dsigma_dOmega(self, E_cm: float, theta: float,
                      m_q: float = 2.16e-3, alpha_s: float = 0.1179,
                      units: str = "pb") -> float:
        M_sq = self.M_squared_avg(E_cm, theta, m_q, alpha_s)
        s = E_cm ** 2
        ds = _xs.dsigma_dOmega_cm(M_sq, s, (m_q, m_q), (m_q, m_q))
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total(self, E_cm: float, cutoff_deg: float = 30.0,
                    m_q: float = 2.16e-3, alpha_s: float = 0.1179,
                    units: str = "pb") -> float:
        s = E_cm ** 2
        tmin = cutoff_deg * np.pi / 180.0
        tmax = np.pi - tmin
        sigma = _xs.sigma_2to2_cm(
            lambda th: self.M_squared_avg(E_cm, th, m_q, alpha_s),
            s, (m_q, m_q), (m_q, m_q),
            theta_min=tmin, theta_max=tmax,
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    def __init__(self):
        self.diagrams = SimpleNamespace(
            t_channel=_Diagram(
                process_name="qq", channel="t_channel",
                expression=_diag_mod.EXPRESSIONS["qq"],
                feynman_rules={"propagator": "gluon",
                               "color_factor": "2/9 (diagonal)",
                               "fermi_sign": "negative (M_t - M_u)"},
                _render_fn=_diag_mod.render_qq_t,
                _tikz_template=_diag_mod.TIKZ_QQ_T,
            ),
            u_channel=_Diagram(
                process_name="qq", channel="u_channel",
                expression=_diag_mod.EXPRESSIONS["qq"],
                feynman_rules={"propagator": "gluon",
                               "color_factor": "2/9 (diagonal)",
                               "color_cross": "-2/27 (interference)"},
                _render_fn=_diag_mod.render_qq_u,
                _tikz_template=_diag_mod.TIKZ_QQ_U,
            ),
        )


# ---------------------------------------------------------------------------
# gg -> gg (4 diagrams: s + t + u + 4-gluon)
# ---------------------------------------------------------------------------

class _GGProcess:
    name = "gg"

    def textbook_M_squared(self, E_cm: float, theta: float,
                           alpha_s: float = 0.1179) -> float:
        """Textbook (9/2) g_s^4 [3 - ut/s^2 - us/t^2 - st/u^2]."""
        return _proc.qcd_gg.textbook_M_squared(E_cm, theta, alpha_s)

    def dsigma_dOmega_textbook(self, E_cm: float, theta: float,
                                alpha_s: float = 0.1179,
                                units: str = "pb") -> float:
        """Textbook differential cross-section using the closed-form |M|^2."""
        M_sq = self.textbook_M_squared(E_cm, theta, alpha_s)
        s = E_cm ** 2
        # gg -> gg: massless initial and final
        ds = _xs.dsigma_dOmega_cm(M_sq, s, (0.0, 0.0), (0.0, 0.0))
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total_textbook(self, E_cm: float, cutoff_deg: float = 30.0,
                              alpha_s: float = 0.1179,
                              units: str = "pb") -> float:
        """Textbook integrated cross-section with cutoff (avoid t,u poles)."""
        s = E_cm ** 2
        tmin = cutoff_deg * np.pi / 180.0
        tmax = np.pi - tmin
        sigma = _xs.sigma_2to2_cm(
            lambda th: self.textbook_M_squared(E_cm, th, alpha_s),
            s, (0.0, 0.0), (0.0, 0.0),
            theta_min=tmin, theta_max=tmax,
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    def __init__(self):
        self.diagrams = SimpleNamespace(
            s_channel=_Diagram(
                process_name="gg", channel="s_channel",
                expression=_diag_mod.EXPRESSIONS["gg"],
                feynman_rules={"vertex": "two 3-gluon",
                               "propagator": "gluon",
                               "color_factor": "f^{abe} f^{cde}"},
                _render_fn=_diag_mod.render_gg_s,
                _tikz_template=_diag_mod.TIKZ_GG_S,
            ),
            four_gluon=_Diagram(
                process_name="gg", channel="four_gluon",
                expression=_diag_mod.EXPRESSIONS["gg"],
                feynman_rules={"vertex": "4-gluon contact",
                               "propagator": "(none, contact term)",
                               "color_structures": "3 distinct (f-products)"},
                _render_fn=_diag_mod.render_gg_4vertex,
                _tikz_template=_diag_mod.TIKZ_GG_4VERTEX,
            ),
        )


# ---------------------------------------------------------------------------
# Singletons
# ---------------------------------------------------------------------------

qqbar = _QQbarProcess()
qq = _QQProcess()
gg = _GGProcess()
