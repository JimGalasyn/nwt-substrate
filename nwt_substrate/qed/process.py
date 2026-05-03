"""
Process objects for the QED shim — uniform interface across compton,
e+e- -> mu+mu-, Møller, Bhabha, muon decay.

Each `Process` exposes:
    .M_squared_avg(...)         spin/color-averaged |M|^2 in QED rel-norm
    .dsigma_dOmega(...)          differential cross-section in pb/sr
    .sigma_total(...)            total cross-section in pb
    .event(...)                  full kinematic + amplitude object (notebook helper)
    .diagrams                    namespace with named Diagram objects
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
from ..amplitudes import diagrams as _draw
from . import diagram as _diag_mod


# Cache the substrate gamma matrices once
_T_OCT = make_octonion_table()
_GAMMAS = make_lorentz_gammas(_T_OCT)


# Conversion: 1 GeV^-2 = 0.3894 mb = 3.894e8 pb
GEV_INV_SQ_TO_PB = 3.894e8


@dataclass
class Event:
    """A single phase-space point with all kinematic + amplitude data."""
    process: str
    s: float
    t: float
    u: float
    M_sq_avg: float
    dsigma_dOmega_pb: float
    sigma_total_pb: Optional[float] = None
    extra: dict = None

    def __repr__(self):
        return (f"Event({self.process!r}, sqrt(s)={np.sqrt(self.s):.3f} GeV, "
                f"|M|^2={self.M_sq_avg:.3e}, "
                f"dsigma/dOmega={self.dsigma_dOmega_pb:.3e} pb/sr)")


def _mandelstam_from_2to2(p1, p2, p3, p4):
    """Compute s, t, u from four 4-momenta (mostly-minus metric)."""
    def dot(a, b):
        return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
    s = float(dot(p1+p2, p1+p2))
    t = float(dot(p1-p3, p1-p3))
    u = float(dot(p1-p4, p1-p4))
    return s, t, u


# ---------------------------------------------------------------------------
# Compton  (lab frame: e- at rest, photon along +z)
# ---------------------------------------------------------------------------

class _ComptonProcess:
    name = "compton"

    def M_squared_avg(self, omega: float, theta: float,
                      m_e: float = 0.510998928e-3) -> float:
        M, _ = _proc.compton.M_squared_avg(omega, theta, m_e, _GAMMAS)
        return M

    def dsigma_dOmega(self, omega: float, theta: float,
                      m_e: float = 0.510998928e-3,
                      units: str = "pb") -> float:
        M_sq, omega_p = _proc.compton.M_squared_avg(omega, theta, m_e, _GAMMAS)
        # Lab frame: dsigma/dOmega = |M|^2 / (64 pi^2 m_e^2) * (omega'/omega)^2
        ds = M_sq * (omega_p / omega) ** 2 / (64.0 * np.pi ** 2 * m_e ** 2)
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total(self, omega: float, m_e: float = 0.510998928e-3,
                    units: str = "pb") -> float:
        sigma = _xs.sigma_compton_lab(
            lambda th: _proc.compton.M_squared_avg(omega, th, m_e, _GAMMAS),
            m_e,
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    @property
    def thomson_limit(self) -> float:
        """Thomson cross-section sigma = (8 pi / 3) r_e^2, in pb."""
        from . import constants as C
        sigma = _xs.thomson_cross_section(C.m_e, alpha=C.alpha)
        return sigma * GEV_INV_SQ_TO_PB

    def event(self, omega: float, theta: float,
              m_e: float = 0.510998928e-3) -> Event:
        p, k, p_p, k_p, omega_p = _proc.compton.kinematics(omega, theta, m_e)
        s, t, u = _mandelstam_from_2to2(p, k, p_p, k_p)
        M_sq, _ = _proc.compton.M_squared_avg(omega, theta, m_e, _GAMMAS)
        ds = self.dsigma_dOmega(omega, theta, m_e)
        return Event(
            process=self.name, s=s, t=t, u=u,
            M_sq_avg=M_sq, dsigma_dOmega_pb=ds,
            extra={"omega_in": omega, "omega_out": omega_p, "theta": theta},
        )

    def __init__(self):
        self.diagrams = SimpleNamespace(
            s_channel=_diag_mod.Diagram(
                process_name="compton", channel="s_channel",
                expression=_draw.EXPRESSIONS["compton"],
                feynman_rules={"propagator": "fermion", "vertex": "QED"},
                _render_fn=lambda ax: _draw.render_compton(ax, "s"),
                _tikz_template=_diag_mod.TIKZ_COMPTON_S,
            ),
            u_channel=_diag_mod.Diagram(
                process_name="compton", channel="u_channel",
                expression=_draw.EXPRESSIONS["compton"],
                feynman_rules={"propagator": "fermion", "vertex": "QED"},
                _render_fn=lambda ax: _draw.render_compton(ax, "u"),
                _tikz_template=_diag_mod.TIKZ_COMPTON_U,
            ),
        )


# ---------------------------------------------------------------------------
# e+ e- -> mu+ mu-
# ---------------------------------------------------------------------------

class _EEMuMuProcess:
    name = "eemumu"

    def M_squared_avg(self, E_cm: float, theta: float,
                      m_e: float = 0.510998928e-3,
                      m_mu: float = 0.10566) -> float:
        M, *_ = _proc.eemumu.M_squared_avg(E_cm, theta, m_e, m_mu, _GAMMAS)
        return M

    def dsigma_dOmega(self, E_cm: float, theta: float,
                      m_e: float = 0.510998928e-3,
                      m_mu: float = 0.10566,
                      units: str = "pb") -> float:
        M_sq = self.M_squared_avg(E_cm, theta, m_e, m_mu)
        s = E_cm ** 2
        ds = _xs.dsigma_dOmega_cm(M_sq, s, (m_e, m_e), (m_mu, m_mu))
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total(self, E_cm: float,
                    m_e: float = 0.510998928e-3,
                    m_mu: float = 0.10566,
                    units: str = "pb") -> float:
        s = E_cm ** 2
        sigma = _xs.sigma_2to2_cm(
            lambda th: self.M_squared_avg(E_cm, th, m_e, m_mu),
            s, (m_e, m_e), (m_mu, m_mu),
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    def event(self, E_cm: float, theta: float,
              m_e: float = 0.510998928e-3,
              m_mu: float = 0.10566) -> Event:
        M_sq, p1, p2, p3, p4 = _proc.eemumu.M_squared_avg(
            E_cm, theta, m_e, m_mu, _GAMMAS)
        s, t, u = _mandelstam_from_2to2(p1, p2, p3, p4)
        return Event(
            process=self.name, s=s, t=t, u=u,
            M_sq_avg=M_sq,
            dsigma_dOmega_pb=self.dsigma_dOmega(E_cm, theta, m_e, m_mu),
            sigma_total_pb=self.sigma_total(E_cm, m_e, m_mu),
            extra={"E_cm": E_cm, "theta": theta},
        )

    def __init__(self):
        self.diagrams = SimpleNamespace(
            s_channel=_diag_mod.Diagram(
                process_name="eemumu", channel="s_channel",
                expression=_draw.EXPRESSIONS["eemumu"],
                feynman_rules={"propagator": "photon", "vertex": "QED"},
                _render_fn=lambda ax: _draw.render_eemumu(ax),
                _tikz_template=_diag_mod.TIKZ_EEMUMU_S,
            ),
        )


# ---------------------------------------------------------------------------
# Møller and Bhabha
# ---------------------------------------------------------------------------

class _MollerProcess:
    name = "moller"

    def M_squared_avg(self, E_cm: float, theta: float,
                      m_e: float = 0.510998928e-3) -> float:
        M, *_ = _proc.moller.M_squared_avg(E_cm, theta, m_e, _GAMMAS)
        return M

    def dsigma_dOmega(self, E_cm: float, theta: float,
                      m_e: float = 0.510998928e-3,
                      units: str = "pb") -> float:
        M_sq = self.M_squared_avg(E_cm, theta, m_e)
        s = E_cm ** 2
        ds = _xs.dsigma_dOmega_cm(M_sq, s, (m_e, m_e), (m_e, m_e))
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total(self, E_cm: float, cutoff_deg: float = 30.0,
                    m_e: float = 0.510998928e-3,
                    units: str = "pb") -> float:
        s = E_cm ** 2
        tmin = cutoff_deg * np.pi / 180.0
        tmax = np.pi - tmin
        sigma = _xs.sigma_2to2_cm(
            lambda th: self.M_squared_avg(E_cm, th, m_e),
            s, (m_e, m_e), (m_e, m_e),
            theta_min=tmin, theta_max=tmax,
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    def __init__(self):
        self.diagrams = SimpleNamespace(
            t_channel=_diag_mod.Diagram(
                process_name="moller", channel="t_channel",
                expression=_draw.EXPRESSIONS["moller"],
                feynman_rules={"propagator": "photon", "vertex": "QED",
                               "fermi_sign": "negative (M_t - M_u)"},
                _render_fn=lambda ax: _draw.render_moller(ax, "t"),
                _tikz_template=_diag_mod.TIKZ_MOLLER_T,
            ),
            u_channel=_diag_mod.Diagram(
                process_name="moller", channel="u_channel",
                expression=_draw.EXPRESSIONS["moller"],
                feynman_rules={"propagator": "photon", "vertex": "QED",
                               "fermi_sign": "negative (M_t - M_u)"},
                _render_fn=lambda ax: _draw.render_moller(ax, "u"),
                _tikz_template=_diag_mod.TIKZ_MOLLER_U,
            ),
        )


class _BhabhaProcess:
    name = "bhabha"

    def M_squared_avg(self, E_cm: float, theta: float,
                      m_e: float = 0.510998928e-3) -> float:
        M, *_ = _proc.bhabha.M_squared_avg(E_cm, theta, m_e, _GAMMAS)
        return M

    def dsigma_dOmega(self, E_cm: float, theta: float,
                      m_e: float = 0.510998928e-3,
                      units: str = "pb") -> float:
        M_sq = self.M_squared_avg(E_cm, theta, m_e)
        s = E_cm ** 2
        ds = _xs.dsigma_dOmega_cm(M_sq, s, (m_e, m_e), (m_e, m_e))
        return ds * GEV_INV_SQ_TO_PB if units == "pb" else ds

    def sigma_total(self, E_cm: float, cutoff_deg: float = 30.0,
                    m_e: float = 0.510998928e-3,
                    units: str = "pb") -> float:
        s = E_cm ** 2
        tmin = cutoff_deg * np.pi / 180.0
        tmax = np.pi - tmin
        sigma = _xs.sigma_2to2_cm(
            lambda th: self.M_squared_avg(E_cm, th, m_e),
            s, (m_e, m_e), (m_e, m_e),
            theta_min=tmin, theta_max=tmax,
        )
        return sigma * GEV_INV_SQ_TO_PB if units == "pb" else sigma

    def __init__(self):
        self.diagrams = SimpleNamespace(
            s_channel=_diag_mod.Diagram(
                process_name="bhabha", channel="s_channel",
                expression=r"$i\mathcal{M}_s = (i e^2 / s) [\bar{v}_2\gamma^\mu u_1][\bar{u}_4\gamma_\mu v_3]$",
                feynman_rules={"propagator": "photon", "vertex": "QED"},
                _render_fn=lambda ax: _draw.render_eemumu(ax),
                _tikz_template=_diag_mod.TIKZ_BHABHA_S,
            ),
            t_channel=_diag_mod.Diagram(
                process_name="bhabha", channel="t_channel",
                expression=r"$i\mathcal{M}_t = -(i e^2 / t) [\bar{u}_3\gamma^\mu u_1][\bar{v}_2\gamma_\mu v_4]$",
                feynman_rules={"propagator": "photon", "vertex": "QED"},
                _render_fn=lambda ax: _draw.render_moller(ax, "t"),
                _tikz_template=_diag_mod.TIKZ_BHABHA_T,
            ),
        )


# ---------------------------------------------------------------------------
# Muon decay
# ---------------------------------------------------------------------------

class _MuonDecayProcess:
    name = "muon_decay"

    def Gamma(self, m_mu: float = 0.10566) -> float:
        """Total muon decay rate in GeV (tree-level Fermi)."""
        return _proc.muon_decay.gamma_via_michel(m_mu)

    def lifetime(self, m_mu: float = 0.10566, units: str = "us") -> float:
        """Muon lifetime in microseconds (tree-level)."""
        Gamma_GeV = self.Gamma(m_mu)
        tau_s = _proc.muon_decay.lifetime_seconds(Gamma_GeV)
        return tau_s * 1e6 if units == "us" else tau_s

    def michel_spectrum(self, x: float, m_mu: float = 0.10566) -> float:
        """dGamma/dx_e where x_e = 2 E_e / m_mu."""
        return _proc.muon_decay.dGamma_dx_michel(x, m_mu)

    def __init__(self):
        self.diagrams = SimpleNamespace(
            tree=_diag_mod.Diagram(
                process_name="muon_decay", channel="tree",
                expression=_draw.EXPRESSIONS["muon_decay"],
                feynman_rules={"propagator": "W (heavy, contracted)",
                               "vertex": "V-A"},
                _render_fn=lambda ax: _draw.render_muon_decay(ax),
                _tikz_template=_diag_mod.TIKZ_MUON_DECAY,
            ),
        )


# ---------------------------------------------------------------------------
# Singletons (so users can write `qed.compton.sigma_total(...)`)
# ---------------------------------------------------------------------------

compton = _ComptonProcess()
eemumu = _EEMuMuProcess()
moller = _MollerProcess()
bhabha = _BhabhaProcess()
muon_decay = _MuonDecayProcess()
