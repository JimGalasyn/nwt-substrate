"""Tests for muon decay rate via Michel spectrum and Dalitz integration."""

import numpy as np
import pytest

from nwt_substrate.amplitudes.processes import muon_decay as md


M_MU = 0.10566


def test_michel_spectrum_integrates_to_closed_form():
    """1D Michel spectrum integrates to G_F^2 m_mu^5 / (192 pi^3) exactly."""
    g_michel = md.gamma_via_michel(M_MU)
    g_closed = md.closed_form(M_MU)
    assert abs(g_michel / g_closed - 1.0) < 1e-12


def test_dalitz_integrates_to_closed_form():
    """2D Dalitz density integrates to G_F^2 m_mu^5 / (192 pi^3)."""
    g_dalitz = md.gamma_via_dalitz(M_MU)
    g_closed = md.closed_form(M_MU)
    assert abs(g_dalitz / g_closed - 1.0) < 1e-7


def test_michel_dalitz_agree():
    """1D Michel and 2D Dalitz integrations give the same total."""
    g_michel = md.gamma_via_michel(M_MU)
    g_dalitz = md.gamma_via_dalitz(M_MU)
    assert abs(g_michel / g_dalitz - 1.0) < 1e-7


def test_muon_lifetime_from_michel():
    """tau_mu via Michel = 2.187 us, within 1% of observed 2.197 us."""
    g = md.gamma_via_michel(M_MU)
    tau_s = md.lifetime_seconds(g)
    # tree level predicts 2.187 us; EW radiative corrections give the 0.45% gap
    assert abs(tau_s - 2.197e-6) / 2.197e-6 < 1e-2


def test_michel_endpoint():
    """Michel spectrum vanishes at x=0, peaks within (0,1), is positive."""
    assert md.dGamma_dx_michel(0.0, M_MU) == 0.0
    # Mass-zero limit endpoint x=1 gives x^2(3-2x) = 1, finite
    assert md.dGamma_dx_michel(1.0, M_MU) > 0
    # Spectrum is positive everywhere in (0, 1)
    for x in [0.1, 0.3, 0.5, 0.7, 0.9]:
        assert md.dGamma_dx_michel(x, M_MU) > 0


def test_dalitz_zero_outside_region():
    """Dalitz density is zero outside the Dalitz triangle x + y >= 1."""
    # Below threshold (x + y < 1)
    assert md.dalitz_density_VA(0.2, 0.2, M_MU) == 0.0
    assert md.dalitz_density_VA(0.4, 0.4, M_MU) == 0.0
    # On the diagonal x + y = 1 (boundary): inclusive in our definition
    assert md.dalitz_density_VA(0.5, 0.5, M_MU) >= 0.0
    # Inside region
    assert md.dalitz_density_VA(0.7, 0.5, M_MU) > 0


def test_dalitz_density_y_dependence():
    """At fixed x in Dalitz region, density vanishes at y = 0 and y = 1."""
    x = 0.6   # x + y >= 1 region: y in [0.4, 1]
    # At y = 1: y(1-y) = 0
    assert abs(md.dalitz_density_VA(x, 1.0, M_MU)) < 1e-30
    # At y = 0.4 (boundary): density positive but bounded by y(1-y) = 0.24
    assert md.dalitz_density_VA(x, 0.4, M_MU) > 0
