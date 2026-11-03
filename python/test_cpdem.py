"""Tests of the Python port, including a cross-check against the Octave toolbox."""
import json
import math
import os

import numpy as np
import pytest

import cpdem

P = dict(k1=37e9, mu1=44e9, sig1=1e-4, k2=2.25e9, mu2=0.0, sig2=5.0)
SIG0 = np.array([0.02, 0.05, 0.1, 0.2, 0.4])
HERE = os.path.dirname(os.path.abspath(__file__))


def test_shape_sphere_limits_and_identity():
    t, f, L = cpdem.shape(1.0)
    assert abs(t - 2/3) < 1e-15 and abs(f + 0.4) < 1e-15 and abs(L - 1/3) < 1e-15
    for a in [1e-3, 0.1, 0.5, 0.9, 1.1, 2, 16.4, 100, 1e3]:
        t, _, L = cpdem.shape(a)
        assert abs(L - (1 - t)) < 1e-10
    for e in (-0.1, 0.1):
        s1 = cpdem.shape(1 + e*(1 - 1e-9)); s2 = cpdem.shape(1 + e*(1 + 1e-9))
        assert np.allclose(s1, s2, rtol=1e-9)


def test_PQ_mbar_closed_forms():
    k0, mu0 = 30e9, 35e9
    p, q = cpdem.PQ(k0, mu0, P["k2"], P["mu2"], 1.0)
    zeta = mu0/6*(9*k0 + 8*mu0)/(k0 + 2*mu0)
    assert abs(p - (k0 + 4/3*mu0)/(P["k2"] + 4/3*mu0)) < 1e-12
    assert abs(q - (mu0 + zeta)/(P["mu2"] + zeta)) < 1e-12
    assert abs(cpdem.mbar(0.1, 5, 1.0) - 3*0.1/(2*0.1 + 5)) < 1e-13
    assert abs(cpdem.mbar(0.3, 0, 1.0) - 1.5) < 1e-13
    # v1.0 toolbox values
    p, q = cpdem.PQ(k0, mu0, P["k2"], P["mu2"], 16.4)
    assert np.allclose([p, q, cpdem.mbar(0.1, 5, 16.4)], [1.7378004867, 2.50792849099, 0.255250509316], rtol=1e-10)


def test_forward_composition_and_monotonic():
    for epar in (16.4, 1.0, 0.1):
        k0, mu0, phi0 = cpdem.forward(**P, sig0=SIG0, epar=epar)
        phi = cpdem.implied_porosity(SIG0, P["sig1"], P["sig2"], epar)
        assert np.allclose(phi, phi0, rtol=1e-7, atol=1e-10)
        k, mu = cpdem.elas_dem(phi, P["k1"], P["mu1"], P["k2"], P["mu2"], epar)
        assert np.allclose(k, k0, rtol=1e-7) and np.allclose(mu, mu0, rtol=1e-7)
        assert np.allclose(cpdem.elec_dem(phi, P["sig1"], P["sig2"], epar), SIG0, rtol=1e-7)
        assert np.all(np.diff(k0) < 0) and np.all(np.diff(mu0) < 0) and np.all(np.diff(phi0) > 0)


def test_mbar_and_R_forms_agree():
    a = cpdem.forward(**P, sig0=SIG0, epar=16.4, form="mbar")
    b = cpdem.forward(**P, sig0=SIG0, epar=16.4, form="R")
    for x, y in zip(a, b):
        assert np.allclose(x, y, rtol=1e-12, atol=1e-14)


def test_bruggeman_closed_form():
    phi = np.array([0.05, 0.1, 0.2, 0.3])
    s = cpdem.elec_dem(phi, P["sig1"], P["sig2"], 1.0)
    lhs = (s - P["sig2"])/(P["sig1"] - P["sig2"])*(P["sig1"]/s)**(1/3)
    assert np.allclose(lhs, 1 - phi, atol=1e-9)


def test_inverse_round_trip():
    k0, mu0, phi0 = cpdem.forward(**P, sig0=SIG0, epar=16.4)
    s, phi, mu = cpdem.inverse(k0, "K", **P, epar=16.4)
    assert np.allclose(s, SIG0, rtol=1e-7) and np.allclose(mu, mu0, rtol=1e-7)
    s, phi, k = cpdem.inverse(mu0, "mu", **P, epar=16.4)
    assert np.allclose(s, SIG0, rtol=1e-7) and np.allclose(k, k0, rtol=1e-7)
    s, _, _ = cpdem.inverse(np.array([P["k1"]*1.01, k0[2], P["k2"]*0.5]), "K", **P, epar=16.4)
    assert np.isnan(s[0]) and np.isnan(s[2]) and abs(s[1] - SIG0[2]) < 1e-7*SIG0[2]


def test_epar_of_porosity():
    law = lambda phi: max(30*phi**0.3, 0.05)
    k0, mu0, phi0 = cpdem.forward(**P, sig0=SIG0, epar=law)
    phi = cpdem.implied_porosity(SIG0, P["sig1"], P["sig2"], law)
    assert np.allclose(phi, phi0, rtol=1e-7, atol=1e-10)
    k, mu = cpdem.elas_dem(phi, P["k1"], P["mu1"], P["k2"], P["mu2"], law)
    assert np.allclose(k, k0, rtol=1e-7) and np.allclose(mu, mu0, rtol=1e-7)
    s, _, _ = cpdem.inverse(k0, "K", **P, epar=law)
    assert np.allclose(s, SIG0, rtol=1e-7)
    with pytest.raises(ValueError):
        cpdem.forward(**P, sig0=0.1, epar=lambda phi: 30*phi**0.3)


def test_calibration_recovers_epars():
    k0, _, _ = cpdem.forward(**P, sig0=SIG0, epar=16.4)
    _, mu0, _ = cpdem.forward(**P, sig0=SIG0, epar=12.8)
    cal = cpdem.calibrate_epar(k0, mu0, SIG0, **P)
    assert abs(cal["eparK"] - 16.4) < 2e-3 and abs(cal["eparMu"] - 12.8) < 2e-3
    assert cal["r2K"] > 0.999999 and cal["r2Mu"] > 0.999999


def test_bounds_contain_model_and_vpvs():
    sig = np.linspace(0.005, 2, 25)
    k_up, k_low, mu_up, mu_low, phi_low, phi_up = cpdem.elec_elas_hs(sig, P["sig1"], P["sig2"], P["k1"], P["k2"], P["mu1"])
    for epar in (0.05, 1.0, 16.4):
        k0, mu0, phi0 = cpdem.forward(**P, sig0=sig, epar=epar)
        assert np.all(k0 <= k_up*(1 + 1e-9)) and np.all(k0 >= k_low*(1 - 1e-9))
        assert np.all(mu0 <= mu_up*(1 + 1e-9)) and np.all(phi0 >= phi_low*(1 - 1e-9)) and np.all(phi0 <= phi_up*(1 + 1e-9))
    out = cpdem.model_vpvs(**P, sig0=SIG0, kepar=16.4, muepar=12.8, density=(2650.0, 1050.0))
    assert np.all((out["vpvs"] > 1.4) & (out["vpvs"] < 2.5))
    assert np.allclose(out["vp"], np.sqrt(out["M"]/out["rho"]))


def test_thermal_electrical():
    lam, phi = cpdem.forward_conduct(7.7, P["sig1"], 0.6, P["sig2"], SIG0, 16.4)
    assert np.allclose(lam, cpdem.elec_dem(phi, 7.7, 0.6, 16.4), rtol=1e-7)
    sig, _ = cpdem.forward_conduct(P["sig1"], 7.7, P["sig2"], 0.6, lam, 16.4)
    assert np.allclose(sig, SIG0, rtol=1e-6)


def test_against_octave_reference():
    """Cross-check with values written by Tests/cpdem_WriteReference.m (GNU Octave)."""
    path = os.path.join(HERE, "reference_values.json")
    if not os.path.exists(path):
        pytest.skip("reference_values.json not present")
    with open(path) as fh:
        ref = json.load(fh)
    prm = ref["params"]
    for case in ref["cases"]:
        epar = case["epar"]
        k0, mu0, phi0 = cpdem.forward(prm["k1"], prm["mu1"], prm["sig1"], prm["k2"], prm["mu2"],
                                      prm["sig2"], np.array(prm["sig0"]), epar)
        assert np.allclose(k0, case["k0"], rtol=1e-7), f"K, EPAR {epar}"
        assert np.allclose(mu0, case["mu0"], rtol=1e-7), f"mu, EPAR {epar}"
        assert np.allclose(phi0, case["phi0"], rtol=1e-6, atol=1e-10), f"phi, EPAR {epar}"
        s, _, _ = cpdem.inverse(np.array(case["k0"]), "K", prm["k1"], prm["mu1"], prm["sig1"],
                                prm["k2"], prm["mu2"], prm["sig2"], epar)
        assert np.allclose(s, prm["sig0"], rtol=1e-6), f"inverse, EPAR {epar}"
        lam, _ = cpdem.forward_conduct(ref["thermal"]["lam1"], prm["sig1"], ref["thermal"]["lam2"],
                                       prm["sig2"], np.array(prm["sig0"]), epar)
        assert np.allclose(lam, case["lam"], rtol=1e-7), f"thermal, EPAR {epar}"
