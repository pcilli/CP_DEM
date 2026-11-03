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
    for a in (1e-11, 1e-9, 1e9, 1e11, 1e13):
        t, f, L = cpdem.shape(a)
        assert all(map(math.isfinite, (t, f, L))) and 0 <= L <= 1
    assert abs(cpdem.shape(1e9)[2] - (math.log(2e9) - 1)/1e18) < 1e-6*(math.log(2e9) - 1)/1e18
    assert cpdem.shape(0.0) == (0.0, 0.0, 1.0) and cpdem.shape(math.inf) == (1.0, -1.0, 0.0)


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


def test_anisotropic_and_aligned():
    phi = np.array([0.05, 0.1, 0.2, 0.3, 0.4])
    assert abs(cpdem.cond_increment(0.3, (5.0, 5.0), 16.4) - (5 - 0.3)*cpdem.mbar(0.3, 5.0, 16.4)) < 1e-14
    sA = cpdem.elec_dem(phi, 1.0, (100.0, 1e-3), 0.1)
    assert np.all((sA > cpdem.elec_dem(phi, 1.0, 1e-3, 0.1)) & (sA < cpdem.elec_dem(phi, 1.0, 100.0, 0.1)))
    spar, sperp = cpdem.aligned_cond_dem(phi, 1.0, 20.0, 1.0)
    assert np.allclose(spar, sperp, rtol=1e-10)
    assert np.allclose((spar - 20)/(1 - 20)*(1/spar)**(1/3), 1 - phi, atol=1e-8)
    spar, sperp = cpdem.aligned_cond_dem(phi, 1.0, 20.0, 1e-7)
    assert np.allclose(spar, 1/((1 - phi) + phi/20), rtol=1e-5)      # series
    assert np.allclose(sperp, (1 - phi) + 20*phi, rtol=1e-5)         # parallel
    spar, sperp = cpdem.aligned_cond_dem(phi, 1.0, 20.0, 1e7)
    assert np.allclose(spar, (1 - phi) + 20*phi, rtol=1e-5)
    assert np.allclose((sperp - 20)/(1 - 20)*(1/sperp)**0.5, 1 - phi, atol=1e-5)  # 2-D Bruggeman
    for epar in (0.1, 16.4):
        spar, sperp = cpdem.aligned_cond_dem(phi, P["sig1"], P["sig2"], epar)
        lpar, lperp = cpdem.aligned_cond_dem(phi, 7.7, 0.6, epar)
        o = cpdem.aligned_forward_conduct(7.7, P["sig1"], 0.6, P["sig2"], spar, epar, "par")
        assert np.allclose(o["sigPerp"], sperp, rtol=1e-6) and np.allclose(o["pPar"], lpar, rtol=1e-6)
        assert np.allclose(o["pPerp"], lperp, rtol=1e-6) and np.allclose(o["phi"], phi, atol=1e-8)


def test_thermal_electrical_paper():
    """Equations 5-9, mirror symmetry, HS bounds and the Section 3 calibration of
    Cilli & Chapman (arXiv:2102.12434) with the paper's Cu-graphite parameters."""
    c = cpdem.published_cu_graphite()
    k1, s1, k2, s2 = c["kappa1"], c["sig1"], c["kappa2"], c["sig2"]
    sig = s1*np.array([0.9, 0.5, 0.2, 0.05]); kap = k1*np.array([0.9, 0.6, 0.3, 0.1])
    rhs = lambda a: cpdem.cond_increment(kap, k2, a)/cpdem.cond_increment(sig, s2, a)
    eq7 = kap*(k2 - kap)*(2*sig + s2)/(sig*(s2 - sig)*(2*kap + k2))
    eq8 = s2*(k2 - kap)*(2*k2 + kap)/(k2*(s2 - sig)*(2*s2 + sig))
    eq9 = (s2 + sig)*(k2 - kap)*(k2 + 5*kap)/((k2 + kap)*(s2 - sig)*(s2 + 5*sig))
    assert np.allclose(rhs(1.0), eq7, rtol=1e-12) and np.allclose(rhs(0.0), eq8, rtol=1e-12)
    assert np.allclose(rhs(math.inf), eq9, rtol=1e-12)
    sig0 = s1*np.array([0.95, 0.7, 0.5, 0.3, 0.1, 0.03])
    k_up, k_low, phi_min, phi_max, k_on1, k_on2 = cpdem.cond_cond_hs(sig0, s1, s2, k1, k2)
    for a in c["alphaCurves"]:
        kap0, _ = cpdem.forward_conduct(k1, s1, k2, s2, sig0, a)
        back, _ = cpdem.forward_conduct(s1, k1, s2, k2, kap0, a)
        assert np.allclose(back, sig0, rtol=1e-7)
        assert np.all(kap0 <= k_up*(1 + 1e-9)) and np.all(kap0 >= k_low*(1 - 1e-7))
    kap_disk, _ = cpdem.forward_conduct(k1, s1, k2, s2, sig0, 0.0)
    assert np.allclose(kap_disk, k_on2, rtol=1e-8)
    phi = np.array([0.1, 0.3, 0.4])
    cal = cpdem.calibrate_conduct(cpdem.elec_dem(phi, k1, k2, c["alpha"]), cpdem.elec_dem(phi, s1, s2, c["alpha"]),
                                  k1, s1, (c["kappaFlakePlane"], c["kappaFlakeAxis"]), (c["sigFlakeAxis"], c["sigFlakePlane"]))
    assert abs(cal["epar"] - c["alpha"]) < 2e-2*c["alpha"] and abs(cal["p2"] - k2) < 2e-2*k2 and abs(cal["sig2"] - s2) < 2e-2*s2


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
        if "alignedPar" in case:
            o = cpdem.aligned_forward_conduct(ref["thermal"]["lam1"], prm["sig1"], ref["thermal"]["lam2"],
                                              prm["sig2"], np.array(prm["sig0"]), epar, "par")
            assert np.allclose(o["pPar"], case["alignedPar"], rtol=1e-6), f"aligned par, EPAR {epar}"
            assert np.allclose(o["pPerp"], case["alignedPerp"], rtol=1e-6), f"aligned perp, EPAR {epar}"
            assert np.allclose(o["sigPerp"], case["alignedSigPerp"], rtol=1e-6), f"aligned sigPerp, EPAR {epar}"
            lamA, _ = cpdem.forward_conduct(ref["thermal"]["lam1"], prm["sig1"], (ref["thermal"]["lam2"], 5*ref["thermal"]["lam2"]),
                                            prm["sig2"], np.array(prm["sig0"]), epar)
            assert np.allclose(lamA, case["anisoInclusion"], rtol=1e-7), f"anisotropic inclusion, EPAR {epar}"
