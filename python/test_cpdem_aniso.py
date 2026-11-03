"""Tests of the anisotropic cross-property DEM mirror, with an Octave cross-check."""
import json
import math
import os

import numpy as np
import pytest

import cpdem
import cpdem_aniso as an

HERE = os.path.dirname(os.path.abspath(__file__))


def test_eshelby_sphere_series_and_berryman():
    nu = 0.25
    S = an.eshelby_spheroid(nu, 1.0)
    a = (7 - 5*nu)/(15*(1 - nu)); b = (5*nu - 1)/(15*(1 - nu)); c = 2*(4 - 5*nu)/(15*(1 - nu))
    assert np.allclose([S[0, 0], S[0, 1], S[0, 2], S[2, 0], S[2, 2], S[3, 3], S[5, 5]], [a, b, b, b, a, c, c], atol=1e-14)
    for e in (-0.02, 0.02):
        assert np.allclose(an.eshelby_spheroid(nu, 1 + e*(1 - 1e-9)), an.eshelby_spheroid(nu, 1 + e*(1 + 1e-9)), rtol=1e-9, atol=1e-12)
    k1, mu1 = 37e9, 44e9
    nu1 = (3*k1 - 2*mu1)/(2*(3*k1 + mu1))
    C = an.ti_stiffness([k1 + 4/3*mu1, k1 - 2/3*mu1, k1 + 4/3*mu1, mu1, mu1])
    for k2, mu2 in ((2.25e9, 0.0), (10e9, 5e9), (80e9, 60e9)):
        C2 = an.ti_stiffness([k2 + 4/3*mu2, k2 - 2/3*mu2, k2 + 4/3*mu2, mu2, mu2])
        for alpha in (0.05, 0.2, 0.9, 1.0, 1.2, 5, 30):
            A = np.linalg.inv(np.eye(6) + an.eshelby_spheroid(nu1, alpha) @ np.linalg.solve(C, C2 - C))
            Aiijj = A[:3, :3].sum(); Aijij = np.trace(A[:3, :3]) + A[3, 3] + A[4, 4] + A[5, 5]
            P, Q = cpdem.PQ(k1, mu1, k2, mu2, alpha)
            assert abs(Aiijj/3 - P) < 1e-10*P and abs((Aijij - Aiijj/3)/5 - Q) < 1e-10*Q


def test_limits_and_section3():
    sig1, K1, MU1, K2 = 1.0, 40e9, 30e9, 2.2e9
    phi = np.array([0.1, 0.2, 0.3])
    sp, sq = an.aniso_elec_dem(phi, sig1, 10.0, 1.0)
    assert np.allclose(sp, cpdem.elec_dem(phi, sig1, 10.0, 1.0), rtol=1e-9) and np.allclose(sp, sq)
    el = an.aniso_elas_dem(phi, K1, MU1, K2, 1.0)
    k, m = cpdem.elas_dem(phi, K1, MU1, K2, 0.0, 1.0)
    assert np.allclose(el["C11"], k + 4/3*m, rtol=1e-9) and np.allclose(el["C44"], m, rtol=1e-9)
    sp, sq = an.aniso_elec_dem(phi, sig1, 10.0, 0.1, S=0.0)
    assert np.allclose(sp, cpdem.elec_dem(phi, sig1, 10.0, 0.1), rtol=1e-9)
    L = cpdem.shape(0.1)[2]
    sp, sq = an.aniso_elec_dem(phi, sig1, 1e-300, 0.1)
    assert np.allclose(sp, (1 - phi)**(1/(1 - L)), rtol=1e-6) and np.allclose(sq, (1 - phi)**(2/(1 + L)), rtol=1e-6)
    paper = {10: (10.0, 3.54, 1.74, 32.39, 54.15, 14.90, 10.06), 1: (10.0, 2.03, 2.03, 41.11, 41.11, 14.74, 14.74),
             0.1: (10.0, 1.43, 3.17, 52.54, 9.88, 2.70, 20.05)}
    for alpha, (sig2, psp, psq, c11, c33, c44, c66) in paper.items():
        sp, sq = an.aniso_elec_dem(np.array([0.3]), sig1, sig2, alpha)
        el = an.aniso_elas_dem(np.array([0.3]), K1, MU1, K2, alpha)
        assert abs(sp[0] - psp) < 0.006 and abs(sq[0] - psq) < 0.006
        assert np.allclose([el["C11"][0], el["C33"][0], el["C44"][0], el["C66"][0]], np.array([c11, c33, c44, c66])*1e9, atol=6e6)


def test_cross_property_and_inverse():
    sig1, K1, MU1, K2 = 1.0, 40e9, 30e9, 2.2e9
    phi = np.array([0.1, 0.2, 0.3])
    for S in (1.0, 0.6):
        for alpha in (10.0, 0.1):
            sp, sq = an.aniso_elec_dem(phi, sig1, 10.0, alpha, S)
            el = an.aniso_elas_dem(phi, K1, MU1, K2, alpha, S)
            f = an.aniso_forward(K1, MU1, sig1, K2, 10.0, sp, alpha, S, "par")
            assert np.allclose(f["sigPerp"], sq, rtol=1e-7) and np.allclose(f["phi"], phi, atol=1e-9)
            for k in ("C11", "C13", "C33", "C44", "C66"):
                assert np.allclose(f[k], el[k], rtol=1e-7)
            g = an.aniso_forward(K1, MU1, sig1, K2, 10.0, sq, alpha, S, "perp")
            assert np.allclose(g["sigPar"], sp, rtol=1e-7) and np.allclose(g["C33"], el["C33"], rtol=1e-7)
    sp, sq = an.aniso_elec_dem(phi, sig1, 0.05, 0.1, 1.0, closure="exact")
    ap, aq = cpdem.aligned_cond_dem(phi, sig1, 0.05, 0.1)
    assert np.allclose(sp, ap, rtol=1e-9) and np.allclose(sq, aq, rtol=1e-9)
    sp, sq = an.aniso_elec_dem(phi, sig1, 10.0, 0.1)
    el = an.aniso_elas_dem(phi, K1, MU1, K2, 0.1)
    for w in ("C33", "C11", "C44", "C66"):
        inv = an.aniso_inverse(el[w], w, K1, MU1, sig1, K2, 10.0, 0.1)
        assert np.allclose(inv["sigPar"], sp, rtol=1e-6) and np.allclose(inv["sigPerp"], sq, rtol=1e-6)
        assert np.allclose(inv["phi"], phi, atol=1e-8) and np.allclose(inv["C13"], el["C13"], rtol=1e-6)
    v = an.ti_velocities(80e9, 20e9, 80e9, 30e9, 30e9, 2500.0)
    assert abs(v["epsilon"]) < 1e-14 and abs(v["gamma"]) < 1e-14 and abs(v["delta"]) < 1e-14


def test_against_octave_reference():
    path = os.path.join(HERE, "reference_aniso.json")
    if not os.path.exists(path):
        pytest.skip("reference_aniso.json not present")
    with open(path) as fh:
        ref = json.load(fh)
    p = ref["params"]
    phi = np.array(p["phi"])
    for case in ref["cases"]:
        sp, sq = an.aniso_elec_dem(phi, p["sig1"], case["sig2"], case["alpha"], case["S"])
        assert np.allclose(sp, case["sigPar"], rtol=1e-7) and np.allclose(sq, case["sigPerp"], rtol=1e-7)
        el = an.aniso_elas_dem(phi, p["k1"], p["mu1"], p["k2"], case["alpha"], case["S"])
        for k in ("C11", "C13", "C33", "C44", "C66"):
            assert np.allclose(el[k], case[k], rtol=1e-7), (k, case)
        f = an.aniso_forward(p["k1"], p["mu1"], p["sig1"], p["k2"], case["sig2"], sp, case["alpha"], case["S"], "par")
        assert np.allclose(f["C33"], case["fwdC33"], rtol=1e-6) and np.allclose(f["sigPerp"], case["fwdSigPerp"], rtol=1e-6)
    assert np.allclose(an.eshelby_spheroid(0.3, 0.2).ravel(order="F"), ref["eshelby_nu03_a02"], rtol=1e-12, atol=1e-14)
    assert np.allclose(an.eshelby_spheroid(0.3, 1.01).ravel(order="F"), ref["eshelby_nu03_a101"], rtol=1e-10, atol=1e-13)
