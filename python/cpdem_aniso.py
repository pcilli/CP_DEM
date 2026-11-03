"""Python mirror of the anisotropic (aligned-spheroid) cross-property DEM extension.

Implements the closure-based model of Cilli, "Anisotropic cross-property differential
effective medium modelling": transversely isotropic electrical DEM with isotropic-host
depolarisation factors (eqs 7-12, partial alignment eqs 54-56), TI elastic DEM with the
isotropic-host spheroidal Eshelby tensor at the Voigt-Reuss-Hill Poisson ratio (eqs
20-39, 57), and the porosity-free cross-property system (eqs 40-52, 63). Mirrors the
MATLAB/Octave functions cpdem_EshelbySpheroid, cpdem_VRH, cpdem_AnisoElecDEM,
cpdem_AnisoElasDEM, cpdem_AnisoForward, cpdem_AnisoInverse and cpdem_TIVelocities.
"""
from __future__ import annotations

import math

import numpy as np

import cpdem
from cpdem import _integrate, aligned_increment, shape

# Series coefficients (a nu + b)/den for each Eshelby component times (1 - nu), e = alpha - 1
_SER = {
    "S1111": [(-5, 7, 15), (-14, 16, 105), (3, -4, 35), (-176, 268, 3465), (260, -445, 9009),
              (-80, 152, 5005), (952, -1988, 109395), (-9728, 22144, 2078505), (4032, -9936, 1616615)],
    "S3333": [(-5, 7, 15), (28, -32, 105), (-18, 16, 105), (352, -224, 3465), (-520, 200, 9009),
              (480, -64, 15015), (-1904, -224, 109395), (19456, 7168, 2078505), (-8064, -4992, 1616615)],
    "S1122": [(5, -1, 15), (14, -4, 105), (-9, 2, 105), (176, -28, 3465), (-260, 25, 9009),
              (240, -8, 15015), (-952, -28, 109395), (9728, 896, 2078505), (-4032, -624, 1616615)],
    "S1133": [(5, -1, 15), (14, 2, 105), (-9, 1, 105), (176, -64, 3465), (-260, 160, 9009),
              (240, -208, 15015), (-952, 1064, 109395), (9728, -13312, 2078505), (-4032, 6528, 1616615)],
    "S3311": [(5, -1, 15), (-28, 2, 105), (18, 1, 105), (-352, -64, 3465), (520, 160, 9009),
              (-480, -208, 15015), (1904, 1064, 109395), (-19456, -13312, 2078505), (8064, 6528, 1616615)],
    "S1212": [(-5, 4, 15), (-14, 10, 105), (9, -7, 105), (-176, 148, 3465), (260, -235, 9009),
              (-240, 232, 15015), (952, -980, 109395), (-9728, 10624, 2078505), (4032, -4656, 1616615)],
    "S1313": [(-5, 4, 15), (7, -5, 105), (-9, 11, 210), (88, -152, 3465), (-130, 290, 9009),
              (120, -328, 15015), (-476, 1540, 109395), (4864, -18176, 2078505), (-2016, 8544, 1616615)],
}


def eshelby_spheroid(nu: float, alpha: float) -> np.ndarray:
    """6x6 Eshelby tensor (engineering-shear Voigt) of a spheroid (axis x3) in an isotropic host."""
    e = alpha - 1.0
    if abs(e) < 0.02:
        def ser(name):
            c = [(a*nu + b)/d for a, b, d in _SER[name]]
            v = c[-1]
            for ci in reversed(c[:-1]):
                v = v*e + ci
            return v/(1 - nu)
        S1111, S3333, S1122, S1133, S3311, S1212, S1313 = (ser(k) for k in
            ("S1111", "S3333", "S1122", "S1133", "S3311", "S1212", "S1313"))
    else:
        q = alpha
        g = shape(alpha)[0]
        d = q*q - 1
        S1111 = 3/(8*(1 - nu))*q*q/d + 1/(4*(1 - nu))*(1 - 2*nu - 9/(4*d))*g
        S3333 = 1/(2*(1 - nu))*(1 - 2*nu + (3*q*q - 1)/d - (1 - 2*nu + 3*q*q/d)*g)
        S1122 = 1/(4*(1 - nu))*(q*q/(2*d) - (1 - 2*nu + 3/(4*d))*g)
        S1133 = -1/(2*(1 - nu))*q*q/d + 1/(4*(1 - nu))*(3*q*q/d - (1 - 2*nu))*g
        S3311 = -1/(2*(1 - nu))*(1 - 2*nu + 1/d) + 1/(2*(1 - nu))*(1 - 2*nu + 3/(2*d))*g
        S1212 = 1/(4*(1 - nu))*(q*q/(2*d) + (1 - 2*nu - 3/(4*d))*g)
        S1313 = 1/(4*(1 - nu))*(1 - 2*nu - (q*q + 1)/d - 0.5*(1 - 2*nu - 3*(q*q + 1)/d)*g)
    S = np.zeros((6, 6))
    S[0, 0] = S[1, 1] = S1111; S[2, 2] = S3333
    S[0, 1] = S[1, 0] = S1122
    S[0, 2] = S[1, 2] = S1133
    S[2, 0] = S[2, 1] = S3311
    S[3, 3] = S[4, 4] = 2*S1313
    S[5, 5] = 2*S1212
    return S


def ti_stiffness(c5) -> np.ndarray:
    C11, C13, C33, C44, C66 = c5
    C12 = C11 - 2*C66
    C = np.zeros((6, 6))
    C[0, 0] = C[1, 1] = C11; C[0, 1] = C[1, 0] = C12
    C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = C13; C[2, 2] = C33
    C[3, 3] = C[4, 4] = C44; C[5, 5] = C66
    return C


def vrh(C: np.ndarray):
    """Voigt-Reuss-Hill K, mu, nu of a TI stiffness (eqs 27-32)."""
    s = np.linalg.inv(C)
    C11, C12, C13, C33, C44, C66 = C[0, 0], C[0, 1], C[0, 2], C[2, 2], C[3, 3], C[5, 5]
    KV = (2*C11 + C33 + 2*C12 + 4*C13)/9
    muV = (2*C11 + C33 - C12 - 2*C13 + 6*C44 + 3*C66)/15
    KR = 1/(2*s[0, 0] + 2*s[0, 1] + 4*s[0, 2] + s[2, 2])
    muR = 15/(8*s[0, 0] + 4*s[2, 2] - 4*s[0, 1] - 8*s[0, 2] + 6*s[3, 3] + 3*s[5, 5])
    K, mu = (KV + KR)/2, (muV + muR)/2
    return K, mu, (3*K - 2*mu)/(2*(3*K + mu))


def aniso_rates(sp, sq, c5, sig2, k2, alpha, S, L, closure="isotropic-host"):
    """(1 - phi) d/dphi of [sigPar, sigPerp, C11, C13, C33, C44, C66]."""
    if closure == "isotropic-host":
        m_par = sp/(sp + L*(sig2 - sp))
        m_perp = sq/(sq + (1 - L)/2*(sig2 - sq))
        m_rand = (2*m_perp + m_par)/3
        dsp = (sig2 - sp)*(m_rand + S*(m_par - m_rand))
        dsq = (sig2 - sq)*(m_rand + S*(m_perp - m_rand))
    elif closure == "exact":
        dsp, dsq = aligned_increment(sp, sq, sig2, alpha)
    else:
        raise ValueError("closure must be 'isotropic-host' or 'exact'")
    C = ti_stiffness(c5)
    _, _, nu = vrh(C)
    if S == 1:
        Scl = eshelby_spheroid(nu, alpha)
    elif S == 0:
        Scl = eshelby_spheroid(nu, 1.0)
    else:
        Scl = (1 - S)*eshelby_spheroid(nu, 1.0) + S*eshelby_spheroid(nu, alpha)
    C2 = np.zeros((6, 6)); C2[:3, :3] = k2
    dC = C2 - C
    A = np.linalg.inv(np.eye(6) + Scl @ np.linalg.solve(C, dC))
    M = dC @ A
    return np.array([dsp, dsq, M[0, 0], M[0, 2], M[2, 2], M[3, 3], M[5, 5]])


def _c0(k1, mu1):
    lam = k1 - 2*mu1/3
    return np.array([lam + 2*mu1, lam, lam + 2*mu1, mu1, mu1])


def aniso_elec_dem(phi, sig1, sig2, alpha, S=1.0, closure="isotropic-host", rtol=1e-10, atol=1e-13):
    """TI electrical DEM (eqs 11-12 with 54-56): returns (sigPar, sigPerp)."""
    phi = np.asarray(phi, float)
    L = shape(alpha)[2]
    def deriv(p, y):
        r = aniso_rates(y[0], y[1], _c0(1.0, 0.5), sig2, 0.1, alpha, S, L, closure)
        return [r[0]/(1 - p), r[1]/(1 - p)]
    Y = _integrate(deriv, phi.ravel(), [sig1, sig1], rtol, atol)
    return Y[:, 0].reshape(phi.shape), Y[:, 1].reshape(phi.shape)


def aniso_elas_dem(phi, k1, mu1, k2, alpha, S=1.0, rtol=1e-10, atol=1e-13):
    """TI elastic DEM (eqs 20-39, 57): returns dict C11, C12, C13, C33, C44, C66."""
    phi = np.asarray(phi, float)
    L = shape(alpha)[2]
    def deriv(p, y):
        r = aniso_rates(1.0, 1.0, y*k1, 1.0, k2, alpha, S, L)
        return r[2:]/(k1*(1 - p))
    Y = _integrate(deriv, phi.ravel(), _c0(k1, mu1)/k1, rtol, atol)*k1
    out = {k: Y[:, i].reshape(phi.shape) for i, k in enumerate(("C11", "C13", "C33", "C44", "C66"))}
    out["C12"] = out["C11"] - 2*out["C66"]
    return out


def aniso_forward(k1, mu1, sig1, k2, sig2, sig0, alpha, S=1.0, which="par", closure="isotropic-host",
                  rtol=1e-10, atol=1e-13):
    """Porosity-free cross-property system (eqs 40-52): TI stiffness and the other
    conductivity component from one measured component. Returns a dict."""
    use_par = which.lower() in ("par", "parallel", "axial", "axis")
    sig0 = np.asarray(sig0, float)
    s = sig0.ravel()
    L = shape(alpha)[2]
    lo, hi = min(sig1, sig2), max(sig1, sig2)
    ok = (s >= lo) & (s <= hi) & (s != sig2)
    c0 = _c0(k1, mu1)/k1

    def deriv(u, y):
        sw, so = sig1*math.exp(u), sig1*math.exp(y[0])
        sp, sq = (sw, so) if use_par else (so, sw)
        r = aniso_rates(sp, sq, y[1:6]*k1, sig2, k2, alpha, S, L, closure)
        rate = r[0]/sp if use_par else r[1]/sq
        other = r[1]/sq if use_par else r[0]/sp
        return np.concatenate(([other], r[2:]/k1, [1 - y[6]]))/rate

    Y = np.full((s.size, 7), np.nan)
    if ok.any():
        Y[ok] = _integrate(deriv, np.log(s[ok]/sig1), np.concatenate(([0.0], c0, [0.0])), rtol, atol)
    other = sig1*np.exp(Y[:, 0])
    out = {k: (Y[:, i + 1]*k1).reshape(sig0.shape) for i, k in enumerate(("C11", "C13", "C33", "C44", "C66"))}
    out["C12"] = out["C11"] - 2*out["C66"]
    sp = np.where(ok, s, np.nan) if use_par else other
    sq = other if use_par else np.where(ok, s, np.nan)
    out["sigPar"], out["sigPerp"], out["phi"] = sp.reshape(sig0.shape), sq.reshape(sig0.shape), Y[:, 6].reshape(sig0.shape)
    return out


def aniso_inverse(c_meas, which, k1, mu1, sig1, k2, sig2, alpha, S=1.0, closure="isotropic-host",
                  rtol=1e-10, atol=1e-13):
    """Inverse cross-property DEM (section 2.7): conductivities from one TI stiffness."""
    names = ["C11", "C13", "C33", "C44", "C66"]
    idx = names.index(which)
    if idx == 1:
        raise ValueError("use C11, C33, C44 or C66")
    others = [i for i in range(5) if i != idx]
    L = shape(alpha)[2]
    c0 = _c0(k1, mu1)
    c_start, c_end = c0[idx], (k2 if idx <= 2 else 0.0)
    c_meas = np.asarray(c_meas, float)
    c = c_meas.ravel()
    ok = (c <= c_start) & (c > c_end)

    def deriv(x, y):
        c5 = np.zeros(5); c5[idx] = c_start + x*k1; c5[others] = y[2:6]*k1
        sp, sq = sig1*math.exp(y[0]), sig1*math.exp(y[1])
        r = aniso_rates(sp, sq, c5, sig2, k2, alpha, S, L, closure)
        rate = r[2 + idx]/k1
        return np.concatenate(([r[0]/sp, r[1]/sq], r[2 + np.array(others)]/k1, [1 - y[6]]))/rate

    Y = np.full((c.size, 7), np.nan)
    if ok.any():
        Y[ok] = _integrate(deriv, (c[ok] - c_start)/k1, np.concatenate(([0.0, 0.0], c0[others]/k1, [0.0])), rtol, atol)
    c5 = np.full((c.size, 5), np.nan)
    c5[:, others] = Y[:, 2:6]*k1
    c5[:, idx] = np.where(ok, c, np.nan)
    out = {k: c5[:, i].reshape(c_meas.shape) for i, k in enumerate(names)}
    out["C12"] = out["C11"] - 2*out["C66"]
    out["sigPar"] = (sig1*np.exp(Y[:, 0])).reshape(c_meas.shape)
    out["sigPerp"] = (sig1*np.exp(Y[:, 1])).reshape(c_meas.shape)
    out["phi"] = Y[:, 6].reshape(c_meas.shape)
    return out


def ti_velocities(C11, C13, C33, C44, C66, rho):
    """Phase velocities along/across the axis and Thomsen parameters."""
    C11, C13, C33, C44, C66 = (np.asarray(x, float) for x in (C11, C13, C33, C44, C66))
    return {"vp0": np.sqrt(C33/rho), "vp90": np.sqrt(C11/rho), "vsv": np.sqrt(C44/rho),
            "vsh90": np.sqrt(C66/rho), "epsilon": (C11 - C33)/(2*C33), "gamma": (C66 - C44)/(2*C44),
            "delta": ((C13 + C44)**2 - (C33 - C44)**2)/(2*C33*(C33 - C44))}
