"""Reference Python port of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

Implements the cross-property differential effective medium (DEM) model of
Cilli and Chapman (2021), which links the elastic moduli and the electrical
conductivity of a two-phase composite through the equivalent pore aspect
ratio (EPAR) without reference to porosity, together with the single-
property DEMs it is built from, the direct inverse, the electrical-elastic
Hashin-Shtrikman bounds and the thermal-electrical cross-property model.

The MATLAB/Octave toolbox in ../Functions is the primary implementation; this
module mirrors it function for function so that the two can be checked
against each other (see test_cpdem.py). Only numpy and scipy are required.

Units: SI is recommended (Pa, S/m, W/(m K)).

References
----------
Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
properties of rocks using cross-property DEM. Geophysical Journal
International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046

Cilli, P.A., and Chapman, M. (2020), The power-law relation between inclusion
aspect ratio and porosity: Implications for electrical and elastic modeling.
JGR Solid Earth, 125, e2019JB019187, DOI:10.1029/2019JB019187

Cilli, P.A., and Chapman, M. (2021), Linking electrical and thermal
conductivity through cross-property inclusion modelling. arXiv:2102.12434

Berryman, J.G. (1980), J. Acoust. Soc. Am. 68, 1820-1831.
Berryman, J.G. (1992), J. Acoust. Soc. Am. 91, 551-571.
Osborn, J.A. (1945), Phys. Rev. 67, 351.
Carcione, J.M., Ursin, B., and Nordskag, J.I. (2007), Geophysics 72, E193-E204.
"""
from __future__ import annotations

import math
from typing import Callable, Union

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

__version__ = "2.0.0"

EPAR = Union[float, Callable[[float], float], np.ndarray]

# 13th-order Taylor coefficients of theta and f about the sphere (e = EPAR - 1)
_CT = [2/3, 4/15, -6/35, 32/315, -40/693, 32/1001, -112/6435, 1024/109395,
       -1152/230945, 2560/969969, -2816/2028117, 12288/16900975,
       -13312/35102025, 28672/145422675]
_CF = [-2/5, -12/35, 2/15, -32/1155, -40/3003, 32/1365, -784/36465,
       1024/62985, -18048/1616615, 7680/1062347, -76032/16900975,
       413696/152108775, -545792/339319575, 1404928/1502700975]


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
def shape(epar: float):
    """theta, f (Berryman 1980) and depolarisation factor L (Osborn 1945).

    EPAR = unique axis / degenerate axis: oblate < 1 < prolate. L is the
    depolarisation factor along the unique axis; the other two are (1-L)/2.
    A 13th-order series is used for |EPAR - 1| < 0.1.
    """
    epar = float(epar)
    if not (math.isfinite(epar) and epar > 0):
        raise ValueError(f"EPAR must be a positive finite number, got {epar}")
    e = epar - 1.0
    if abs(e) < 0.1:
        theta = _CT[-1]
        f = _CF[-1]
        for ct, cf in zip(reversed(_CT[:-1]), reversed(_CF[:-1])):
            theta = theta*e + ct
            f = f*e + cf
        return theta, f, 1.0 - theta
    if epar < 1.0:
        a, c = 1.0, epar
        theta = (a*a*c)/(a*a - c*c)**1.5*(math.acos(c/a) - (c/a)*math.sqrt(1 - c*c/(a*a)))
        f = c*c/(a*a - c*c)*(3*theta - 2)
        M = 1.0/epar
        L = M*M/(M*M - 1)*(1 - math.asin(math.sqrt(M*M - 1)/M)/math.sqrt(M*M - 1))
    else:
        a, c = epar, 1.0
        theta = (a*c*c)/(a*a - c*c)**1.5*((a/c)*math.sqrt(a*a/(c*c) - 1) - math.acosh(a/c))
        f = a*a/(a*a - c*c)*(2 - 3*theta)
        M = epar
        s = math.sqrt(M*M - 1)
        L = 1/(M*M - 1)*(M/(2*s)*math.log((M + s)/(M - s)) - 1)
    return theta, f, L


def PQ(k0, mu0, k2, mu2, epar):
    """Berryman's (1980) P and Q for a spheroidal inclusion (k2, mu2) in a host (k0, mu0)."""
    T, f, _ = shape(epar)
    k0 = np.asarray(k0, float); mu0 = np.asarray(mu0, float)
    nu0 = (3*k0 - 2*mu0)/(2*(3*k0 + mu0))
    A = mu2/mu0 - 1
    B = (k2/k0 - mu2/mu0)/3
    R = (1 - 2*nu0)/(2*(1 - nu0))
    F1 = 1 + A*(1.5*(f + T) - R*(1.5*f + 2.5*T - 4/3))
    F2 = (1 + A*(1 + 1.5*(f + T) - (R/2)*(3*f + 5*T)) + B*(3 - 4*R)
          + (A/2)*(A + 3*B)*(3 - 4*R)*(f + T - R*(f - T + 2*T*T)))
    F3 = 1 + A*(1 - (f + 1.5*T) + R*(f + T))
    F4 = 1 + (A/4)*(f + 3*T - R*(f - T))
    F5 = A*(-f + R*(f + T - 4/3)) + B*T*(3 - 4*R)
    F6 = 1 + A*(1 + f - R*(f + T)) + B*(1 - T)*(3 - 4*R)
    F7 = 2 + (A/4)*(3*f + 9*T - R*(3*f + 5*T)) + B*T*(3 - 4*R)
    F8 = A*(1 - 2*R + (f/2)*(R - 1) + (T/2)*(5*R - 3)) + B*(1 - T)*(3 - 4*R)
    F9 = A*((R - 1)*f - R*T) + B*T*(3 - 4*R)
    p = F1/F2
    q = (2/F3 + 1/F4 + (F4*F5 + F6*F7 - F8*F9)/(F2*F4))/5
    return p, q


def mbar(sig, sig2, epar):
    """m-bar of Cilli & Chapman (2021): (1/3) sig sum_j 1/(sig + L_j (sig2 - sig))."""
    _, _, L = shape(epar)
    sig = np.asarray(sig, float)
    return sig/3*(1/(sig - L*(sig - sig2)) + 4/(sig + sig2 + L*(sig - sig2)))


def R_berryman(sig, sig2, epar):
    """Berryman's (1995) R = m-bar/(3 sig)."""
    return mbar(sig, sig2, epar)/(3*np.asarray(sig, float))


def _eval_epar(epar: EPAR, phi: float) -> float:
    if callable(epar):
        a = float(epar(min(max(phi, 0.0), 1.0)))
        if not (math.isfinite(a) and a > 0):
            raise ValueError("EPAR(phi) must return a positive finite value; "
                             "a power law needs a floor, e.g. max(a*phi**b, 1e-3)")
        return a
    return float(epar)


# --------------------------------------------------------------------------
# Integration helper
# --------------------------------------------------------------------------
def _integrate(fun, targets, y0, rtol=1e-10, atol=1e-13, method="DOP853"):
    """Integrate dy/dx = fun(x, y) from x = 0 to every target; rows of the result."""
    targets = np.atleast_1d(np.asarray(targets, float)).ravel()
    y0 = np.atleast_1d(np.asarray(y0, float))
    Y = np.full((targets.size, y0.size), np.nan)
    Y[targets == 0] = y0
    for sgn in (1, -1):
        idx = np.flatnonzero(np.sign(targets) == sgn)
        if idx.size == 0:
            continue
        xs = np.unique(targets[idx])
        if sgn < 0:
            xs = xs[::-1]
        sol = solve_ivp(fun, (0.0, xs[-1]), y0, t_eval=xs, rtol=rtol, atol=atol,
                        method=method)
        if not sol.success or sol.t.size != xs.size:
            raise RuntimeError(f"ODE integration failed: {sol.message}")
        pos = {x: i for i, x in enumerate(xs)}
        Y[idx] = sol.y.T[[pos[x] for x in targets[idx]]]
    return Y


def _in_range(x, a, b):
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    return (x >= lo) & (x <= hi) & (x != b)


# --------------------------------------------------------------------------
# Cross-property DEM: elastic <-> electrical
# --------------------------------------------------------------------------
def forward(k1, mu1, sig1, k2, mu2, sig2, sig0, epar: EPAR, rtol=1e-10, atol=1e-13):
    """Bulk modulus, shear modulus and implied porosity from conductivity.

    Integrates dK/dsig = (K2-K) P/((sig2-sig) m-bar), dmu/dsig = (mu2-mu) Q/
    ((sig2-sig) m-bar) and dphi/dsig = (1-phi)/((sig2-sig) m-bar) from
    (sig1, k1, mu1, 0) to each sig0, in the scaled variables u = ln(sig/sig1),
    y = [K/k1, mu/k1, phi]. Returns arrays shaped like sig0 (NaN outside
    (sig1, sig2)). Parameters may be per-sample arrays; EPAR may be a scalar,
    a per-sample array or a function of porosity.
    """
    sig0 = np.asarray(sig0, float)
    shp = sig0.shape
    s0 = sig0.ravel()
    n = s0.size
    prm = [np.broadcast_to(np.asarray(v, float).ravel() if np.ndim(v) else v, (n,)).astype(float)
           for v in (k1, mu1, sig1, k2, mu2, sig2)]
    per_sample_epar = (not callable(epar)) and np.ndim(epar) > 0
    if per_sample_epar:
        ev = np.asarray(epar, float).ravel()
        if ev.size != n:
            raise ValueError("per-sample EPAR must have one value per measurement")
    ok = _in_range(s0, prm[2], prm[5])
    shared = all(np.ndim(v) == 0 for v in (k1, mu1, sig1, k2, mu2, sig2)) and not per_sample_epar
    k0 = np.full(n, np.nan); mu0 = np.full(n, np.nan); phi0 = np.full(n, np.nan)

    def make(k1_, k2_, mu2_, sig1_, sig2_, e):
        def deriv(u, y):
            sig = sig1_*math.exp(u)
            K, mu, phi = y[0]*k1_, y[1]*k1_, y[2]
            a = _eval_epar(e, phi)
            p, q = PQ(K, mu, k2_, mu2_, a)
            g = sig/((sig2_ - sig)*mbar(sig, sig2_, a))
            return [(k2_ - K)*p*g/k1_, (mu2_ - mu)*q*g/k1_, (1 - phi)*g]
        return deriv

    if shared:
        if ok.any():
            Y = _integrate(make(k1, k2, mu2, sig1, sig2, epar), np.log(s0[ok]/sig1),
                           [1.0, mu1/k1, 0.0], rtol, atol)
            k0[ok] = Y[:, 0]*k1; mu0[ok] = Y[:, 1]*k1; phi0[ok] = Y[:, 2]
    else:
        for i in np.flatnonzero(ok):
            e = ev[i] if per_sample_epar else epar
            Y = _integrate(make(prm[0][i], prm[3][i], prm[4][i], prm[2][i], prm[5][i], e),
                           math.log(s0[i]/prm[2][i]), [1.0, prm[1][i]/prm[0][i], 0.0], rtol, atol)
            k0[i] = Y[0, 0]*prm[0][i]; mu0[i] = Y[0, 1]*prm[0][i]; phi0[i] = Y[0, 2]
    return k0.reshape(shp), mu0.reshape(shp), phi0.reshape(shp)


def inverse(m0, which, k1, mu1, sig1, k2, mu2, sig2, epar: EPAR, rtol=1e-10, atol=1e-13):
    """Conductivity from a bulk ('K') or shear ('mu') modulus by direct integration.

    Returns (sig0, phi0, other modulus) shaped like m0; NaN where the modulus
    lies outside the interval between the background and inclusion values.
    Scalar parameters only (loop for per-sample parameters).
    """
    use_k = str(which).lower() in ("k", "bulk", "0")
    m0 = np.asarray(m0, float)
    shp = m0.shape
    m = m0.ravel()
    m1, m2 = (k1, k2) if use_k else (mu1, mu2)
    ok = _in_range(m, m1, m2)
    sig0 = np.full(m.size, np.nan); phi0 = np.full(m.size, np.nan); other = np.full(m.size, np.nan)

    def deriv(x, y):
        sig = sig1*math.exp(y[0])
        phi = y[2]
        a = _eval_epar(epar, phi)
        if use_k:
            K, mu = k1 + x*k1, y[1]*k1
            p, q = PQ(K, mu, k2, mu2, a)
            h = k1/((k2 - K)*p)
            d_other = (mu2 - mu)*q/k1*h
        else:
            mu, K = mu1 + x*k1, y[1]*k1
            p, q = PQ(K, mu, k2, mu2, a)
            h = k1/((mu2 - mu)*q)
            d_other = (k2 - K)*p/k1*h
        return [(sig2 - sig)*mbar(sig, sig2, a)/sig*h, d_other, (1 - phi)*h]

    if ok.any():
        y0 = [0.0, mu1/k1, 0.0] if use_k else [0.0, 1.0, 0.0]
        Y = _integrate(deriv, (m[ok] - m1)/k1, y0, rtol, atol)
        sig0[ok] = sig1*np.exp(Y[:, 0]); other[ok] = Y[:, 1]*k1; phi0[ok] = Y[:, 2]
    return sig0.reshape(shp), phi0.reshape(shp), other.reshape(shp)


# --------------------------------------------------------------------------
# Single-property DEMs and diagnostics
# --------------------------------------------------------------------------
def elec_dem(phi, sig1, sig2, epar: EPAR, rtol=1e-10, atol=1e-13):
    """Electrical DEM: (1-phi) dsig/dphi = (sig2 - sig) m-bar, sig(0) = sig1."""
    phi = np.asarray(phi, float)
    def deriv(p, y):
        s = sig1*math.exp(y[0])
        return [(sig2 - s)*mbar(s, sig2, _eval_epar(epar, p))/(s*(1 - p))]
    Y = _integrate(deriv, phi.ravel(), [0.0], rtol, atol)
    return (sig1*np.exp(Y[:, 0])).reshape(phi.shape)


def elas_dem(phi, k1, mu1, k2, mu2, epar: EPAR, rtol=1e-10, atol=1e-13):
    """Elastic DEM of Berryman (1992): (1-phi) dK/dphi = (K2-K) P, (1-phi) dmu/dphi = (mu2-mu) Q."""
    phi = np.asarray(phi, float)
    def deriv(p, y):
        K, mu = y[0]*k1, y[1]*k1
        P, Q = PQ(K, mu, k2, mu2, _eval_epar(epar, p))
        return [(k2 - K)*P/(k1*(1 - p)), (mu2 - mu)*Q/(k1*(1 - p))]
    Y = _integrate(deriv, phi.ravel(), [1.0, mu1/k1], rtol, atol)
    return (Y[:, 0]*k1).reshape(phi.shape), (Y[:, 1]*k1).reshape(phi.shape)


def implied_porosity(sig0, sig1, sig2, epar: EPAR, rtol=1e-10, atol=1e-13):
    """Porosity the electrical DEM implies for each conductivity."""
    sig0 = np.asarray(sig0, float)
    s = sig0.ravel()
    ok = _in_range(s, sig1, sig2)
    phi = np.full(s.size, np.nan)
    def deriv(u, y):
        sig = sig1*math.exp(u)
        return [sig*(1 - y[0])/((sig2 - sig)*mbar(sig, sig2, _eval_epar(epar, y[0])))]
    if ok.any():
        phi[ok] = _integrate(deriv, np.log(s[ok]/sig1), [0.0], rtol, atol)[:, 0]
    return phi.reshape(sig0.shape)


def forward_conduct(p1, sig1, p2, sig2, sig0, epar: EPAR, rtol=1e-10, atol=1e-13):
    """Cross-property DEM between two conductivity-like properties (e.g. thermal from electrical)."""
    sig0 = np.asarray(sig0, float)
    s = sig0.ravel()
    ok = _in_range(s, sig1, sig2)
    p0 = np.full(s.size, np.nan); phi0 = np.full(s.size, np.nan)
    def deriv(u, y):
        sig = sig1*math.exp(u)
        p, phi = p1*math.exp(y[0]), y[1]
        a = _eval_epar(epar, phi)
        g = sig/((sig2 - sig)*mbar(sig, sig2, a))
        return [(p2 - p)*mbar(p, p2, a)/p*g, (1 - phi)*g]
    if ok.any():
        Y = _integrate(deriv, np.log(s[ok]/sig1), [0.0, 0.0], rtol, atol)
        p0[ok] = p1*np.exp(Y[:, 0]); phi0[ok] = Y[:, 1]
    return p0.reshape(sig0.shape), phi0.reshape(sig0.shape)


def elec_elas_hs(sig, sig_sol, sig_fl, k_sol, k_fl, mu_sol):
    """Electrical-elastic Hashin-Shtrikman bounds (Carcione et al. 2007).

    Returns (k_up, k_low, mu_up, mu_low, phi_low, phi_up).
    """
    if not sig_fl > sig_sol:
        raise ValueError("requires sig_fl > sig_sol")
    sig = np.asarray(sig, float)
    phi_low = (sig_sol - sig)/(sig_sol - sig_fl)*3*sig_fl/(sig + 2*sig_fl)
    phi_up = (sig_sol - sig)/(sig_sol - sig_fl)*(sig_fl + 2*sig_sol)/(sig + 2*sig_sol)
    x = 4/3*mu_sol
    k_up = 1/((1 - phi_low)/(k_sol + x) + phi_low/(k_fl + x)) - x
    k_low = 1/(phi_up/k_fl + (1 - phi_up)/k_sol)
    z = mu_sol/6*(9*k_sol + 8*mu_sol)/(k_sol + 2*mu_sol)
    mu_up = 1/((1 - phi_low)/(mu_sol + z) + phi_low/z) - z
    return k_up, k_low, mu_up, np.zeros_like(sig), phi_low, phi_up


# --------------------------------------------------------------------------
# Calibration and velocities
# --------------------------------------------------------------------------
def calibrate_epar(k0, mu0, sig0, k1, mu1, sig1, k2, mu2, sig2, bounds=(1.0, 100.0),
                   norm=2, relative=False, weights=None, log_search=True, xtol=1e-6):
    """Calibrate alpha_K and alpha_mu (either k0 or mu0 may be None). Returns a dict."""
    sig0 = np.asarray(sig0, float).ravel()
    w = np.ones_like(sig0) if weights is None else np.asarray(weights, float).ravel()
    lb, ub = bounds
    out = {"eparK": np.nan, "eparMu": np.nan, "misfitK": np.nan, "misfitMu": np.nan,
           "kModel": None, "muModel": None, "r2K": np.nan, "r2Mu": np.nan}

    def solve(meas, use_k):
        meas = np.asarray(meas, float).ravel()
        def model(e):
            kk, mm, _ = forward(k1, mu1, sig1, k2, mu2, sig2, sig0, e)
            return kk if use_k else mm
        def cost(e):
            r = meas - model(e)
            if relative:
                r = r/meas
            r = w*r
            return np.abs(r).sum() if norm == 1 else math.sqrt((r*r).sum())
        if log_search:
            res = minimize_scalar(lambda lx: cost(math.exp(lx)), bounds=(math.log(lb), math.log(ub)),
                                  method="bounded", options={"xatol": xtol})
            e = math.exp(res.x)
        else:
            res = minimize_scalar(cost, bounds=(lb, ub), method="bounded", options={"xatol": xtol})
            e = res.x
        mdl = model(e)
        sst = ((meas - meas.mean())**2).sum()
        r2 = 1 - ((meas - mdl)**2).sum()/sst if sst > 0 else np.nan
        return e, res.fun, mdl, r2

    if k0 is not None:
        out["eparK"], out["misfitK"], out["kModel"], out["r2K"] = solve(k0, True)
    if mu0 is not None:
        out["eparMu"], out["misfitMu"], out["muModel"], out["r2Mu"] = solve(mu0, False)
    return out


def model_vpvs(k1, mu1, sig1, k2, mu2, sig2, sig0, kepar, muepar, density=None):
    """Vp, Vs and Vp/Vs from conductivity (bulk modulus with alpha_K, shear with alpha_mu).

    density: None -> Gardner (rho = 310 Vp^0.25 in SI), a float/array of measured
    densities, or a tuple (rho1, rho2) for rho = (1-phi) rho1 + phi rho2 with the
    mean implied porosity of the two paths. Returns a dict.
    """
    k, _, phiK = forward(k1, mu1, sig1, k2, mu2, sig2, sig0, kepar)
    _, mu, phiMu = forward(k1, mu1, sig1, k2, mu2, sig2, sig0, muepar)
    M = k + 4/3*mu
    if density is None:
        a, b = 310.0, 0.25
        vp = (M/a)**(1/(2 + b))
        rho = a*vp**b
    elif isinstance(density, tuple):
        phi = (phiK + phiMu)/2
        rho = (1 - phi)*density[0] + phi*density[1]
        vp = np.sqrt(M/rho)
    else:
        rho = np.broadcast_to(np.asarray(density, float), np.shape(k)).astype(float)
        vp = np.sqrt(M/rho)
    vs = np.sqrt(mu/rho)
    return {"k": k, "mu": mu, "phiK": phiK, "phiMu": phiMu, "rho": rho, "vp": vp,
            "vs": vs, "vpvs": vp/vs, "M": M}
