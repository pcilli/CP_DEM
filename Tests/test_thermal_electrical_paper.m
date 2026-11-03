function test_thermal_electrical_paper()
% The thermal-electrical cross-property DEM of Cilli & Chapman
% (arXiv:2102.12434): equations 5-9, the mirror symmetry of the forward
% and inverse models, the end-member limits, the Hashin-Shtrikman bounds
% and the three-parameter calibration of Section 3, using the paper's
% copper-graphite parameters.

c = cpdem_PublishedCuGraphite();
k1 = c.kappa1; s1 = c.sig1; k2 = c.kappa2; s2 = c.sig2;

% --- Equations 7-9: closed-form right-hand sides for spheres, disks and
% needles equal the general right-hand side (equation 5) at those EPARs
sig = s1*[0.9 0.5 0.2 0.05];
kap = k1*[0.9 0.6 0.3 0.1];
rhs = @(k, s, a) cpdem_CondIncrement(k, k2, a)./cpdem_CondIncrement(s, s2, a);
eq7 = kap.*(k2 - kap).*(2*sig + s2)./(sig.*(s2 - sig).*(2*kap + k2));
eq8 = s2*(k2 - kap).*(2*k2 + kap)./(k2*(s2 - sig).*(2*s2 + sig));
eq9 = (s2 + sig).*(k2 - kap).*(k2 + 5*kap)./((k2 + kap).*(s2 - sig).*(s2 + 5*sig));
cpdem_AssertClose(rhs(kap, sig, 1),   eq7, 1e-12, 0, 'equation 7 (spheres)');
cpdem_AssertClose(rhs(kap, sig, 0),   eq8, 1e-12, 0, 'equation 8 (disks)');
cpdem_AssertClose(rhs(kap, sig, Inf), eq9, 1e-12, 0, 'equation 9 (needles)');

% --- Equation 2: m-bar as printed in the paper
L = 1 - 2/3;  % sphere
for a = [0.097 1 5]
    [~, ~, L] = cpdem_Shape(a);
    s = 0.3; 
    paper = (1/3)*s*(4/(s + s2 + L*(s - s2)) + 1/(s - L*(s - s2)));
    cpdem_AssertClose(mbar(s, s2, a), paper, 1e-14, 0, sprintf('equation 2 at EPAR %g', a));
end

% --- Equation 5 (thermal from electrical) and equation 6 (electrical from
% thermal) are inverses of each other; curves of Figure 1a and 1b mirror
sig0 = s1*[0.95 0.7 0.5 0.3 0.1 0.03];
for a = c.alphaCurves
    kap0 = cpdem_ForwardConduct(k1, s1, k2, s2, sig0, a);
    back = cpdem_ForwardConduct(s1, k1, s2, k2, kap0, a);
    cpdem_AssertClose(back, sig0, 1e-7, 0, sprintf('equations 5 and 6 are reciprocal, EPAR %g', a));
    assert(all(kap0 < k1 & kap0 > k2) && all(diff(kap0) < 0), 'thermal conductivity between end members and monotonic');
end

% --- End-member limits: the model reaches (sig2, kappa2) as phi -> 1
kapNear = cpdem_ForwardConduct(k1, s1, k2, s2, s2*(1 + 1e-6), c.alpha);
cpdem_AssertClose(kapNear, k2, 1e-3, 0, 'high volume fraction limit');
cpdem_AssertClose(cpdem_ForwardConduct(k1, s1, k2, s2, s1, c.alpha), k1, 0, 0, 'low volume fraction limit');

% --- Hashin-Shtrikman bounds (Carcione stratagem) contain every modelled
% trend, and the disk trend lies on the bound built on the inclusion phase,
% the single-property lower bound (Section 4, Norris 1985)
[kUp, kLow, phiMin, phiMax, kOn1, kOn2] = cpdem_CondCondHS(sig0, s1, s2, k1, k2);
assert(all(phiMin <= phiMax) && all(phiMin >= 0) && all(phiMax <= 1), 'volume-fraction interval');
for a = c.alphaCurves
    kap0 = cpdem_ForwardConduct(k1, s1, k2, s2, sig0, a);
    assert(all(kap0 <= kUp*(1 + 1e-9)) && all(kap0 >= kLow*(1 - 1e-9)), sprintf('within HS bounds, EPAR %g', a));
end
kapDisk = cpdem_ForwardConduct(k1, s1, k2, s2, sig0, 0);
cpdem_AssertClose(kapDisk, kOn2, 1e-8, 0, 'disks on the Hashin-Shtrikman curve built on the inclusion phase');
% The single-property disk DEM equals the single-property bound built on
% the inclusion phase (Norris 1985) for both properties
ph = [0.1 0.3 0.5];
cpdem_AssertClose(cpdem_ElecDEM(ph, s1, s2, 0), s2 + (1 - ph)./(1./(s1 - s2) + ph./(3*s2)), 1e-8, 0, 'disk DEM = HS on phase 2, sigma');
cpdem_AssertClose(cpdem_ElecDEM(ph, k1, k2, 0), k2 + (1 - ph)./(1./(k1 - k2) + ph./(3*k2)), 1e-8, 0, 'disk DEM = HS on phase 2, kappa');
% Swapping the properties bounds sigma given kappa
[sUp, sLow] = cpdem_CondCondHS(kapDisk, k1, k2, s1, s2);
assert(all(sig0 <= sUp*(1 + 1e-9)) && all(sig0 >= sLow*(1 - 1e-9)), 'bounds in the other direction');
% Agreement with the electrical-elastic routine's porosity bounds when the
% inclusion is the more conductive phase
[~, ~, ~, ~, pL, pU] = cpdem_ElecElasHS([0.02 0.1], 1e-4, 5, 37e9, 2.25e9, 44e9);
[~, ~, qMin, qMax] = cpdem_CondCondHS([0.02 0.1], 1e-4, 5, 1, 2);
cpdem_AssertClose([qMin qMax], [pL pU], 1e-12, 0, 'porosity bounds agree with cpdem_ElecElasHS');

% --- Section 3 calibration: three parameters recovered from synthetic
% "measurements" at the paper's volume fractions (phi <= 0.4) generated
% with the paper's solved parameters
phi = c.phiCalibration(c.phiCalibration > 0);
sigMeas = cpdem_ElecDEM(phi, s1, s2, c.alpha);
kapMeas = cpdem_ElecDEM(phi, k1, k2, c.alpha);
cal = cpdem_CalibrateConduct(kapMeas, sigMeas, k1, s1, ...
    [c.kappaFlakePlane c.kappaFlakeAxis], [c.sigFlakeAxis c.sigFlakePlane], 'Refine', 1);
cpdem_AssertClose([cal.epar cal.p2 cal.sig2], [c.alpha k2 s2], 2e-2, 0, 'three-parameter calibration');
assert(cal.r2 > 0.9999, 'calibration fit');
% Bounds are honoured when they exclude the truth
cal2 = cpdem_CalibrateConduct(kapMeas, sigMeas, k1, s1, [30 50], c.sig2, 'GridSize', [5 3 1], 'Refine', 1, 'Limits', false, 'TolX', 1e-3);
assert(cal2.p2 >= 30 && cal2.p2 <= 50 && cal2.sig2 == c.sig2, 'bounds honoured');

end
