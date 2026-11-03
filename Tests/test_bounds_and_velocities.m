function test_bounds_and_velocities()
% The electrical-elastic Hashin-Shtrikman bounds contain the DEM
% predictions and the implied porosity, and the Vp/Vs model is consistent
% with the moduli it is built from.

p = cpdem_TestParams();
sig = linspace(0.005, 2, 25);
[kUp, kLow, muUp, muLow, phiLow, phiUp] = cpdem_ElecElasHS(sig, p.sig1, p.sig2, p.k1, p.k2, p.mu1);

assert(all(phiLow <= phiUp) && all(phiLow >= 0) && all(phiUp <= 1), 'porosity bounds ordered');
assert(all(kLow <= kUp) && all(muLow <= muUp) && all(muLow == 0), 'elastic bounds ordered');

for epar = [0.05 1 16.4]
    [k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, sig, epar);
    tol = 1e-9;
    assert(all(k0 <= kUp*(1 + tol)) && all(k0 >= kLow*(1 - tol)), sprintf('K within HS bounds, EPAR %g', epar));
    assert(all(mu0 <= muUp*(1 + tol)) && all(mu0 >= muLow), sprintf('mu within HS bounds, EPAR %g', epar));
    assert(all(phi0 >= phiLow*(1 - tol)) && all(phi0 <= phiUp*(1 + tol)), sprintf('porosity within bounds, EPAR %g', epar));
end

% Bounds at the end members
[kUp, kLow, ~, ~, phiLow, phiUp] = cpdem_ElecElasHS(p.sig1, p.sig1, p.sig2, p.k1, p.k2, p.mu1);
cpdem_AssertClose([kUp kLow phiLow phiUp], [p.k1 p.k1 0 0], 1e-12, 0, 'bounds at the mineral end');

% Fluid less conductive than solid is rejected
threw = false;
try
    cpdem_ElecElasHS(0.1, 1, 0.01, p.k1, p.k2, p.mu1);
catch
    threw = true;
end
assert(threw, 'sigFl <= sigSol must throw');

% Vp/Vs model: moduli match cpdem_Forward with each EPAR, velocities match
% the density used, Gardner is self-consistent.
out = cpdem_ModelVpVs(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK, p.eparMu, ...
    struct('method', 'implied', 'rho1', 2650, 'rho2', 1050));
[k0, ~, phiK] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK);
[~, mu0, phiMu] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparMu);
cpdem_AssertClose(out.k, k0, 1e-12, 0, 'VpVs bulk modulus');
cpdem_AssertClose(out.mu, mu0, 1e-12, 0, 'VpVs shear modulus');
rho = (1 - (phiK + phiMu)/2)*2650 + (phiK + phiMu)/2*1050;
cpdem_AssertClose(out.rho, rho, 1e-12, 0, 'implied density');
cpdem_AssertClose(out.vp, sqrt((k0 + 4/3*mu0)./rho), 1e-12, 0, 'Vp');
cpdem_AssertClose(out.vs, sqrt(mu0./rho), 1e-12, 0, 'Vs');
cpdem_AssertClose(out.vpvs, out.vp./out.vs, 1e-12, 0, 'Vp/Vs');
assert(all(out.vpvs > 1.4 & out.vpvs < 2.5), 'plausible Vp/Vs for sandstone');

outG = cpdem_ModelVpVs(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK, p.eparMu);
cpdem_AssertClose(outG.rho, 310*outG.vp.^0.25, 1e-12, 0, 'Gardner density');
cpdem_AssertClose(outG.rho.*outG.vp.^2, outG.M, 1e-12, 0, 'Gardner Vp consistent with M');

outM = cpdem_ModelVpVs(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK, p.eparMu, ...
    struct('method', 'measured', 'rho', 2300));
cpdem_AssertClose(outM.vs, sqrt(mu0/2300), 1e-12, 0, 'measured density');

end
