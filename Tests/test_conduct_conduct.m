function test_conduct_conduct()
% Thermal-electrical cross-property DEM is consistent with the electrical
% DEM evaluated for each property at the common implied porosity.

p = cpdem_TestParams();
lam1 = 7.7;  % quartz thermal conductivity, W/m/K
lam2 = 0.6;  % water
for epar = [0.2 1 16.4]
    [lam, phi] = cpdem_ForwardConduct(lam1, p.sig1, lam2, p.sig2, p.sig0, epar);
    phiE = cpdem_ImpliedPorosity(p.sig0, p.sig1, p.sig2, epar);
    cpdem_AssertClose(phi, phiE, 1e-7, 1e-10, 'common porosity');
    lamE = cpdem_ElecDEM(phi, lam1, lam2, epar);
    cpdem_AssertClose(lam, lamE, 1e-7, 0, sprintf('thermal conductivity, EPAR %g', epar));
    assert(all(lam < lam1 & lam > lam2) && all(diff(lam) < 0), 'thermal conductivity decreases with porosity');
end

% Swapping the roles of the two properties inverts the relation
[lam, ~] = cpdem_ForwardConduct(lam1, p.sig1, lam2, p.sig2, p.sig0, 16.4);
[sig, ~] = cpdem_ForwardConduct(p.sig1, lam1, p.sig2, lam2, lam, 16.4);
cpdem_AssertClose(sig, p.sig0, 1e-6, 0, 'round trip');

end
