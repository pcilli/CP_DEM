function test_anisotropic_conductivity()
% Anisotropic inclusions and aligned inclusions: reduction to the isotropic
% model, the exact laminate and fibre limits, and consistency between the
% volume-fraction DEM and the cross-property integration.

p = cpdem_TestParams();
phi = [0.05 0.1 0.2 0.3 0.4];

% --- Anisotropic inclusion material, randomly oriented ------------------
% Equal components reduce to the isotropic form (and to m-bar)
for epar = [0.1 1 16.4]
    s = 0.3;
    cpdem_AssertClose(cpdem_CondIncrement(s, [5 5], epar), (5 - s)*mbar(s, 5, epar), 1e-14, 0, 'isotropic reduction');
    cpdem_AssertClose(cpdem_ElecDEM(phi, p.sig1, [p.sig2 p.sig2], epar), ...
                      cpdem_ElecDEM(phi, p.sig1, p.sig2, epar), 1e-10, 0, 'ElecDEM anisotropic = isotropic');
end
% An anisotropic inclusion gives a result between the two isotropic ones
sA = cpdem_ElecDEM(phi, 1, [100 1e-3], 0.1);
sLo = cpdem_ElecDEM(phi, 1, 1e-3, 0.1);
sHi = cpdem_ElecDEM(phi, 1, 100, 0.1);
assert(all(sA > sLo & sA < sHi), 'anisotropic inclusion between isotropic end members');
% Forward conduct and implied porosity accept anisotropic inclusions
[lam, ph] = cpdem_ForwardConduct(7.7, p.sig1, [0.6 0.6], p.sig2, p.sig0, 16.4);
[lam0, ph0] = cpdem_ForwardConduct(7.7, p.sig1, 0.6, p.sig2, p.sig0, 16.4);
cpdem_AssertClose([lam ph], [lam0 ph0], 1e-10, 0, 'ForwardConduct anisotropic = isotropic');
phA = cpdem_ImpliedPorosity(sA, 1, [100 1e-3], 0.1);
cpdem_AssertClose(phA, phi, 1e-7, 1e-10, 'implied fraction with anisotropic inclusion');

% --- Aligned inclusions --------------------------------------------------
sig1 = 1; sig2 = 20;
% Spheres: both components equal the isotropic Bruggeman closed form
[spar, sperp] = cpdem_AlignedCondDEM(phi, sig1, sig2, 1);
cpdem_AssertClose(spar, sperp, 1e-10, 0, 'aligned spheres isotropic');
cpdem_AssertClose((spar - sig2)./(sig1 - sig2).*(sig1./spar).^(1/3), 1 - phi, 1e-8, 0, 'aligned spheres Bruggeman');

% Thin aligned disks: series average along the axis, parallel in the plane
[spar, sperp] = cpdem_AlignedCondDEM(phi, sig1, sig2, 1e-7);
cpdem_AssertClose(spar, 1./((1 - phi)/sig1 + phi/sig2), 1e-5, 0, 'aligned disks: series');
cpdem_AssertClose(sperp, (1 - phi)*sig1 + phi*sig2, 1e-5, 0, 'aligned disks: parallel');

% Long aligned needles: parallel average along the axis, 2-D Bruggeman in
% the plane, (sig - sig2)/(sig1 - sig2) (sig1/sig)^(1/2) = 1 - phi
[spar, sperp] = cpdem_AlignedCondDEM(phi, sig1, sig2, 1e7);
cpdem_AssertClose(spar, (1 - phi)*sig1 + phi*sig2, 1e-5, 0, 'aligned needles: parallel');
cpdem_AssertClose((sperp - sig2)./(sig1 - sig2).*(sig1./sperp).^(1/2), 1 - phi, 1e-5, 0, 'aligned needles: 2-D Bruggeman');

% Anisotropic host and inclusion are accepted; ordering is physical
[spar, sperp] = cpdem_AlignedCondDEM(phi, [1 2], [50 0.5], 0.2);
assert(all(isfinite([spar sperp])) && all(sperp > 2*(1 - phi)) && all(spar > 0), 'TI host, anisotropic inclusion');

% --- Aligned cross-property -------------------------------------------
lam1 = 7.7; lam2 = 0.6; s1 = p.sig1; s2 = p.sig2;
for epar = [0.1 16.4]
    [spar, sperp] = cpdem_AlignedCondDEM(phi, s1, s2, epar);
    [lpar, lperp] = cpdem_AlignedCondDEM(phi, lam1, lam2, epar);
    o = cpdem_AlignedForwardConduct(lam1, s1, lam2, s2, spar, epar, 'par');
    cpdem_AssertClose([o.sigPerp o.pPar o.pPerp o.phi], [sperp lpar lperp phi], 1e-6, 1e-9, sprintf('aligned cross-property from par, EPAR %g', epar));
    o = cpdem_AlignedForwardConduct(lam1, s1, lam2, s2, sperp, epar, 'perp');
    cpdem_AssertClose([o.sigPar o.pPar o.pPerp o.phi], [spar lpar lperp phi], 1e-6, 1e-9, sprintf('aligned cross-property from perp, EPAR %g', epar));
end
% With spheres the aligned cross-property model equals the isotropic one
[spar, ~] = cpdem_AlignedCondDEM(phi, s1, s2, 1);
o = cpdem_AlignedForwardConduct(lam1, s1, lam2, s2, spar, 1, 'par');
cpdem_AssertClose(o.pPar, cpdem_ForwardConduct(lam1, s1, lam2, s2, spar, 1), 1e-7, 0, 'spheres: aligned = isotropic');

end
