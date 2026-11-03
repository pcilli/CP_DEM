function test_forward_composition()
% The cross-property DEM must equal the composition of the electrical DEM
% (conductivity -> porosity) and the elastic DEM (porosity -> moduli), for
% prolate, spherical and oblate pores, and the single integration of a
% whole data set must match sample-by-sample integration.

p = cpdem_TestParams();
for epar = [16.4 1 0.1]
    [k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, epar);

    phi = cpdem_ImpliedPorosity(p.sig0, p.sig1, p.sig2, epar);
    cpdem_AssertClose(phi, phi0, 1e-7, 1e-9, sprintf('implied porosity, EPAR %g', epar));

    [k, mu] = cpdem_ElasDEM(phi, p.k1, p.mu1, p.k2, p.mu2, epar);
    cpdem_AssertClose(k, k0, 1e-7, 0, sprintf('bulk modulus composition, EPAR %g', epar));
    cpdem_AssertClose(mu, mu0, 1e-7, 0, sprintf('shear modulus composition, EPAR %g', epar));

    sig = cpdem_ElecDEM(phi, p.sig1, p.sig2, epar);
    cpdem_AssertClose(sig, p.sig0, 1e-7, 0, sprintf('electrical DEM round trip, EPAR %g', epar));

    % Sample-by-sample integration vs single integration
    for ii = 1:numel(p.sig0)
        [ks, ms, ps] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0(ii), epar);
        cpdem_AssertClose([ks ms ps], [k0(ii) mu0(ii) phi0(ii)], 1e-8, 0, 'single vs batch');
    end

    % Monotonic in conductivity (the v1.0 solver noise broke this)
    assert(all(diff(k0) < 0) && all(diff(mu0) < 0) && all(diff(phi0) > 0), ...
        'moduli must decrease and porosity increase with conductivity');
end

% Equations 18 & 19 (m-bar form) and 20 & 21 (R form) of the paper give the
% same result
[kA, muA, phiA] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 16.4, struct('Form', 'mbar'));
[kB, muB, phiB] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 16.4, struct('Form', 'R'));
cpdem_AssertClose([kB muB phiB], [kA muA phiA], 1e-12, 0, 'm-bar and R forms');
threw = false;
try
    cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 0.1, 16.4, struct('Form', 'x'));
catch
    threw = true;
end
assert(threw, 'unknown Form must throw');

% Output keeps the shape of sig0, and sig0 == sig1 returns the background
[k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, [p.sig1 0.1; 0.2 0.3], 16.4);
assert(isequal(size(k0), [2 2]) && isequal(size(mu0), [2 2]) && isequal(size(phi0), [2 2]), 'shape');
cpdem_AssertClose([k0(1,1) mu0(1,1) phi0(1,1)], [p.k1 p.mu1 0], 0, 0, 'background at sig1');

% Bruggeman closed form for spheres: (sig - sig2)/(sig1 - sig2) (sig1/sig)^(1/3) = 1 - phi
phi = [0.05 0.1 0.2 0.3];
s = cpdem_ElecDEM(phi, p.sig1, p.sig2, 1);
cpdem_AssertClose((s - p.sig2)./(p.sig1 - p.sig2).*(p.sig1./s).^(1/3), 1 - phi, 1e-8, 0, 'Bruggeman closed form');

% Archie's law for insulating inclusions in a conducting host (Sen et al. 1981)
m = mbar(1, 0, 0.2);
s = cpdem_ElecDEM(phi, 1, 1e-300, 0.2);   % (sig2 = 0 is not allowed; use a negligible value)
cpdem_AssertClose(s, (1 - phi).^m, 1e-6, 0, 'Archie limit');

end
