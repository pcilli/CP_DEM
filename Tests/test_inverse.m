function test_inverse()
% Direct ODE inversion recovers the conductivity from either modulus to
% solver tolerance, carries the other modulus and porosity along the path,
% and handles out-of-range measurements. The v1.0 wrapper agrees.

p = cpdem_TestParams();
[k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK);

[sK, phiK, muK] = cpdem_Inverse(k0, 'K', p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK);
cpdem_AssertClose(sK, p.sig0, 1e-7, 0, 'sigma_K round trip');
cpdem_AssertClose(muK, mu0, 1e-7, 0, 'shear modulus along K path');
cpdem_AssertClose(phiK, phi0, 1e-7, 1e-10, 'porosity along K path');

[sMu, phiMu, kMu] = cpdem_Inverse(mu0, 'mu', p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK);
cpdem_AssertClose(sMu, p.sig0, 1e-7, 0, 'sigma_mu round trip');
cpdem_AssertClose(kMu, k0, 1e-7, 0, 'bulk modulus along mu path');

% v1.0 modFlag convention and shape preservation
s = cpdem_Inverse(reshape(k0, [], 1), 0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK);
assert(iscolumn(s), 'shape preserved');
cpdem_AssertClose(s(:), p.sig0(:), 1e-7, 0, 'modFlag 0');

% Out-of-range moduli -> NaN with a warning; in range unaffected
ws = warning('off', 'CPDEM:outOfRange');
s = cpdem_Inverse([p.k1*1.01, k0(2), p.k2*0.5], 'K', p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK);
warning(ws);
assert(isnan(s(1)) && isnan(s(3)), 'NaN for unreachable moduli');
cpdem_AssertClose(s(2), p.sig0(2), 1e-7, 0, 'in-range sample unaffected');

% v1.0 wrapper: same answers, bounds honoured
sw = cpdem_InvertConduct(k0, mu0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK, p.sig1, p.sig2, 0, 2, 0);
cpdem_AssertClose(sw, p.sig0, 1e-7, 0, 'cpdem_InvertConduct K');
sw = cpdem_InvertConduct(k0, mu0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK, p.sig1, p.sig2, 1, 1, 0);
cpdem_AssertClose(sw, p.sig0, 1e-7, 0, 'cpdem_InvertConduct mu');
sw = cpdem_InvertConduct([p.k1*1.1 p.k2*0.5], [1 1], p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK, p.sig1, p.sig2, 0, 2, 0);
cpdem_AssertClose(sw, [p.sig1 p.sig2], 0, 0, 'wrapper clips unreachable moduli to the bounds');

% Per-sample parameters (v1.0 paramFlag = 1 style) round trip
k1v = p.k1*[1 1.02 0.98 1 1.01];
[k0v, ~] = cpdem_Forward(k1v, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK);
sv = cpdem_Inverse(k0v, 'K', k1v, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.eparK);
cpdem_AssertClose(sv, p.sig0, 1e-7, 0, 'per-sample parameters round trip');

end
