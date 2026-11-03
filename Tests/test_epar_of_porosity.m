function test_epar_of_porosity()
% EPAR given as a function of porosity (power law of Cilli & Chapman 2020)
% is honoured consistently by the forward, inverse and single-property DEMs.

p = cpdem_TestParams();
law = @(phi) max(30*phi.^0.3, 0.05);   % floored power law (see cpdem_EvalEPAR)

[k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, law);
phi = cpdem_ImpliedPorosity(p.sig0, p.sig1, p.sig2, law);
cpdem_AssertClose(phi, phi0, 1e-7, 1e-10, 'implied porosity with EPAR(phi)');

[k, mu] = cpdem_ElasDEM(phi, p.k1, p.mu1, p.k2, p.mu2, law);
cpdem_AssertClose([k mu], [k0 mu0], 1e-7, 0, 'elastic DEM with EPAR(phi)');

sig = cpdem_ElecDEM(phi, p.sig1, p.sig2, law);
cpdem_AssertClose(sig, p.sig0, 1e-7, 0, 'electrical DEM with EPAR(phi)');

s = cpdem_Inverse(k0, 'K', p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, law);
cpdem_AssertClose(s, p.sig0, 1e-7, 0, 'inverse with EPAR(phi)');

% A constant law equals a constant EPAR
[kc, muc] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, @(phi) 16.4);
[kk, mm] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 16.4);
cpdem_AssertClose([kc muc], [kk mm], 1e-12, 0, 'constant law');

% Bad law is rejected
threw = false;
try
    cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 0.1, @(phi) -1);
catch
    threw = true;
end
assert(threw, 'negative EPAR from a function handle must throw');

end
