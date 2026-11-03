function test_calibrate_epar()
% EPAR calibration recovers the EPARs used to generate synthetic data,
% through both the new and the v1.0 interface, for shared and per-sample
% calibration and for both norms.

p = cpdem_TestParams();
[k0, ~] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparK);
[~, mu0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, p.eparMu);

cal = cpdem_CalibrateEPAR(k0, mu0, p.sig0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2);
cpdem_AssertClose([cal.eparK cal.eparMu], [p.eparK p.eparMu], 1e-4, 0, 'calibrated EPARs');
assert(cal.r2K > 0.999999 && cal.r2Mu > 0.999999, 'R^2 of exact data');
cpdem_AssertClose(cal.kModel(:), k0(:), 1e-5, 0, 'modelled K');
cpdem_AssertClose(cal.muModel(:), mu0(:), 1e-5, 0, 'modelled mu');

% One modulus only, l1 norm, linear search, custom bounds via name/value
cal = cpdem_CalibrateEPAR(k0, [], p.sig0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, ...
    cpdem_Options('Norm', 1, 'LogSearch', false, 'Bounds', [5 40]));
cpdem_AssertClose(cal.eparK, p.eparK, 1e-4, 0, 'l1, linear search');
assert(isnan(cal.eparMu) && isempty(cal.muModel), 'skipped shear modulus');

% Per-sample calibration (each sample generated with its own EPAR)
eparS = [10 12 14 16 18];
[kS, ~] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, eparS);
cal = cpdem_CalibrateEPAR(kS, [], p.sig0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, ...
    struct('PerSample', true, 'Bounds', [2 50]));
cpdem_AssertClose(cal.eparK(:), eparS(:), 1e-4, 0, 'per-sample EPARs');

% v1.0 interface
e = cpdem_SolveEPAR(k0, mu0, p.sig0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 0, 2, 0, 2, 50);
cpdem_AssertClose(e, p.eparK, 1e-4, 0, 'cpdem_SolveEPAR K');
[e, misfit] = cpdem_SolveEPAR(k0, mu0, p.sig0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 1, 1, 0, 2, 50);
cpdem_AssertClose(e, p.eparMu, 1e-4, 0, 'cpdem_SolveEPAR mu');
assert(misfit >= 0 && misfit < 1e-3*mean(mu0), 'misfit of exact data');

% Invalid flags are errors, not silent misbehaviour
threw = false;
try
    cpdem_SolveEPAR(k0, mu0, p.sig0, p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 2, 2, 0, 2, 50);
catch
    threw = true;
end
assert(threw, 'invalid modFlag must throw');

end
