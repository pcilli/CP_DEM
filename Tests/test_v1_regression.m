function test_v1_regression()
% The v2.0 forward model agrees with v1.0 to the accuracy v1.0 had.
%
% The reference values were produced with the v1.0 code (the bogus fourth
% argument to ode45 removed, as MATLAB silently ignored it, so that the
% solver ran with its default RelTol = 1e-3) in GNU Octave. v1.0 therefore
% carries solver noise of up to a few percent (it is not even monotonic in
% conductivity, see test_forward_composition); v2.0 is converged to ~1e-9,
% so the comparison is made at the 5 % level. The independent Python port
% (python/test_cpdem.py) checks v2.0 to much tighter tolerance.

p = cpdem_TestParams();

v1 = struct();
v1.k_16p4  = [2.867370895e+10 2.685919369e+10 2.447090003e+10 2.174589621e+10 1.810405698e+10];
v1.mu_16p4 = [2.982343247e+10 2.696852459e+10 2.356126868e+10 1.9980613e+10 1.541706247e+10];
v1.k_12p8  = [2.593278515e+10 2.357875012e+10 2.129167722e+10 1.865675263e+10 1.522033227e+10];
v1.mu_12p8 = [2.557290781e+10 2.257919733e+10 1.98239771e+10 1.657102609e+10 1.222967791e+10];
v1.k_0p1   = [7669514667 6456594845 5718044900 5066364117 4460774275];
v1.mu_0p1  = [4377053648 3043265735 2276941905 1641793564 1094711824];

[k, mu] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 16.4);
cpdem_AssertClose(k, v1.k_16p4, 5e-2, 0, 'K, EPAR 16.4');
cpdem_AssertClose(mu, v1.mu_16p4, 5e-2, 0, 'mu, EPAR 16.4');
[k, mu] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 12.8);
cpdem_AssertClose(k, v1.k_12p8, 5e-2, 0, 'K, EPAR 12.8');
cpdem_AssertClose(mu, v1.mu_12p8, 5e-2, 0, 'mu, EPAR 12.8');
[k, mu] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 0.1);
cpdem_AssertClose(k, v1.k_0p1, 5e-2, 0, 'K, EPAR 0.1');
cpdem_AssertClose(mu, v1.mu_0p1, 5e-2, 0, 'mu, EPAR 0.1');

% v1.0 Hashin-Shtrikman values (unchanged algorithm)
[kUp, kLow, muUp] = cpdem_ElecElasHS(p.sig0, p.sig1, p.sig2, p.k1, p.k2, p.mu1);
cpdem_AssertClose(kUp,  [3.66759e+10 3.61939e+10 3.54076e+10 3.38955e+10 3.10932e+10], 1e-5, 0, 'HS kUp');
cpdem_AssertClose(kLow, [2.2817e+09 2.26257e+09 2.25622e+09 2.25304e+09 2.25146e+09], 1e-5, 0, 'HS kLow');
cpdem_AssertClose(muUp, [4.34544e+10 4.2649e+10 4.13509e+10 3.89083e+10 3.45593e+10], 1e-5, 0, 'HS muUp');

end
