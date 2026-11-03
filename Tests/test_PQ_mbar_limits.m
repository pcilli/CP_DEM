function test_PQ_mbar_limits()
% P, Q and m-bar reduce to the known closed forms for spheres and to
% Archie's cementation exponent for insulating inclusions, and reproduce
% the v1.0 values away from the sphere.

p = cpdem_TestParams();
k0 = 30e9; mu0 = 35e9;

% Sphere: P = (K + 4/3 mu)/(K2 + 4/3 mu), Q = (mu + zeta)/(mu2 + zeta)
[P1, Q1] = cpdem_PQ(k0, mu0, p.k2, p.mu2, 1);
zeta = mu0/6*(9*k0 + 8*mu0)/(k0 + 2*mu0);
cpdem_AssertClose(P1, (k0 + 4/3*mu0)/(p.k2 + 4/3*mu0), 1e-12, 0, 'sphere P');
cpdem_AssertClose(Q1, (mu0 + zeta)/(p.mu2 + zeta), 1e-12, 0, 'sphere Q');

% Sphere m-bar is Bruggeman's 3 sig/(2 sig + sig2)
cpdem_AssertClose(mbar(0.1, 5, 1), 3*0.1/(2*0.1 + 5), 1e-12, 0, 'sphere m-bar');

% Insulating inclusions: m-bar = (1/3) sum 1/(1 - L_j), 3/2 for spheres
cpdem_AssertClose(mbar(0.3, 0, 1), 1.5, 1e-12, 0, 'Archie m for spheres');
[~, ~, L] = cpdem_Shape(16.4);
cpdem_AssertClose(mbar(0.3, 0, 16.4), (1/3)*(1/(1 - L) + 4/(1 + L)), 1e-12, 0, 'Archie m for needles');
% ... and it is independent of conductivity in that limit
cpdem_AssertClose(mbar([0.01 0.1 1], 0, 0.2), mbar(0.5, 0, 0.2)*[1 1 1], 1e-12, 0, 'Archie m independent of sig');

% R = m-bar/(3 sig)
cpdem_AssertClose(R(0.1, 5, 16.4), mbar(0.1, 5, 16.4)/(3*0.1), 1e-14, 0, 'R vs m-bar');

% v1.0 regression values (prolate and oblate, away from the sphere)
cpdem_AssertClose([P(k0,mu0,p.k2,p.mu2,16.4) Q(k0,mu0,p.k2,p.mu2,16.4) mbar(0.1,5,16.4)], ...
    [1.7378004867 2.50792849099 0.255250509316], 1e-10, 0, 'v1.0 values at 16.4');
cpdem_AssertClose([P(k0,mu0,p.k2,p.mu2,0.1) Q(k0,mu0,p.k2,p.mu2,0.1) mbar(0.1,5,0.1)], ...
    [4.01804751418 4.82942412355 0.158881189164], 1e-10, 0, 'v1.0 values at 0.1');

% Vectorised evaluation along a path agrees with scalar calls
ks = [37 30 20]*1e9; mus = [44 35 20]*1e9;
[Pv, Qv] = cpdem_PQ(ks, mus, p.k2, p.mu2, 12.8);
for ii = 1:3
    [Ps, Qs] = cpdem_PQ(ks(ii), mus(ii), p.k2, p.mu2, 12.8);
    cpdem_AssertClose([Pv(ii) Qv(ii)], [Ps Qs], 1e-14, 0, 'vectorised P,Q');
end

end
