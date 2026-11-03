function p = cpdem_TestParams()
% Shared constituent properties for the test suite (quartz / brine, SI)
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

p.k1   = 37e9;    % quartz bulk modulus, Pa
p.mu1  = 44e9;    % quartz shear modulus, Pa
p.sig1 = 1e-4;    % small but non-zero mineral conductivity, S/m
p.k2   = 2.25e9;  % brine bulk modulus, Pa
p.mu2  = 0;       % inviscid fluid
p.sig2 = 5;       % brine conductivity, S/m
p.sig0 = [0.02 0.05 0.1 0.2 0.4];  % effective conductivities, S/m
p.eparK  = 16.4;  % Cilli & Chapman (2021), all sandstones of Han et al. (2011)
p.eparMu = 12.8;

end
