function p = cpdem_ExampleParameters()
% Constituent properties shared by the example scripts (SI units)
%
%   p = cpdem_ExampleParameters()
%
% Representative values for a quartz sandstone saturated with brine. They
% are NOT the calibration values of Cilli and Chapman (2021); replace them
% with the properties appropriate to your own rocks and fluids. The
% toolbox has only been used with moduli in Pa and conductivity in S/m.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

% Background (mineral) phase: quartz
p.k1   = 37e9;    % bulk modulus, Pa
p.mu1  = 44e9;    % shear modulus, Pa
p.sig1 = 1e-4;    % conductivity, S/m. Must be > 0: a DEM with a perfectly
                  % insulating host never conducts, however much fluid is
                  % added, so a small finite mineral conductivity is needed.
p.rho1 = 2650;    % density, kg/m^3

% Inclusion (pore fluid) phase: brine
p.k2   = 2.25e9;  % bulk modulus, Pa
p.mu2  = 0;       % shear modulus, Pa (inviscid)
p.sig2 = 5;       % conductivity, S/m
p.rho2 = 1050;    % density, kg/m^3

% Published EPARs (Cilli & Chapman 2021, all sandstones of Han et al. 2011)
e = cpdem_PublishedEPAR('all');
p.eparK  = e.eparK;   % 16.4
p.eparMu = e.eparMu;  % 12.8

end
