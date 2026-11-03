function d = cpdem_SyntheticData(p, noiseLevel)
% Synthetic joint elastic-electrical data set for the example scripts
%
%   d = cpdem_SyntheticData(p)             2 % noise
%   d = cpdem_SyntheticData(p, noiseLevel) fractional noise level
%
% Generates a SYNTHETIC set of "measurements" from the cross-property DEM
% with the published EPARs in p (bulk modulus with alpha_K, shear modulus
% with alpha_mu), then adds deterministic pseudo-random noise so that the
% examples give identical numbers on MATLAB and Octave. It stands in for a
% joint laboratory data set such as Han et al. (2011), which is not
% redistributed with the toolbox.
%
% Fields of d: sig0, res0 (resistivity), k0, mu0, phiTrue (porosity the
% model passed through), rho (density from phiTrue), and the noise-free
% kTrue, muTrue.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 2
    noiseLevel = 0.02;
end

n = 24;
% Conductivities spread evenly in log resistivity between 2 and 60 ohm m
d.res0 = logspace(log10(2), log10(60), n);
d.sig0 = 1./d.res0;

[d.kTrue, ~, phiK]   = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, d.sig0, p.eparK);
[~, d.muTrue, phiMu] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, d.sig0, p.eparMu);
d.phiTrue = (phiK + phiMu)/2;
d.rho = (1 - d.phiTrue)*p.rho1 + d.phiTrue*p.rho2;

% Deterministic pseudo-noise (identical on every platform)
ii = 1:n;
nk  = cos(2.399*ii + 0.3).*sin(0.731*ii);
nmu = sin(1.871*ii + 1.1).*cos(0.517*ii);
d.k0  = d.kTrue .*(1 + noiseLevel*nk /max(abs(nk)));
d.mu0 = d.muTrue.*(1 + noiseLevel*nmu/max(abs(nmu)));
d.noiseLevel = noiseLevel;

end
