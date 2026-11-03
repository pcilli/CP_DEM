function [k, mu] = cpdem_ElasDEM(phi, k1, mu1, k2, mu2, epar, opts)
% Elastic DEM: effective bulk and shear moduli as a function of porosity
%
%   [k, mu] = cpdem_ElasDEM(phi, k1, mu1, k2, mu2, epar)
%   [k, mu] = cpdem_ElasDEM(phi, k1, mu1, k2, mu2, epar, opts)
%
% Integrates the differential effective medium equations of Berryman
% (1992) for a host (k1, mu1) containing randomly oriented spheroidal
% inclusions (k2, mu2, aspect ratio EPAR), the elastic "half" of the
% cross-property DEM of Cilli and Chapman (2021):
%
%   (1 - phi) dK /d phi = (K2  - K ) P(K, mu, K2, mu2, EPAR)
%   (1 - phi) dmu/d phi = (mu2 - mu) Q(K, mu, K2, mu2, EPAR)
%
% with K(0) = k1, mu(0) = mu1, and P, Q from Berryman (1980).
%
% Inputs
%   phi      : porosity, scalar or array, 0 <= phi < 1
%   k1, mu1  : host moduli (> 0)
%   k2, mu2  : inclusion moduli (k2 > 0, mu2 >= 0)
%   epar     : aspect ratio (scalar) or function handle epar(phi)
%   opts     : optional cpdem_Options struct
%
% Outputs (shape of phi): effective bulk and shear moduli.
%
% References:
% Berryman, J. G. (1992). Single-scattering approximations for coefficients
% in Biot's equations of poroelasticity. The Journal of the Acoustical
% Society of America, 91(2), 551-571, DOI:10.1121/1.402518
%
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 7 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
cpdem_CheckPositive('k1', k1, 'mu1', mu1, 'k2', k2, 'mu2>=0', mu2, 'phi>=0', phi);
if any(phi(:) >= 1)
    error('CPDEM:badInput', 'Porosity must be < 1.');
end
if ~all([isscalar(k1), isscalar(mu1), isscalar(k2), isscalar(mu2)])
    error('CPDEM:badInput', 'Moduli must be scalars.');
end

shape = size(phi);
Y = cpdem_Integrate(@(p, y) deriv(p, y, k1, k2, mu2, epar), phi(:), [1; mu1/k1], opts);
k  = reshape(Y(:, 1)*k1, shape);
mu = reshape(Y(:, 2)*k1, shape);

end

function dy = deriv(phi, y, k1, k2, mu2, epar)
K  = y(1)*k1;
mu = y(2)*k1;
[p, q] = cpdem_PQ(K, mu, k2, mu2, cpdem_EvalEPAR(epar, phi));
dy = [ (k2 - K)*p/(k1*(1 - phi));
       (mu2 - mu)*q/(k1*(1 - phi)) ];
end
