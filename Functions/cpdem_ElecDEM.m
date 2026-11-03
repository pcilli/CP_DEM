function [sig, mb] = cpdem_ElecDEM(phi, sig1, sig2, epar, opts)
% Electrical DEM: effective conductivity as a function of porosity
%
%   sig = cpdem_ElecDEM(phi, sig1, sig2, epar)
%   [sig, mb] = cpdem_ElecDEM(phi, sig1, sig2, epar, opts)
%
% Integrates the differential effective medium equation for the electrical
% conductivity of a host (sig1) containing randomly oriented spheroidal
% inclusions (sig2, aspect ratio EPAR), the electrical "half" of the
% cross-property DEM of Cilli and Chapman (2021):
%
%   (1 - phi) d sig/d phi = (sig2 - sig) m-bar(sig, sig2, EPAR)
%
% with sig(0) = sig1. For spheres this is Bruggeman's asymmetric DEM with
% closed form (sig - sig2)/(sig1 - sig2) (sig1/sig)^(1/3) = 1 - phi; for
% insulating inclusions (sig2 = 0) it is Archie's law, sig = sig1 (1 -
% phi)^m-bar, with the aspect-ratio dependent cementation exponent of Sen
% et al. (1981) and Mendelson and Cohen (1982).
%
% Phase roles. In the cross-property DEM of Cilli and Chapman (2021) the
% mineral is the host and the brine-filled pores are the inclusions
% (sig1 = mineral, sig2 = fluid, phi = porosity). The electrical model of
% Cilli and Chapman (2020) and of Sen et al. (1981) instead embeds
% insulating grains in a conducting fluid host; to evaluate that form call
%   sig = cpdem_ElecDEM(1 - phi, sigFluid, sigGrain, grainEPAR)
% i.e. with the inclusion fraction 1 - phi and the grain aspect ratio. A
% porosity-dependent aspect ratio, as in the power law of Cilli and Chapman
% (2020), is passed as a function handle of the INCLUSION fraction (the
% state carried along the path): epar = @(f) law(1 - f) in the grain-
% inclusion form, epar = @(phi) law(phi) in the pore-inclusion form.
%
% Inputs
%   phi  : porosity, scalar or array, 0 <= phi < 1
%   sig1 : host conductivity (> 0)
%   sig2 : inclusion conductivity (> 0); or [sig2_axis sig2_plane] for an
%          anisotropic inclusion material (see cpdem_CondIncrement)
%   epar : aspect ratio (scalar) or function handle epar(phi)
%   opts : optional cpdem_Options struct
%
% Outputs (shape of phi)
%   sig  : effective conductivity
%   mb   : m-bar evaluated at (sig, sig2, EPAR(phi))
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Cilli, P. A., and Chapman, M. (2020). The power-law relation between
% inclusion aspect ratio and porosity: Implications for electrical and
% elastic modeling. Journal of Geophysical Research: Solid Earth, 125,
% e2019JB019187, DOI:10.1029/2019JB019187
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 5 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
cpdem_CheckPositive('sig1', sig1, 'sig2', sig2, 'phi>=0', phi);
if any(phi(:) >= 1)
    error('CPDEM:badInput', 'Porosity must be < 1.');
end
if ~isscalar(sig1) || numel(sig2) > 2
    error('CPDEM:badInput', 'sig1 must be a scalar and sig2 a scalar or [sig_axis sig_plane].');
end

shape = size(phi);
fun = @(p, y) cpdem_CondIncrement(sig1*exp(y), sig2, cpdem_EvalEPAR(epar, p)) ...
              /(sig1*exp(y)*(1 - p));
Y = cpdem_Integrate(fun, phi(:), 0, opts);
sig = reshape(sig1*exp(Y(:, 1)), shape);

if nargout > 1
    if ~isscalar(sig2)
        error('CPDEM:badInput', 'm-bar is only defined for isotropic inclusions.');
    end
    mb = zeros(shape);
    for ii = 1:numel(phi)
        mb(ii) = mbar(sig(ii), sig2, cpdem_EvalEPAR(epar, phi(ii)));
    end
end

end
