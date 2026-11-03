function [p0, phi0] = cpdem_ForwardConduct(p1, sig1, p2, sig2, sig0, epar, opts)
% Cross-property DEM between two conductivity-like properties
%
%   p0 = cpdem_ForwardConduct(p1, sig1, p2, sig2, sig0, epar)
%   [p0, phi0] = cpdem_ForwardConduct(p1, sig1, p2, sig2, sig0, epar, opts)
%
% Models one transport property that obeys Laplace's equation (thermal
% conductivity, dielectric permittivity, magnetic permeability, ...) from
% another (here written as electrical conductivity) using the chain rule
% on the two DEM equations, which share the geometric factor m-bar
% (Cilli & Chapman 2021, arXiv:2102.12434):
%
%   d p/d sig = (p2 - p) m-bar(p, p2, EPAR) / ((sig2 - sig) m-bar(sig, sig2, EPAR))
%
% with p(sig1) = p1. As for the elastic-electrical model, porosity is a
% dummy variable and the relation depends only on the pore aspect ratio and
% the end-member properties.
%
% Inputs
%   p1, sig1 : background (mineral) values of the two properties
%   p2, sig2 : inclusion (fluid) values
%   sig0     : measured values of the "known" property, scalar or array,
%              between sig1 and sig2
%   epar     : aspect ratio (scalar) or function handle of porosity
%   opts     : optional cpdem_Options struct
%
% Outputs (shape of sig0): the modelled property p0 and implied porosity.
%
% Reference:
% Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal
% conductivity through cross-property inclusion modelling.
% arXiv:2102.12434.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 7 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
cpdem_CheckPositive('p1', p1, 'sig1', sig1, 'p2', p2, 'sig2', sig2, 'sig0', sig0);
if ~all([isscalar(p1), isscalar(sig1), isscalar(p2), isscalar(sig2)])
    error('CPDEM:badInput', 'End-member properties must be scalars.');
end

shape = size(sig0);
n = numel(sig0);
sig0 = sig0(:);
p0 = nan(n, 1);
phi0 = nan(n, 1);

bad = ~cpdem_InRange(sig0, repmat(sig1, n, 1), repmat(sig2, n, 1));
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d values lie outside (sig1, sig2); NaN returned for them.', nnz(bad), n);
end

if any(~bad)
    Y = cpdem_Integrate(@(u, y) deriv(u, y, p1, p2, sig1, sig2, epar), ...
                        log(sig0(~bad)/sig1), [0; 0], opts);
    p0(~bad) = p1*exp(Y(:, 1));
    phi0(~bad) = Y(:, 2);
end

p0 = reshape(p0, shape);
phi0 = reshape(phi0, shape);

end

function dy = deriv(u, y, p1, p2, sig1, sig2, epar)
% u = ln(sig/sig1); y = [ln(p/p1); phi]
sig = sig1*exp(u);
p   = p1*exp(y(1));
phi = y(2);
alpha = cpdem_EvalEPAR(epar, phi);
g = sig/((sig2 - sig)*mbar(sig, sig2, alpha));
dy = [ (p2 - p)*mbar(p, p2, alpha)/p*g;
       (1 - phi)*g ];
end
