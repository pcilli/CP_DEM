function [p0, phi0] = cpdem_ForwardConduct(p1, sig1, p2, sig2, sig0, epar, opts)
% Cross-property DEM between two conductivity-like properties
%
%   p0 = cpdem_ForwardConduct(p1, sig1, p2, sig2, sig0, epar)
%   [p0, phi0] = cpdem_ForwardConduct(p1, sig1, p2, sig2, sig0, epar, opts)
%
% Implements the thermal-electrical cross-property DEM of Cilli and
% Chapman (2021, arXiv:2102.12434). Their electrical DEM (equations 1-2)
% and thermal DEM (equations 3-4), with the same geometric factor
% (m-bar for sigma, t-bar for kappa, identical functions), are combined by
% the chain rule into their equation 5,
%
%   d kappa/d sigma = (kappa2 - kappa) t-bar / ((sigma2 - sigma) m-bar),
%   kappa(sigma = sigma1) = kappa1,
%
% which this function integrates with p standing for kappa (or any other
% property obeying Laplace's equation: permittivity, magnetic permeability,
% diffusion constant). Their equation 6, electrical conductivity from
% thermal, is the reciprocal and is obtained by swapping the roles of the
% two properties in the argument list:
%
%   sig0 = cpdem_ForwardConduct(sig1, kappa1, sig2, kappa2, kappa0, epar).
%
% As for the elastic-electrical model, the inclusion volume fraction is a
% dummy variable and the relation depends only on the inclusion aspect
% ratio and the end-member properties. The exact special cases of the
% paper, spheres (EPAR = 1, equation 7), disks (EPAR = 0, equation 8) and
% needles (EPAR = Inf, equation 9), are evaluated for those EPAR values.
%
% Inputs
%   p1, sig1 : background (mineral) values of the two properties
%   p2, sig2 : inclusion values; either may be [axis plane] for an
%              inclusion material that is itself anisotropic (e.g. the
%              through-plane and in-plane conductivities of graphite
%              flakes), see cpdem_CondIncrement
%   sig0     : measured values of the "known" property, scalar or array,
%              between sig1 and sig2
%   epar     : aspect ratio (scalar, 0 and Inf allowed) or function handle
%              of inclusion fraction
%   opts     : optional cpdem_Options struct
%
% Outputs (shape of sig0): the modelled property p0 and implied porosity.
%
% See also cpdem_CalibrateConduct (the three-parameter inversion of the
% paper), cpdem_CondCondHS (the thermal-electrical Hashin-Shtrikman
% bounds) and cpdem_PublishedCuGraphite (the paper's parameters).
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
if ~isscalar(p1) || ~isscalar(sig1) || numel(p2) > 2 || numel(sig2) > 2
    error('CPDEM:badInput', ...
        'Host properties must be scalars; inclusion properties scalars or [axis plane] pairs.');
end

shape = size(sig0);
n = numel(sig0);
sig0 = sig0(:);
p0 = nan(n, 1);
phi0 = nan(n, 1);

lo = min([sig1 sig2(:).']); hi = max([sig1 sig2(:).']);
bad = ~(sig0 >= lo & sig0 <= hi) | any(sig0 == sig2(:).', 2);
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d values lie outside (sig1, sig2); NaN returned for them.', nnz(bad), n);
end

if any(~bad)
    if isnumeric(epar)
        % Constant aspect ratio: evaluate the geometry once and use the
        % inlined increments (several times faster in interpreted code).
        [~, ~, L] = cpdem_Shape(epar);
        [sa, sp] = axisPlane(sig2);
        [pa, pp] = axisPlane(p2);
        fun = @(u, y) derivFast(u, y, p1, sig1, sa, sp, pa, pp, L, (1 - L)/2);
    else
        fun = @(u, y) deriv(u, y, p1, p2, sig1, sig2, epar);
    end
    Y = cpdem_Integrate(fun, log(sig0(~bad)/sig1), [0; 0], opts);
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
g = sig/cpdem_CondIncrement(sig, sig2, alpha);     % = sig/((sig2 - sig) m-bar) when isotropic
dy = [ cpdem_CondIncrement(p, p2, alpha)/p*g;
       (1 - phi)*g ];
end

function dy = derivFast(u, y, p1, sig1, sa, sp, pa, pp, L, Lp)
% Same as deriv with the geometry precomputed: G = (1/3) sig sum_j
% (sig2_j - sig)/(sig + L_j (sig2_j - sig)), which is (sig2 - sig) m-bar
% for an isotropic inclusion.
sig = sig1*exp(u);
p   = p1*exp(y(1));
phi = y(2);
Gs = (sig/3)*((sa - sig)/(sig + L*(sa - sig)) + 2*(sp - sig)/(sig + Lp*(sp - sig)));
Gp = (p/3)*((pa - p)/(p + L*(pa - p)) + 2*(pp - p)/(p + Lp*(pp - p)));
g = sig/Gs;
dy = [ Gp/p*g;
       (1 - phi)*g ];
end

function [a, b] = axisPlane(v)
if isscalar(v)
    a = v; b = v;
else
    a = v(1); b = v(2);
end
end
