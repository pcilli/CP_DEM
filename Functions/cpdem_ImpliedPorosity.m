function phi = cpdem_ImpliedPorosity(sig0, sig1, sig2, epar, opts)
% Porosity implied by the electrical DEM for a measured conductivity
%
%   phi = cpdem_ImpliedPorosity(sig0, sig1, sig2, epar)
%   phi = cpdem_ImpliedPorosity(sig0, sig1, sig2, epar, opts)
%
% Inverts the electrical DEM (see cpdem_ElecDEM) by integrating
%
%   d phi/d sig = (1 - phi) / ((sig2 - sig) m-bar(sig, sig2, EPAR))
%
% from sig = sig1 (phi = 0) to each measured conductivity sig0. Porosity is
% a dummy variable in the cross-property DEM of Cilli and Chapman (2021)
% and does not enter its predictions; this function is a diagnostic
% extension of the toolbox that reports the porosity the DEM path "passes
% through" to reach a measurement. For a real rock it need not equal the
% true porosity (EPAR is an effective parameter), so treat it as a model
% quantity, e.g. for the 'implied' density option of cpdem_ModelVpVs.
%
% Inputs
%   sig0 : measured conductivity, scalar or array, between sig1 and sig2
%   sig1 : host (mineral) conductivity (> 0)
%   sig2 : inclusion (fluid) conductivity
%   epar : aspect ratio (scalar), per-sample array, or function handle of
%          porosity
%   opts : optional cpdem_Options struct
%
% Output (shape of sig0): implied porosity; NaN where sig0 is outside
% (sig1, sig2).
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 5 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
cpdem_CheckPositive('sig0', sig0, 'sig1', sig1, 'sig2', sig2);
if ~isscalar(sig1) || numel(sig2) > 2
    error('CPDEM:badInput', 'sig1 must be a scalar and sig2 a scalar or [sig_axis sig_plane].');
end

shape = size(sig0);
n = numel(sig0);
sig0 = sig0(:);
phi = nan(n, 1);

lo = min([sig1 sig2(:).']); hi = max([sig1 sig2(:).']);
bad = ~(sig0 >= lo & sig0 <= hi) | any(sig0 == sig2(:).', 2);
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d conductivities lie outside (sig1, sig2); NaN returned for them.', nnz(bad), n);
end

fun = @(u, p, e) sig1*exp(u)*(1 - p)/cpdem_CondIncrement(sig1*exp(u), sig2, cpdem_EvalEPAR(e, p));

if isnumeric(epar) && ~isscalar(epar)
    if numel(epar) ~= n
        error('CPDEM:badEPAR', 'EPAR must be a scalar, a function handle, or have one value per measurement.');
    end
    for ii = find(~bad).'
        Y = cpdem_Integrate(@(u, p) fun(u, p, epar(ii)), log(sig0(ii)/sig1), 0, opts);
        phi(ii) = Y(1);
    end
elseif any(~bad)
    Y = cpdem_Integrate(@(u, p) fun(u, p, epar), log(sig0(~bad)/sig1), 0, opts);
    phi(~bad) = Y(:, 1);
end

phi = reshape(phi, shape);

end
