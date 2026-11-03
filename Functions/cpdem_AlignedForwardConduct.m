function out = cpdem_AlignedForwardConduct(p1, sig1, p2, sig2, sig0, epar, which, opts)
% Anisotropic cross-property DEM between two conductivity-like properties
%
%   out = cpdem_AlignedForwardConduct(p1, sig1, p2, sig2, sig0, epar, which)
%   out = cpdem_AlignedForwardConduct(..., which, opts)
%
% Cross-property modelling for a composite whose spheroidal inclusions are
% aligned, so that both properties (e.g. thermal conductivity p and
% electrical conductivity sig) are transversely isotropic with a common
% symmetry axis. The four DEM equations of cpdem_AlignedCondDEM (two
% components of each property) share the inclusion fraction, which the
% chain rule eliminates: with one measured component of sig as the
% independent variable, the other component of sig and both components of
% p follow without reference to the volume fraction, as in the isotropic
% model of Cilli and Chapman (2021, arXiv:2102.12434).
%
% Inputs
%   p1, sig1 : host values of the two properties; each a scalar or
%              [par perp] for a transversely isotropic host
%   p2, sig2 : inclusion values; each a scalar or [axis plane] for an
%              anisotropic inclusion material
%   sig0     : measured values of the chosen component of sig (array)
%   epar     : aspect ratio (scalar) or function handle of inclusion fraction
%   which    : 'par' if sig0 is measured along the symmetry axis,
%              'perp' if in the plane normal to it
%   opts     : optional cpdem_Options struct
%
% Output struct with fields (each shaped like sig0)
%   pPar, pPerp     : the modelled property along and across the axis
%   sigPar, sigPerp : both components of sig (one of them is sig0)
%   phi             : inclusion fraction along the path (diagnostic)
%
% This is an extension of the toolbox beyond the published isotropic
% model; see cpdem_AlignedIncrement for the derivation and its checks.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 8 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
cpdem_CheckPositive('p1', p1, 'sig1', sig1, 'p2', p2, 'sig2', sig2, 'sig0', sig0);
switch lower(which)
    case {'par', 'parallel', 'axis'}
        usePar = true;
    case {'perp', 'perpendicular', 'plane'}
        usePar = false;
    otherwise
        error('CPDEM:badInput', 'which must be ''par'' or ''perp''.');
end
[s1par, s1perp] = split2(sig1);
[p1par, p1perp] = split2(p1);
if usePar
    sRef = s1par;
else
    sRef = s1perp;
end

shape = size(sig0);
n = numel(sig0);
sig0 = sig0(:);

% Reachable range of the measured component: between the host value and
% the corresponding inclusion value(s).
s2 = sig2(:).';
lo = min([sRef s2]); hi = max([sRef s2]);
bad = ~(sig0 >= lo & sig0 <= hi) | (sig0 == sRef & false);
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d values lie outside the range spanned by the host and inclusion; NaN returned for them.', nnz(bad), n);
end

pPar = nan(n, 1); pPerp = nan(n, 1); sPar = nan(n, 1); sPerp = nan(n, 1); phi = nan(n, 1);
if any(~bad)
    % state y = [ln(sig_other/sig_other1); ln(p_par/p1par); ln(p_perp/p1perp); phi]
    fun = @(u, y) deriv(u, y, usePar, s1par, s1perp, p1par, p1perp, sig2, p2, epar);
    Y = cpdem_Integrate(fun, log(sig0(~bad)/sRef), [0; 0; 0; 0], opts);
    if usePar
        sPar(~bad)  = sig0(~bad);
        sPerp(~bad) = s1perp*exp(Y(:, 1));
    else
        sPerp(~bad) = sig0(~bad);
        sPar(~bad)  = s1par*exp(Y(:, 1));
    end
    pPar(~bad)  = p1par*exp(Y(:, 2));
    pPerp(~bad) = p1perp*exp(Y(:, 3));
    phi(~bad)   = Y(:, 4);
end

out = struct('pPar', reshape(pPar, shape), 'pPerp', reshape(pPerp, shape), ...
             'sigPar', reshape(sPar, shape), 'sigPerp', reshape(sPerp, shape), ...
             'phi', reshape(phi, shape));

end

function dy = deriv(u, y, usePar, s1par, s1perp, p1par, p1perp, sig2, p2, epar)
if usePar
    spar  = s1par*exp(u);
    sperp = s1perp*exp(y(1));
else
    sperp = s1perp*exp(u);
    spar  = s1par*exp(y(1));
end
ppar  = p1par*exp(y(2));
pperp = p1perp*exp(y(3));
phi   = y(4);
alpha = cpdem_EvalEPAR(epar, phi);

[Gspar, Gsperp] = cpdem_AlignedIncrement(spar, sperp, sig2, alpha);
[Gppar, Gpperp] = cpdem_AlignedIncrement(ppar, pperp, p2, alpha);

% d ln(sig_measured)/d phi = G_measured/(sig_measured (1 - phi)); divide
% every other rate by it (the (1 - phi) factors cancel).
if usePar
    rate = Gspar/spar;
    dOther = Gsperp/sperp;
else
    rate = Gsperp/sperp;
    dOther = Gspar/spar;
end
dy = [dOther; Gppar/ppar; Gpperp/pperp; (1 - phi)]/rate;
end

function [a, b] = split2(v)
if isscalar(v)
    a = v; b = v;
elseif numel(v) == 2
    a = v(1); b = v(2);
else
    error('CPDEM:badInput', 'Properties must be scalars or 2-vectors.');
end
end
