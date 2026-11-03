function cal = cpdem_CalibrateEPAR(k0, mu0, sig0, k1,mu1,sig1,k2,mu2,sig2, varargin)
% Calibrate the bulk and shear modulus EPARs to joint elastic-electrical data
%
%   cal = cpdem_CalibrateEPAR(k0, mu0, sig0, k1,mu1,sig1,k2,mu2,sig2)
%   cal = cpdem_CalibrateEPAR(k0, mu0, sig0, k1,mu1,sig1,k2,mu2,sig2, opts)
%   cal = cpdem_CalibrateEPAR(..., 'Name', value, ...)
%
% Finds the equivalent pore aspect ratios alpha_K and alpha_mu that
% minimise the misfit between the measured bulk and shear moduli and those
% forward modelled from the measured conductivities with cpdem_Forward
% (Cilli & Chapman 2021). Either k0 or mu0 may be empty to skip that
% modulus. The search is a bounded scalar minimisation (fminbnd) over
% opts.Bounds, in log(EPAR) by default, so the whole data set costs a
% single ODE integration per trial EPAR when the parameters are shared.
%
% Inputs
%   k0, mu0, sig0 : measured bulk modulus, shear modulus (Pa) and
%                   conductivity (S/m); arrays of equal size (k0 or mu0
%                   may be [])
%   k1,mu1,sig1   : background (mineral) properties (scalars or per-sample)
%   k2,mu2,sig2   : inclusion (fluid) properties (scalars or per-sample)
%   opts          : cpdem_Options struct and/or name/value pairs, e.g.
%                   cpdem_CalibrateEPAR(..., 'Bounds', [2 50], 'Norm', 1)
%                   Fields used: Bounds, LogSearch, Norm, Relative,
%                   Weights, TolX, Display, PerSample (see below) and the
%                   solver tolerances.
%
% Output struct cal with fields
%   eparK, eparMu     : calibrated EPARs (NaN if that modulus was skipped)
%   misfitK, misfitMu : minimised misfit (in the chosen norm)
%   kModel, muModel   : modelled moduli at the calibrated EPARs
%   resK, resMu       : residuals (measured - modelled)
%   rmsK, rmsMu       : root-mean-square residuals (Pa)
%   mapeK, mapeMu     : mean absolute percentage error (%)
%   r2K, r2Mu         : coefficient of determination
%   bounds, norm, relative : settings used
%
% Per-sample EPARs: pass the option PerSample = true to calibrate an EPAR
% for every sample individually (eparK and eparMu then have the size of
% sig0). Cilli and Chapman (2021) report alpha_K = 16.4 and alpha_mu =
% 12.8 for the brine-saturated sandstones of Han et al. (2011) at 60 MPa,
% with the model working better for prolate (EPAR > 1) than oblate pores;
% the misfit may have several local minima if oblate EPARs are included
% in the bounds.
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Han, T., Best, A.I., Sothcott, J. and MacGregor, L.M. (2011), Joint
% elastic-electrical properties of reservoir sandstones and their
% relationships with petrophysical parameters. Geophysical Prospecting,
% 59, 518-535, DOI:10.1111/j.1365-2478.2010.00940.x
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

opts = cpdem_Options(varargin{:});
perSample = logical(opts.PerSample);

sig0 = sig0(:);
n = numel(sig0);
if isempty(k0) && isempty(mu0)
    error('CPDEM:badInput', 'At least one of k0 and mu0 must be given.');
end
if ~isempty(k0) && numel(k0) ~= n
    error('CPDEM:badInput', 'k0 must have one value per conductivity measurement.');
end
if ~isempty(mu0) && numel(mu0) ~= n
    error('CPDEM:badInput', 'mu0 must have one value per conductivity measurement.');
end
k0 = k0(:); mu0 = mu0(:);

w = opts.Weights;
if isempty(w)
    w = ones(n, 1);
elseif numel(w) ~= n
    error('CPDEM:badInput', 'Weights must have one value per measurement.');
end
w = w(:);

lb = opts.Bounds(1); ub = opts.Bounds(2);
if ~(lb > 0 && ub > lb)
    error('CPDEM:badInput', 'Bounds must be [lb ub] with 0 < lb < ub.');
end

cal = struct('eparK', NaN, 'eparMu', NaN, 'misfitK', NaN, 'misfitMu', NaN, ...
             'kModel', [], 'muModel', [], 'resK', [], 'resMu', [], ...
             'rmsK', NaN, 'rmsMu', NaN, 'mapeK', NaN, 'mapeMu', NaN, ...
             'r2K', NaN, 'r2Mu', NaN, 'bounds', opts.Bounds, 'norm', opts.Norm, ...
             'relative', opts.Relative);

ctx = struct('opts', opts, 'perSample', perSample, 'n', n, 'sig0', sig0, 'w', w, ...
             'lb', lb, 'ub', ub, 'k1', k1, 'mu1', mu1, 'sig1', sig1, ...
             'k2', k2, 'mu2', mu2, 'sig2', sig2);

% Bulk modulus EPAR
if ~isempty(k0)
    [e, m, model] = solveEPAR(ctx, k0, true);
    cal.eparK = e; cal.misfitK = m; cal.kModel = model;
    cal.resK = k0 - model;
    [cal.rmsK, cal.mapeK, cal.r2K] = fitStats(k0, model);
end
% Shear modulus EPAR
if ~isempty(mu0)
    [e, m, model] = solveEPAR(ctx, mu0, false);
    cal.eparMu = e; cal.misfitMu = m; cal.muModel = model;
    cal.resMu = mu0 - model;
    [cal.rmsMu, cal.mapeMu, cal.r2Mu] = fitStats(mu0, model);
end

end

function [e, m, model] = solveEPAR(ctx, meas, useK)
fopts = optimset('TolX', ctx.opts.TolX, 'Display', ctx.opts.Display);
if ctx.perSample
    n = ctx.n;
    e = nan(n, 1); m = nan(n, 1); model = nan(n, 1);
    for ii = 1:n
        prm = pickParams(ctx, ii);
        cost = @(x) misfit(ctx, prm, meas(ii), ctx.sig0(ii), x, useK, ctx.w(ii));
        [e(ii), m(ii)] = minimise(ctx, cost, fopts);
        model(ii) = fwd(ctx, prm, ctx.sig0(ii), e(ii), useK);
    end
else
    prm = pickParams(ctx, []);
    cost = @(x) misfit(ctx, prm, meas, ctx.sig0, x, useK, ctx.w);
    [e, m] = minimise(ctx, cost, fopts);
    model = fwd(ctx, prm, ctx.sig0, e, useK);
end
end

function [e, m] = minimise(ctx, cost, fopts)
if ctx.opts.LogSearch
    [x, m] = fminbnd(@(lx) cost(exp(lx)), log(ctx.lb), log(ctx.ub), fopts);
    e = exp(x);
else
    [e, m] = fminbnd(cost, ctx.lb, ctx.ub, fopts);
end
end

function c = misfit(ctx, prm, meas, s, e, useK, wt)
model = fwd(ctx, prm, s, e, useK);
r = meas - model;
if ctx.opts.Relative
    r = r./meas;
end
r = wt.*r;
if ctx.opts.Norm == 1
    c = sum(abs(r));
elseif ctx.opts.Norm == 2
    c = sqrt(sum(r.^2));
else
    error('CPDEM:badOptions', 'Norm must be 1 or 2.');
end
end

function model = fwd(ctx, prm, s, e, useK)
[kk, mm] = cpdem_Forward(prm.k1, prm.mu1, prm.sig1, prm.k2, prm.mu2, prm.sig2, s, e, ctx.opts);
if useK
    model = kk;
else
    model = mm;
end
end

function prm = pickParams(ctx, ii)
% Parameters for sample ii (per-sample calibration) or all (ii empty)
names = {'k1', 'mu1', 'sig1', 'k2', 'mu2', 'sig2'};
prm = struct();
for jj = 1:numel(names)
    v = ctx.(names{jj});
    if ~isempty(ii) && ~isscalar(v)
        v = v(ii);
    end
    prm.(names{jj}) = v;
end
end

function [rms, mape, r2] = fitStats(meas, model)
r = meas - model;
rms = sqrt(mean(r.^2));
mape = 100*mean(abs(r./meas));
sst = sum((meas - mean(meas)).^2);
if sst > 0
    r2 = 1 - sum(r.^2)/sst;
else
    r2 = NaN;
end
end
