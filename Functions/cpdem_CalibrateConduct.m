function cal = cpdem_CalibrateConduct(pMeas, sigMeas, p1, sig1, p2Bounds, sig2Bounds, varargin)
% Calibrate the thermal-electrical cross-property DEM (aspect ratio and
% effective inclusion conductivities) to measurements, as in Section 3 of
% Cilli and Chapman (arXiv:2102.12434)
%
%   cal = cpdem_CalibrateConduct(pMeas, sigMeas, p1, sig1, p2Bounds, sig2Bounds)
%   cal = cpdem_CalibrateConduct(..., 'Name', value, ...)
%
% The paper's copper-graphite inclusions are transversely isotropic, so
% with randomly oriented inclusions treated as isotropic there are three
% unknowns: the inclusion aspect ratio and the effective thermal and
% electrical conductivities of the inclusion phase, the latter bounded by
% the values measured along the flake's long and short axes. All three are
% solved simultaneously by minimising the l2 norm of the misfit between
% measured and modelled thermal conductivity (equation 5, thermal from
% electrical), with the aspect ratio searched over [0, Inf] (disks to
% spheres to needles). The paper used the samples with inclusion fraction
% up to 0.4, as inclusion models assume a dilute dispersion.
%
% Inputs
%   pMeas, sigMeas : measured values of the modelled property p (kappa in
%                    the paper) and of the known property sig (electrical
%                    conductivity), arrays of equal size
%   p1, sig1       : host (matrix) values
%   p2Bounds       : [lo hi] admissible range of the effective inclusion
%                    value of p (e.g. the flake's two principal values); a
%                    scalar fixes it
%   sig2Bounds     : likewise for sig
%   options        : 'Norm' (2), 'Relative' (false), 'Weights', 'Display',
%                    'TolX' (1e-5 on the transformed parameters, i.e. about
%                    1e-5 relative in the aspect ratio), 'AlphaStarts'
%                    (initial aspect ratios of the multi-start Nelder-Mead
%                    search, default [0.1 1]; the exact limits 0 and Inf
%                    are always tried as well), 'Limits' (true: include
%                    those limits), 'AlphaBounds' ([1e-4 1e4], the
%                    interior of the aspect-ratio search in log space; the
%                    exact limits 0 and Inf complete the paper's [0, Inf])
%                    and the solver tolerances (the profile stage uses
%                    RelTol 1e-8 unless RelTol is given; the refinement
%                    uses the default 1e-10, which the narrow solution
%                    basin needs)
%
% Output struct cal with fields
%   epar, p2, sig2 : solved aspect ratio and effective inclusion values
%   misfit         : minimised misfit
%   pModel         : modelled p at sigMeas
%   rms, mape, r2  : fit statistics
%   starts         : misfit after refinement of each start (to judge
%                    uniqueness)
%   profile        : misfit minimised over p2 and sig2 at each aspect
%                    ratio of the grid (fields alpha, misfit): the
%                    landscape the paper's search over [0, Inf] explores
%
% With few samples the misfit has several local minima (needle-like
% aspect ratios with the inclusion values at their bounds fit almost as
% well as the true solution, whose basin is narrow in all three
% parameters), so the search is global: the misfit is profiled over a
% log-spaced grid of aspect ratios (option 'GridSize', default [9 3 3]:
% number of aspect ratios, and the starting grid of the 2-D search over
% p2 and sig2 at each), optimising the two inclusion values at each
% aspect ratio, and the profile is then minimised over log aspect ratio
% with a bounded one-dimensional search (fminbnd) in the bracket around
% its best grid point (and, if 'Refine' > 1, around the next best ones).
% With few samples the three parameters trade off strongly and the
% solution is a sharp, V-shaped minimum of the profile, where a direct
% three-parameter search (Nelder-Mead or Gauss-Newton) stalls in the
% valley; the profile search is robust to this. The two-dimensional
% sub-problem at fixed aspect ratio is solved by a Gauss-Newton
% (Levenberg-Marquardt) iteration for the l2 norm and by fminsearch for
% the l1 norm. Bounded parameters are handled with sine transforms and
% the aspect ratio is searched in log space, so nothing beyond core
% MATLAB or Octave is needed; the exact disk (0) and needle (Inf) limits
% are included and returned if they fit best. 'AlphaStarts' adds extra
% aspect ratios to the grid.
%
% Reference:
% Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal
% conductivity through cross-property inclusion modelling. arXiv:2102.12434
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

alphaStarts = [];
tryLimits = true;
alphaBounds = [1e-4 1e4];
gridSize = [9 3 3];
nRefine = 2;
args = varargin;
keep = true(size(args));
for ii = 1:2:numel(args)
    if ~ischar(args{ii}), continue; end
    switch lower(args{ii})
        case 'alphastarts', alphaStarts = args{ii+1}; keep([ii ii+1]) = false;
        case 'limits',      tryLimits = logical(args{ii+1}); keep([ii ii+1]) = false;
        case 'alphabounds', alphaBounds = args{ii+1}; keep([ii ii+1]) = false;
        case 'gridsize',    gridSize = args{ii+1}; keep([ii ii+1]) = false;
        case 'refine',      nRefine = args{ii+1}; keep([ii ii+1]) = false;
    end
end
opts = cpdem_Options('TolX', 1e-5, args{keep});
% The profile stage can use a looser solver tolerance; the refinement must
% not, because the solver's own noise (about RelTol in the modelled
% values) floors the misfit and blurs the narrow basin of the solution.
coarseOpts = opts;
if ~any(strcmpi(args(1:2:end), 'RelTol'))
    coarseOpts.RelTol = 1e-8; coarseOpts.AbsTol = 1e-11;
end

pMeas = pMeas(:); sigMeas = sigMeas(:);
if numel(pMeas) ~= numel(sigMeas)
    error('CPDEM:badInput', 'pMeas and sigMeas must have the same number of elements.');
end
cpdem_CheckPositive('pMeas', pMeas, 'sigMeas', sigMeas, 'p1', p1, 'sig1', sig1, ...
    'p2Bounds', p2Bounds, 'sig2Bounds', sig2Bounds);
w = opts.Weights;
if isempty(w), w = ones(size(pMeas)); end
w = w(:);

pB = sort(p2Bounds(:).');   if isscalar(pB), pB = [pB pB]; end
sB = sort(sig2Bounds(:).'); if isscalar(sB), sB = [sB sB]; end

ctx = struct('pMeas', pMeas, 'sigMeas', sigMeas, 'p1', p1, 'sig1', sig1, ...
             'pB', pB, 'sB', sB, 'w', w, 'opts', opts, 'la', log(alphaBounds(1)), ...
             'lb', log(alphaBounds(2)));
ctxCoarse = ctx;
ctxCoarse.opts = coarseOpts;

fopts = optimset('TolX', opts.TolX, 'TolFun', opts.TolX*1e-3, 'MaxFunEvals', 1500, ...
                 'MaxIter', 1500, 'Display', opts.Display);

% --- Profile of the misfit over aspect ratio: at each alpha of a
% log-spaced grid (plus the exact limits), optimise p2 and sig2 (2-D)
alphas = sort([logspace(log10(alphaBounds(1)), log10(alphaBounds(2)), gridSize(1)) alphaStarts(:).']);
if tryLimits
    alphas = [0 alphas Inf];
end
% Starting grid of the 2-D search, uniform in p2 and sig2 (sine transform),
% kept just inside the bounds where the transform's derivative vanishes
xp = asin(linspace(-0.95, 0.95, gridSize(2))); if pB(1) == pB(2), xp = 0; end
xs = asin(linspace(-0.95, 0.95, gridSize(3))); if sB(1) == sB(2), xs = 0; end
[XP2, XS2] = ndgrid(xp, xs);
profile = nan(size(alphas));
x2 = cell(size(alphas));
for ia = 1:numel(alphas)
    a = alphas(ia);
    gc = nan(size(XP2));
    for ii = 1:numel(XP2)
        gc(ii) = cost([NaN; XP2(ii); XS2(ii)], ctxCoarse, a);
    end
    [cg, kk] = min(gc(:));
    x0 = [XP2(kk); XS2(kk)];
    % Warm start from the previous aspect ratio's solution when it is better
    if ia > 1 && isfinite(profile(ia-1)) && cost([NaN; x2{ia-1}], ctxCoarse, a) < cg
        x0 = x2{ia-1};
    end
    [x2{ia}, profile(ia)] = solve2D(x0, ctxCoarse, a, fopts);
end

% --- Minimise the profile over log alpha in the bracket(s) around its
% best interior grid point(s), with the 2-D sub-problem solved at the
% user's (fine) solver tolerance and warm-started from the grid solution
interior = find(isfinite(alphas) & alphas > 0);
[~, order] = sort(profile(interior));
best = struct('cost', Inf, 'x', [], 'epar', NaN);
starts = nan(0, 1);
bopts = optimset('TolX', opts.TolX, 'Display', opts.Display);
for ii = 1:min(nRefine, numel(order))
    ia = interior(order(ii));
    lo = max(alphaBounds(1), alphas(max(ia - 1, 1)));
    hi = min(alphaBounds(2), alphas(min(ia + 1, numel(alphas))));
    if ~isfinite(hi), hi = alphaBounds(2); end
    if lo == 0, lo = alphaBounds(1); end
    warm = x2{ia};
    prof = @(la) profileCost(exp(la), warm, ctx, fopts);
    [laBest, c] = fminbnd(prof, log(lo), log(hi), bopts);
    [c, xBest] = profileCost(exp(laBest), warm, ctx, fopts);
    starts(end+1) = c; %#ok<AGROW>
    if c < best.cost
        best.cost = c; best.x = [NaN; xBest]; best.epar = exp(laBest);
    end
end

% --- The exact disk and needle limits of the paper (alpha = 0, Inf):
% polish their 2-D solutions at the fine tolerance and keep them if best
if tryLimits
    for ia = find(~isfinite(alphas) | alphas == 0)
        a = alphas(ia);
        [x, c] = solve2D(x2{ia}, ctx, a, fopts);
        starts(end+1) = c; %#ok<AGROW>
        if c < best.cost
            best.cost = c; best.x = [NaN; x]; best.epar = a;
        end
    end
end

[p2, sig2] = unpack(best.x, ctx);
cal.epar   = best.epar;
cal.p2     = p2;
cal.sig2   = sig2;
cal.misfit = best.cost;
cal.pModel = cpdem_ForwardConduct(p1, sig1, p2, sig2, sigMeas, cal.epar, opts);
r = pMeas - cal.pModel;
cal.rms  = sqrt(mean(r.^2));
cal.mape = 100*mean(abs(r./pMeas));
sst = sum((pMeas - mean(pMeas)).^2);
if sst > 0, cal.r2 = 1 - sum(r.^2)/sst; else, cal.r2 = NaN; end
cal.starts = starts;
cal.profile = struct('alpha', alphas, 'misfit', profile);
cal.bounds = struct('p2', pB, 'sig2', sB, 'alpha', alphaBounds);

end

function a = alphaOfX(x, ctx)
% Aspect ratio from the unbounded search variable: log-uniform over AlphaBounds
a = exp(ctx.la + (ctx.lb - ctx.la)*(1 + sin(x))/2);
end

function x = xOfAlpha(a, ctx)
x = asin(min(max(2*(log(a) - ctx.la)/(ctx.lb - ctx.la) - 1, -1), 1));
end

function [p2, sig2] = unpack(x, ctx)
p2   = ctx.pB(1) + (ctx.pB(2) - ctx.pB(1))*(1 + sin(x(2)))/2;
sig2 = ctx.sB(1) + (ctx.sB(2) - ctx.sB(1))*(1 + sin(x(3)))/2;
end

function r = resid(x, ctx, alphaFixed)
% Weighted (optionally relative) residual vector; NaN where unreachable.
% x = [x_alpha; x_p2; x_sig2], or [x_p2; x_sig2] when alphaFixed is given
% (x_alpha may then be NaN in a 3-vector).
if isempty(alphaFixed)
    alpha = alphaOfX(x(1), ctx);
    xr = x(2:3);
else
    alpha = alphaFixed;
    xr = x(end-1:end);
end
[p2, sig2] = unpack([NaN; xr(:)], ctx);
ws = warning('off', 'CPDEM:outOfRange');
try
    model = cpdem_ForwardConduct(ctx.p1, ctx.sig1, p2, sig2, ctx.sigMeas, alpha, ctx.opts);
catch
    model = nan(size(ctx.sigMeas));
end
warning(ws);
r = ctx.pMeas - model;
if ctx.opts.Relative
    r = r./ctx.pMeas;
end
r = ctx.w.*r;
end

function c = cost(x, ctx, alphaFixed)
if nargin < 3
    alphaFixed = [];
end
r = resid(x, ctx, alphaFixed);
if any(isnan(r))
    c = 1e300;     % measurement unreachable with these trial end members
elseif ctx.opts.Norm == 1
    c = sum(abs(r));
else
    c = sqrt(sum(r.^2));
end
end

function [c, x] = profileCost(alpha, warm, ctx, fopts)
% Misfit minimised over p2 and sig2 at a fixed aspect ratio
[x, c] = solve2D(warm, ctx, alpha, fopts);
end

function [x, c] = solve2D(x0, ctx, alphaFixed, fopts)
% Two-dimensional sub-problem at fixed aspect ratio: Levenberg-Marquardt
% for the l2 norm, fminsearch otherwise. x0 = [x_p2; x_sig2].
x0 = x0(:);
if ctx.opts.Norm ~= 2
    [x, c] = fminsearch(@(x) cost(x, ctx, alphaFixed), x0, fopts);
    return
end
x = x0;
r = resid(x, ctx, alphaFixed);
c = sqrt(sum(r.^2));
if any(isnan(r))
    c = 1e300;
    return
end
lambda = 1e-3;
h = 1e-4;
n = numel(x);
for it = 1:60
    % Central-difference Jacobian
    J = zeros(numel(r), n);
    for jj = 1:n
        dx = zeros(n, 1); dx(jj) = h;
        rp = resid(x + dx, ctx, alphaFixed);
        rm = resid(x - dx, ctx, alphaFixed);
        J(:, jj) = (rp - rm)/(2*h);
    end
    if any(isnan(J(:)))
        break
    end
    % A parameter sitting exactly on a bound (sin = +-1) has a zero
    % column; nudge it inside and recompute once.
    dead = find(all(J == 0, 1) & abs(abs(sin(x.')) - 1) < 1e-9);
    if ~isempty(dead)
        x(dead) = sign(x(dead)).*(pi/2 - 1e-2);
        r = resid(x, ctx, alphaFixed);
        c = sqrt(sum(r.^2));
        continue
    end
    % Fixed parameters (equal bounds) have zero columns: solve on the
    % active ones only.
    act = any(J ~= 0, 1);
    if ~any(act)
        break
    end
    A = J(:, act)'*J(:, act);
    g = J(:, act)'*r;
    scale = max(trace(A)/nnz(act), 1e-300);
    improved = false;
    while lambda < 1e10
        step = zeros(n, 1);
        step(act) = -(A + lambda*diag(diag(A)) + 1e-12*scale*eye(nnz(act)))\g;
        rn = resid(x + step, ctx, alphaFixed);
        cn = sqrt(sum(rn.^2));
        if ~any(isnan(rn)) && cn < c
            x = x + step; r = rn; c = cn;
            lambda = max(lambda/3, 1e-12);
            improved = true;
            break
        end
        lambda = lambda*10;
    end
    if ~improved || norm(step) < ctx.opts.TolX || c < 1e-10*max(1, sqrt(sum((ctx.w.*ctx.pMeas).^2)))
        break
    end
end
end
