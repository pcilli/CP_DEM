function Y = cpdem_Integrate(odefun, targets, y0, opts)
% Integrate an ODE system from x = 0 to a set of target values of x
%
%   Y = cpdem_Integrate(odefun, targets, y0, opts)
%
% Solves dy/dx = odefun(x, y), y(0) = y0, and returns the state at each
% element of TARGETS as rows of Y (numel(targets) x numel(y0)). Targets may
% be unsorted, repeated, zero, or of either sign: the positive and negative
% targets are each handled with a single call to the ODE solver by passing
% the sorted unique targets as the solver's output grid. This is what
% allows the toolbox to model a whole data set with one integration.
%
% opts is a cpdem_Options struct (RelTol, AbsTol, MaxStep, Solver).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

targets = targets(:);
y0 = y0(:);
n = numel(targets);
Y = nan(n, numel(y0));

isZero = (targets == 0);
if any(isZero)
    Y(isZero, :) = repmat(y0.', nnz(isZero), 1);
end
if any(~isfinite(targets))
    error('CPDEM:badInput', 'Integration targets must be finite.');
end

odeopts = odeset('RelTol', opts.RelTol, 'AbsTol', opts.AbsTol);
if isfield(opts, 'MaxStep') && ~isempty(opts.MaxStep)
    odeopts = odeset(odeopts, 'MaxStep', opts.MaxStep);
end

for sgn = [1 -1]
    idx = find(sign(targets) == sgn);
    if isempty(idx)
        continue
    end
    xs = unique(targets(idx));         % ascending
    if sgn < 0
        xs = flipud(xs);               % increasing distance from 0
    end
    if numel(xs) == 1
        tspan = [0; xs/2; xs];
        keep  = 3;
    else
        tspan = [0; xs];
        keep  = 2:numel(tspan);
    end

    [X, YY] = feval(opts.Solver, odefun, tspan, y0, odeopts);

    if numel(X) ~= numel(tspan) || abs(X(end) - tspan(end)) > 1e-12*max(1, abs(tspan(end)))
        error('CPDEM:integrationFailed', ...
            ['The ODE solver did not reach the requested end point ', ...
             '(stopped at x = %g of %g). The requested property value may ', ...
             'be outside the range the model can reach with these ', ...
             'constituent properties, or the tolerances may be too tight.'], ...
            X(end), tspan(end));
    end

    Ys = YY(keep, :);
    [~, loc] = ismember(targets(idx), xs);
    Y(idx, :) = Ys(loc, :);
end

end
