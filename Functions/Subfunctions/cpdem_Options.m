function opts = cpdem_Options(varargin)
% Default solver/optimiser options for the Cross-Property DEM Toolbox
%
%   opts = cpdem_Options()                    default options
%   opts = cpdem_Options(s)                   defaults overridden by struct s
%   opts = cpdem_Options('Name', value, ...)  defaults overridden by pairs
%   opts = cpdem_Options(s, 'Name', value)    both
%
% Fields (defaults in brackets):
%   RelTol   [1e-10]   relative tolerance of the ODE solver
%   AbsTol   [1e-13]   absolute tolerance (the states are O(1): moduli are
%                      scaled by the background bulk modulus, conductivity
%                      is carried as a logarithm, and porosity is in [0,1))
%   MaxStep  []        maximum step in the independent variable (empty:
%                      solver default)
%   Solver   ['ode45'] name of the ODE solver (e.g. 'ode45', 'ode23',
%                      'ode113' on MATLAB)
%   Bounds   [1 100]   EPAR search interval for calibration (prolate by
%                      default; Cilli & Chapman 2021 found prolate EPARs fit
%                      sandstones best, but [0.01 100] is permissible)
%   LogSearch [true]   search for EPAR in log space (scale aware)
%   Norm     [2]       1 -> l1 misfit, 2 -> l2 misfit in calibration
%   Relative [false]   use relative (fractional) residuals in calibration
%   Weights  []        per-sample weights for calibration (empty: equal)
%   TolX     [1e-6]    fminbnd tolerance on (log) EPAR
%   Display  ['off']   'off' or 'iter' for the optimiser
%   PerSample [false]  calibrate an EPAR for every sample individually
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

opts = struct( ...
    'RelTol',    1e-10, ...
    'AbsTol',    1e-13, ...
    'MaxStep',   [], ...
    'Solver',    'ode45', ...
    'Bounds',    [1 100], ...
    'LogSearch', true, ...
    'Norm',      2, ...
    'Relative',  false, ...
    'Weights',   [], ...
    'TolX',      1e-6, ...
    'Display',   'off', ...
    'PerSample', false);

args = varargin;
if ~isempty(args) && (isstruct(args{1}) || isempty(args{1}))
    s = args{1};
    args(1) = [];
    if isstruct(s)
        fn = fieldnames(s);
        for ii = 1:numel(fn)
            opts = setopt(opts, fn{ii}, s.(fn{ii}));
        end
    end
end

if mod(numel(args), 2) ~= 0
    error('CPDEM:badOptions', 'Options must be given as name/value pairs.');
end
for ii = 1:2:numel(args)
    opts = setopt(opts, args{ii}, args{ii+1});
end

end

function opts = setopt(opts, name, value)
known = fieldnames(opts);
hit = strcmpi(known, name);
if ~any(hit)
    error('CPDEM:badOptions', 'Unknown option "%s". Known options: %s.', ...
        name, strjoin(known', ', '));
end
opts.(known{hit}) = value;
end
