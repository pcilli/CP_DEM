function out = cpdem_AnisoInverse(cMeas, which, k1, mu1, sig1, k2, sig2, alpha, S, varargin)
% Anisotropic cross-property DEM: conductivity tensor from one TI stiffness
%
%   out = cpdem_AnisoInverse(C33, 'C33', k1, mu1, sig1, k2, sig2, alpha)
%   out = cpdem_AnisoInverse(cMeas, which, k1, mu1, sig1, k2, sig2, alpha, S, opts)
%
% Inverse problem of Cilli, "Anisotropic cross-property differential
% effective medium modelling" (section 2.7): the chain rule is inverted and
% a monotonic stiffness component is used as the independent variable, e.g.
% eq. 63,  dsig_par/dC33 = (sig2 - sig_par) mbar_par / M33,  with the
% remaining stiffnesses, the other conductivity component and the
% inclusion fraction integrated alongside from the matrix values.
%
% Inputs
%   cMeas  : measured values of the chosen stiffness (Pa), array
%   which  : 'C11', 'C33', 'C44' or 'C66' (all decrease monotonically with
%            inclusion fraction for fluid inclusions)
%   k1, mu1, sig1, k2, sig2, alpha, S, opts : as in cpdem_AnisoForward
%
% Output struct (fields shaped like cMeas): sigPar, sigPerp, C11, C12,
% C13, C33, C44, C66 and phi. Measurements outside the range between the
% matrix value and the fluid limit (k2 for C11/C33, 0 for C44/C66) are
% unreachable and return NaN with a warning.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

if nargin < 9 || isempty(S), S = 1; end
[opts, closure] = cpdem_AnisoOptions(varargin);
cpdem_CheckPositive('cMeas>=0', cMeas, 'k1', k1, 'mu1', mu1, 'sig1', sig1, 'k2', k2, 'sig2>=0', sig2);
cpdem_AnisoCheck(alpha, S, closure);
idx = find(strcmpi(which, {'C11', 'C13', 'C33', 'C44', 'C66'}));
if isempty(idx) || idx == 2
    error('CPDEM:badInput', 'which must be ''C11'', ''C33'', ''C44'' or ''C66''.');
end

[~, ~, L] = cpdem_Shape(alpha);
prm = struct('sig2', sig2, 'k2', k2, 'alpha', alpha, 'S', S, 'L', L, 'closure', closure);
lam = k1 - 2*mu1/3;
c0 = [lam + 2*mu1; lam; lam + 2*mu1; mu1; mu1];
cStart = c0(idx);
if idx <= 3, cEnd = k2; else, cEnd = 0; end

shape = size(cMeas);
n = numel(cMeas);
cMeas = cMeas(:);
bad = ~(cMeas <= cStart & cMeas > cEnd);
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d stiffnesses lie outside the range between the matrix value and the fluid limit; NaN returned for them.', nnz(bad), n);
end

% state y = [ln sigPar/sig1; ln sigPerp/sig1; other four c/k1; phi];
% independent variable x = (C_which - C_which,0)/k1
others = setdiff(1:5, idx);
Y = nan(n, 7);
if any(~bad)
    fun = @(x, y) deriv(x, y, prm, idx, others, cStart, sig1, k1);
    Y(~bad, :) = cpdem_Integrate(fun, (cMeas(~bad) - cStart)/k1, [0; 0; c0(others)/k1; 0], opts);
end

c5 = nan(n, 5);
c5(:, others) = Y(:, 3:6)*k1;
c5(:, idx) = cMeas; c5(bad, idx) = NaN;
out.sigPar  = reshape(sig1*exp(Y(:,1)), shape);
out.sigPerp = reshape(sig1*exp(Y(:,2)), shape);
out.C11 = reshape(c5(:,1), shape); out.C13 = reshape(c5(:,2), shape);
out.C33 = reshape(c5(:,3), shape); out.C44 = reshape(c5(:,4), shape);
out.C66 = reshape(c5(:,5), shape); out.C12 = out.C11 - 2*out.C66;
out.phi = reshape(Y(:,7), shape);

end

function dy = deriv(x, y, prm, idx, others, cStart, sig1, k1)
c5 = zeros(5, 1);
c5(idx) = cStart + x*k1;
c5(others) = y(3:6)*k1;
sp = sig1*exp(y(1)); sq = sig1*exp(y(2));
phi = y(7);
r = cpdem_AnisoRates(sp, sq, c5, prm);
rate = r(2 + idx)/k1;              % (1 - phi) dx/dphi
dy = [r(1)/sp; r(2)/sq; r(2 + others)/k1; (1 - phi)]/rate;
end
