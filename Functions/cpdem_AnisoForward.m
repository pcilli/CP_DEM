function out = cpdem_AnisoForward(k1, mu1, sig1, k2, sig2, sig0, alpha, S, which, varargin)
% Anisotropic cross-property DEM: TI stiffness from one conductivity component
%
%   out = cpdem_AnisoForward(k1, mu1, sig1, k2, sig2, sig0, alpha)
%   out = cpdem_AnisoForward(k1, mu1, sig1, k2, sig2, sig0, alpha, S, which, opts)
%
% Integrates the porosity-free system of Cilli, "Anisotropic cross-property
% differential effective medium modelling" (eqs 40-44 or 46-50 together with
% eq. 52): with the measured conductivity component as the independent
% variable, the other component and the five independent transversely
% isotropic stiffnesses follow by the chain rule,
%
%   dC_IJ/dsig_p = M_IJ /((sig2 - sig_p) mbar_p),   IJ in {11, 13, 33, 44, 66},
%   dsig_q/dsig_p = (sig2 - sig_q) mbar_q /((sig2 - sig_p) mbar_p),
%
% from the isotropic matrix at sig_p = sig1. Both DEM halves use the
% isotropic-host closures of the paper (see cpdem_AnisoRates).
%
% Inputs
%   k1, mu1, sig1 : matrix bulk modulus, shear modulus (Pa), conductivity (S/m)
%   k2, sig2      : fluid bulk modulus and conductivity (mu2 = 0)
%   sig0          : measured values of the chosen conductivity component (array)
%   alpha         : spheroid aspect ratio c/a (prolate > 1, oblate < 1)
%   S             : alignment order parameter in [0, 1] (default 1)
%   which         : 'par' (default) if sig0 is the axial conductivity,
%                   'perp' if it is the transverse one
%   opts          : cpdem_Options struct or name/value pairs; the extra
%                   option 'Closure' = 'exact' uses the exact coaxial
%                   electrical solution instead of the paper's closure
%
% Output struct (fields shaped like sig0): C11, C12, C13, C33, C44, C66
% (Pa), sigPar, sigPerp (S/m) and phi, the inclusion fraction the path
% passes through (diagnostic).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

if nargin < 8 || isempty(S), S = 1; end
if nargin < 9 || isempty(which), which = 'par'; end
[opts, closure] = cpdem_AnisoOptions(varargin);
cpdem_CheckPositive('k1', k1, 'mu1', mu1, 'sig1', sig1, 'k2', k2, 'sig2>=0', sig2, 'sig0', sig0);
cpdem_AnisoCheck(alpha, S, closure);
usePar = parsePar(which);

[~, ~, L] = cpdem_Shape(alpha);
prm = struct('sig2', sig2, 'k2', k2, 'alpha', alpha, 'S', S, 'L', L, 'closure', closure);
lam = k1 - 2*mu1/3;
c0 = [lam + 2*mu1; lam; lam + 2*mu1; mu1; mu1]/k1;

shape = size(sig0);
n = numel(sig0);
sig0 = sig0(:);
lo = min(sig1, sig2); hi = max(sig1, sig2);
bad = ~(sig0 >= lo & sig0 <= hi) | sig0 == sig2;
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d conductivities lie outside the interval between the matrix and fluid values; NaN returned for them.', nnz(bad), n);
end

Y = nan(n, 7);
if any(~bad)
    % state y = [ln(sig_other/sig1); c5/k1; phi], independent variable u = ln(sig_p/sig1)
    fun = @(u, y) deriv(u, y, prm, usePar, sig1, k1);
    Y(~bad, :) = cpdem_Integrate(fun, log(sig0(~bad)/sig1), [0; c0; 0], opts);
end
other = sig1*exp(Y(:, 1));
if usePar
    sPar = sig0; sPerp = other;
else
    sPerp = sig0; sPar = other;
end
sPar(bad) = NaN; sPerp(bad) = NaN;

out.C11 = reshape(Y(:,2)*k1, shape); out.C13 = reshape(Y(:,3)*k1, shape);
out.C33 = reshape(Y(:,4)*k1, shape); out.C44 = reshape(Y(:,5)*k1, shape);
out.C66 = reshape(Y(:,6)*k1, shape); out.C12 = out.C11 - 2*out.C66;
out.sigPar = reshape(sPar, shape); out.sigPerp = reshape(sPerp, shape);
out.phi = reshape(Y(:,7), shape);

end

function dy = deriv(u, y, prm, usePar, sig1, k1)
sw = sig1*exp(u);          % the measured component
so = sig1*exp(y(1));       % the other component
if usePar
    sp = sw; sq = so;
else
    sq = sw; sp = so;
end
phi = y(7);
r = cpdem_AnisoRates(sp, sq, y(2:6)*k1, prm);
if usePar
    rate = r(1)/sp; dOther = r(2)/sq;
else
    rate = r(2)/sq; dOther = r(1)/sp;
end
dy = [dOther; r(3:7)/k1; (1 - phi)]/rate;
end

function usePar = parsePar(which)
switch lower(which)
    case {'par', 'parallel', 'axial', 'axis'}
        usePar = true;
    case {'perp', 'perpendicular', 'transverse', 'plane'}
        usePar = false;
    otherwise
        error('CPDEM:badInput', 'which must be ''par'' or ''perp''.');
end
end
