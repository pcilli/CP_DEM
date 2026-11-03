function [sigPar, sigPerp] = cpdem_AnisoElecDEM(phi, sig1, sig2, alpha, S, varargin)
% Transversely isotropic electrical DEM of the anisotropic cross-property paper
%
%   [sigPar, sigPerp] = cpdem_AnisoElecDEM(phi, sig1, sig2, alpha)
%   [sigPar, sigPerp] = cpdem_AnisoElecDEM(phi, sig1, sig2, alpha, S, opts)
%   [sigPar, sigPerp] = cpdem_AnisoElecDEM(..., S, 'Name', value, ...)
%
% Integrates eqs 11-12 of Cilli, "Anisotropic cross-property differential
% effective medium modelling" (aligned spheroids of aspect ratio alpha in an
% isotropic matrix, isotropic-host depolarisation closure), with the
% partial-alignment surrogate of eqs 54-56 for an order parameter S:
%
%   (1 - phi) dsig_par /dphi = (sig2 - sig_par ) sig_par /(sig_par  + L (sig2 - sig_par)),
%   (1 - phi) dsig_perp/dphi = (sig2 - sig_perp) sig_perp/(sig_perp + (1-L)/2 (sig2 - sig_perp)),
%
% from sig_par = sig_perp = sig1 at phi = 0. S = 1 (default) is perfect
% alignment; S = 0 reproduces the isotropic randomly oriented DEM
% (cpdem_ElecDEM) exactly.
%
% Inputs: phi (array, 0 <= phi < 1), sig1, sig2 (scalars), alpha (scalar
% aspect ratio c/a), S in [0, 1], opts: cpdem_Options struct or name/value
% pairs; the extra option 'Closure' selects 'isotropic-host' (the paper's
% closure, default) or 'exact' (the exact coaxial solution via the affine
% stretch, equivalent to cpdem_AlignedCondDEM; S must then be 1).
%
% Outputs (shape of phi): axial and transverse conductivities.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

if nargin < 5 || isempty(S), S = 1; end
[opts, closure] = cpdem_AnisoOptions(varargin);
cpdem_CheckPositive('sig1', sig1, 'sig2>=0', sig2, 'phi>=0', phi);
cpdem_AnisoCheck(alpha, S, closure);
if any(phi(:) >= 1), error('CPDEM:badInput', 'Porosity must be < 1.'); end

[~, ~, L] = cpdem_Shape(alpha);
prm = struct('sig2', sig2, 'alpha', alpha, 'S', S, 'L', L, 'closure', closure);
fun = @(p, y) elecOnly(p, y, prm);
Y = cpdem_Integrate(fun, phi(:), [sig1; sig1], opts);
sigPar  = reshape(Y(:, 1), size(phi));
sigPerp = reshape(Y(:, 2), size(phi));

end

function dy = elecOnly(phi, y, prm)
sp = y(1); sq = y(2);
switch prm.closure
    case 'isotropic-host'
        L = prm.L;
        mPar  = sp/(sp + L*(prm.sig2 - sp));
        mPerp = sq/(sq + (1 - L)/2*(prm.sig2 - sq));
        mRand = (2*mPerp + mPar)/3;
        dy = [(prm.sig2 - sp)*(mRand + prm.S*(mPar - mRand));
              (prm.sig2 - sq)*(mRand + prm.S*(mPerp - mRand))]/(1 - phi);
    case 'exact'
        [Gpar, Gperp] = cpdem_AlignedIncrement(sp, sq, prm.sig2, prm.alpha);
        dy = [Gpar; Gperp]/(1 - phi);
end
end
