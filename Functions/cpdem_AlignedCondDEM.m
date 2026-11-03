function [sigPar, sigPerp] = cpdem_AlignedCondDEM(phi, sig1, sig2, epar, opts)
% Anisotropic conductivity DEM for aligned spheroidal inclusions
%
%   [sigPar, sigPerp] = cpdem_AlignedCondDEM(phi, sig1, sig2, epar)
%   [sigPar, sigPerp] = cpdem_AlignedCondDEM(phi, sig1, sig2, epar, opts)
%
% Integrates the differential effective medium equations for a host
% containing spheroidal inclusions of aspect ratio EPAR whose unique axes
% are all parallel, so that the effective medium is transversely
% isotropic with conductivity sigPar along the common axis and sigPerp in
% the plane normal to it:
%
%   (1 - phi) d sig_i/d phi = sig_i (sig2_i - sig_i) / (sig_i + L_i(EPAR') (sig2_i - sig_i)),
%   EPAR' = EPAR sqrt(sig_perp/sig_par),  i = par, perp,
%
% with the depolarisation factors evaluated in the affinely stretched
% frame in which the current host is isotropic (see
% cpdem_AlignedIncrement). This applies to any property obeying Laplace's
% equation (electrical or thermal conductivity, permittivity, ...).
%
% Inputs
%   phi  : inclusion volume fraction, scalar or array, 0 <= phi < 1
%   sig1 : host conductivity: scalar (isotropic) or [sig1_par sig1_perp]
%   sig2 : inclusion conductivity: scalar, or [sig2_axis sig2_plane] for
%          an anisotropic inclusion material (e.g. graphite flakes)
%   epar : aspect ratio (scalar) or function handle of inclusion fraction
%   opts : optional cpdem_Options struct
%
% Outputs (shape of phi): the two principal effective conductivities.
%
% This is an extension of the toolbox beyond Cilli and Chapman (2021,
% arXiv:2102.12434), whose model assumes randomly oriented inclusions; it
% uses only the classical affine-transformation argument for ellipsoidal
% inclusions in an anisotropic matrix.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 5 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
cpdem_CheckPositive('sig1', sig1, 'sig2', sig2, 'phi>=0', phi);
if any(phi(:) >= 1)
    error('CPDEM:badInput', 'Inclusion fraction must be < 1.');
end
[s1par, s1perp] = splitHost(sig1);

shape = size(phi);
Y = cpdem_Integrate(@(p, y) deriv(p, y, s1par, s1perp, sig2, epar), phi(:), [0; 0], opts);
sigPar  = reshape(s1par*exp(Y(:, 1)), shape);
sigPerp = reshape(s1perp*exp(Y(:, 2)), shape);

end

function dy = deriv(phi, y, s1par, s1perp, sig2, epar)
spar  = s1par*exp(y(1));
sperp = s1perp*exp(y(2));
[Gpar, Gperp] = cpdem_AlignedIncrement(spar, sperp, sig2, cpdem_EvalEPAR(epar, phi));
dy = [Gpar/spar; Gperp/sperp]/(1 - phi);
end

function [spar, sperp] = splitHost(sig1)
if isscalar(sig1)
    spar = sig1; sperp = sig1;
elseif numel(sig1) == 2
    spar = sig1(1); sperp = sig1(2);
else
    error('CPDEM:badInput', 'Host conductivity must be a scalar or [sig_par sig_perp].');
end
end
