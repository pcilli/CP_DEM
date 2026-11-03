function [k0, mu0, phi0] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0,epar,opts)
% Forward model elastic moduli from electrical conductivity
%
%   [k0, mu0]       = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0,epar)
%   [k0, mu0, phi0] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0,epar,opts)
%
% Numerically integrates the cross-property DEM equations of Cilli and
% Chapman (2021) (equations 18 & 19, written with m-bar; equivalently
% equations 20 & 21 written with Berryman's R) from the background
% (mineral) properties at sig = sig1 to the measured conductivity sig0:
%
%   dK /d sig = (K2  - K ) P / ((sig2 - sig) m-bar)
%   dmu/d sig = (mu2 - mu) Q / ((sig2 - sig) m-bar)
%
% with K(sig1) = k1 and mu(sig1) = mu1. Porosity is a dummy variable of the
% model and plays no part in the result. As a diagnostic (not part of the
% paper's model) the electrical DEM porosity is integrated alongside,
% dphi/dsig = (1 - phi)/((sig2 - sig) m-bar), and returned as the optional
% third output: it is the porosity the DEM path passes through to reach a
% measurement, which need not equal the true porosity of a real rock.
%
% Inputs (SI units are recommended: moduli in Pa, conductivity in S/m)
%   k1, mu1, sig1 : background (mineral) bulk modulus, shear modulus and
%                   conductivity. sig1 must be > 0: DEM with a perfectly
%                   insulating host never conducts.
%   k2, mu2, sig2 : inclusion (pore fluid) bulk modulus, shear modulus
%                   (0 for a fluid) and conductivity.
%   sig0          : measured effective conductivity; scalar or array of any
%                   shape. Values must lie between sig1 and sig2.
%   epar          : equivalent pore aspect ratio (EPAR): unique axis length
%                   over degenerate axis length (oblate < 1 < prolate). It
%                   may be a scalar, an array matching sig0 (per-sample
%                   EPAR), or a function handle of porosity, epar(phi), e.g.
%                   the power law of Cilli & Chapman (2020).
%   opts          : optional cpdem_Options struct (solver tolerances, and
%                   Form = 'mbar' or 'R' to evaluate equations 18 & 19 or
%                   20 & 21 of the paper; they are identical).
%
% Each of k1, mu1, sig1, k2, mu2, sig2 may also be an array matching sig0
% (sample-dependent parameters); this replaces v1.0's "paramFlag".
%
% Outputs have the shape of sig0. When every parameter is shared by all
% measurements the whole set is computed with a single ODE integration.
%
% Note: a different EPAR (the bulk modulus EPAR, alpha_K, or the shear
% modulus EPAR, alpha_mu) is in general needed to predict each modulus
% accurately for real rocks (Cilli & Chapman 2021). With alpha_K only k0 is
% physically meaningful; with alpha_mu only mu0 is. See
% cpdem_CalibrateEPAR and cpdem_ModelVpVs.
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Rewritten for v2.0 (vectorised, tight tolerances,
% integration in log conductivity, porosity output, EPAR(phi)).

if nargin < 8
    error('CPDEM:badInput', 'cpdem_Forward needs at least 8 inputs.');
end
if nargin < 9 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end

cpdem_CheckPositive('k1', k1, 'mu1', mu1, 'sig1', sig1, 'k2', k2, ...
    'mu2>=0', mu2, 'sig2', sig2, 'sig0', sig0);

shape = size(sig0);
n = numel(sig0);
sig0 = sig0(:);
[k1v, mu1v, sig1v, k2v, mu2v, sig2v] = cpdem_Broadcast(n, k1, mu1, sig1, k2, mu2, sig2);

% Per-sample EPAR?
if isnumeric(epar) && ~isscalar(epar)
    if numel(epar) ~= n
        error('CPDEM:badEPAR', ...
            'EPAR must be a scalar, a function handle, or have one value per measurement.');
    end
    eparv = epar(:);
    perSampleEPAR = true;
else
    perSampleEPAR = false;
end

% Measurements outside (sig1, sig2) cannot be reached by the model.
bad = ~cpdem_InRange(sig0, sig1v, sig2v);
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d conductivities lie outside the interval (sig1, sig2) spanned by the model; NaN returned for them.', ...
        nnz(bad), n);
end

k0   = nan(n, 1);
mu0  = nan(n, 1);
phi0 = nan(n, 1);

shared = all([isscalar(k1), isscalar(mu1), isscalar(sig1), isscalar(k2), ...
              isscalar(mu2), isscalar(sig2)]) && ~perSampleEPAR;

if shared
    ok = ~bad;
    if any(ok)
        u = log(sig0(ok)/sig1);
        Y = cpdem_Integrate(@(u, y) cpdem_Deriv(u, y, k1, k2, mu2, sig1, sig2, epar, opts.Form), ...
                            u, [1; mu1/k1; 0], opts);
        k0(ok)   = Y(:, 1)*k1;
        mu0(ok)  = Y(:, 2)*k1;
        phi0(ok) = Y(:, 3);
    end
else
    for ii = find(~bad).'
        if perSampleEPAR
            e = eparv(ii);
        else
            e = epar;
        end
        u = log(sig0(ii)/sig1v(ii));
        Y = cpdem_Integrate(@(u, y) cpdem_Deriv(u, y, k1v(ii), k2v(ii), mu2v(ii), ...
                                               sig1v(ii), sig2v(ii), e, opts.Form), ...
                            u, [1; mu1v(ii)/k1v(ii); 0], opts);
        k0(ii)   = Y(1)*k1v(ii);
        mu0(ii)  = Y(2)*k1v(ii);
        phi0(ii) = Y(3);
    end
end

k0   = reshape(k0, shape);
mu0  = reshape(mu0, shape);
phi0 = reshape(phi0, shape);

end
