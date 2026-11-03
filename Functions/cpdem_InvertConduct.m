function sigInv = cpdem_InvertConduct(k0,mu0,k1,mu1,sig1,k2,mu2,sig2,epar,lb,ub,modFlag,normFlag,paramFlag,opts)
% Invert for electrical conductivity from elastic measurements
%
%   sigInv = cpdem_InvertConduct(k0,mu0,k1,mu1,sig1,k2,mu2,sig2,epar,lb,ub,modFlag,normFlag,paramFlag)
%
% Numerically invert equations 18 & 19, or 20 & 21 of Cilli and Chapman
% (2021) for conductivity. This is the v1.0 interface, retained for
% backward compatibility. Since v2.0 the inversion is performed directly by
% cpdem_Inverse (integrating the cross-property equations with the modulus
% as the independent variable), which is exact to solver tolerance rather
% than a bounded search; results are then clipped to [lb, ub] and
% measurements outside the range of the model are returned at the nearest
% bound, as the v1.0 search would have done.
%
% Inputs
%   k0, mu0       : measured bulk and shear moduli (scalar or vectors)
%   k1, mu1, sig1 : background (mineral) properties
%   k2, mu2, sig2 : inclusion (fluid) properties
%   epar          : EPAR (use alpha_K with modFlag 0, alpha_mu with 1)
%   lb, ub        : conductivity bounds (typically sig1 and sig2)
%   modFlag       : 0 -> invert the bulk modulus (sigma_K);
%                   1 -> invert the shear modulus (sigma_mu)
%   normFlag      : accepted for compatibility (each sample is inverted
%                   exactly, so no norm is involved)
%   paramFlag     : 0 -> scalar parameters, 1 -> per-sample parameter
%                   vectors (now detected automatically)
%   opts          : optional cpdem_Options struct
%
% Output: inverted conductivity, same size as the measured modulus.
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Revised for v2.0.

if nargin < 14
    error('CPDEM:badInput', 'cpdem_InvertConduct needs 14 inputs (see help).');
end
if nargin < 15 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end
if ~any(normFlag == [1 2])
    error('CPDEM:badInput', 'Invalid normFlag. Set to 1 or 2.');
end
if ~any(paramFlag == [0 1])
    error('CPDEM:badInput', 'Invalid paramFlag. Set to 0 or 1.');
end
useK = cpdem_WhichModulus(modFlag);

if useK
    m0 = k0; m1 = k1; m2 = k2;
else
    m0 = mu0; m1 = mu1; m2 = mu2;
end

ws = warning('off', 'CPDEM:outOfRange');
sigInv = cpdem_Inverse(m0, useK, k1, mu1, sig1, k2, mu2, sig2, epar, opts);
warning(ws);

% Emulate the bounded search of v1.0 for measurements the model cannot
% reach: a modulus at or beyond the background value maps to the
% background-side bound, one at or beyond the inclusion value to the
% inclusion-side bound.
n = numel(m0);
[m1v, m2v] = cpdem_Broadcast(n, m1, m2);
m0v = m0(:);
out = isnan(sigInv(:));
if any(out)
    atBackground = out & (abs(m0v - m1v) <= abs(m0v - m2v));
    sigInv(reshape(atBackground, size(sigInv))) = lb;
    sigInv(reshape(out & ~atBackground, size(sigInv))) = ub;
end
sigInv = min(max(sigInv, lb), ub);

end
