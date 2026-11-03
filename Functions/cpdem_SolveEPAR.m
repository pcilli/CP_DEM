function [epar, misfit, cal] = cpdem_SolveEPAR(k0,mu0,sig0,k1,mu1,sig1,k2,mu2,sig2,modFlag,normFlag,paramFlag,lb,ub,opts)
% Solve for EPAR using elastic moduli and electrical conductivity
%
%   epar = cpdem_SolveEPAR(k0,mu0,sig0,k1,mu1,sig1,k2,mu2,sig2,modFlag,normFlag,paramFlag,lb,ub)
%   [epar, misfit, cal] = cpdem_SolveEPAR(..., lb, ub, opts)
%
% Solve equations 18 & 19, or 20 & 21 of Cilli and Chapman (2021) for EPAR
% with all other parameters and measurements known. This is the v1.0
% interface, retained for backward compatibility; it now calls
% cpdem_CalibrateEPAR, which is much faster (one ODE integration per trial
% EPAR for the whole data set) and more accurate (tight solver tolerances,
% search in log EPAR).
%
% Inputs
%   k0, mu0, sig0 : measured bulk modulus, shear modulus and conductivity
%   k1, mu1, sig1 : background (mineral) properties
%   k2, mu2, sig2 : inclusion (fluid) properties
%   modFlag       : 0 -> solve for the bulk modulus EPAR (alpha_K)
%                   1 -> solve for the shear modulus EPAR (alpha_mu)
%   normFlag      : 1 -> l1 norm of residuals, 2 -> l2 norm
%   paramFlag     : 0 -> parameters are scalars; 1 -> parameters are
%                   per-sample vectors (now detected automatically; the
%                   flag is accepted for compatibility)
%   lb, ub        : EPAR search bounds
%   opts          : optional cpdem_Options struct
%
% Outputs
%   epar   : best-fitting EPAR
%   misfit : minimised misfit
%   cal    : full calibration struct from cpdem_CalibrateEPAR
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Revised for v2.0.

if nargin < 14
    error('CPDEM:badInput', 'cpdem_SolveEPAR needs 14 inputs (see help).');
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

opts.Bounds = [lb ub];
opts.Norm = normFlag;

if useK
    cal = cpdem_CalibrateEPAR(k0, [], sig0, k1, mu1, sig1, k2, mu2, sig2, opts);
    epar = cal.eparK;
    misfit = cal.misfitK;
else
    cal = cpdem_CalibrateEPAR([], mu0, sig0, k1, mu1, sig1, k2, mu2, sig2, opts);
    epar = cal.eparMu;
    misfit = cal.misfitMu;
end

end
