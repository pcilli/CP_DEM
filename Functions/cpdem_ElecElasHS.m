function [k_Up,k_Low,mu_Up,mu_Low,phi_Low,phi_Up] = cpdem_ElecElasHS(sig,sigSol,sigFl,kSol,kFl,muSol)
% Calculate electrical-elastic Hashin-Shtrikman bounds
%
%   [k_Up,k_Low,mu_Up,mu_Low] = cpdem_ElecElasHS(sig,sigSol,sigFl,kSol,kFl,muSol)
%   [k_Up,k_Low,mu_Up,mu_Low,phi_Low,phi_Up] = cpdem_ElecElasHS(...)
%
% Calculate electrical-elastic Hashin-Shtrikman bounds using the method of
% Carcione et al. (2007). The measured conductivity is first bracketed by
% the Hashin-Shtrikman conductivity bounds, which are inverted for the
% minimum and maximum porosity consistent with it,
%
%   phi_Low = (sigSol - sig)/(sigSol - sigFl) * 3 sigFl/(sig + 2 sigFl)
%   phi_Up  = (sigSol - sig)/(sigSol - sigFl) * (sigFl + 2 sigSol)/(sig + 2 sigSol)
%
% (from the upper and lower conductivity bounds respectively), and the
% elastic Hashin-Shtrikman bounds are then evaluated at those porosities:
% the upper elastic bound at phi_Low and the lower elastic bound at phi_Up.
% The fluid phase is assumed inviscid (muFl = 0), so the lower shear
% modulus bound is zero and the lower bulk modulus bound is the Reuss
% average.
%
% Inputs (sig may be an array; the others are scalars)
%   sig    : measured effective conductivity, sigSol < sig < sigFl
%   sigSol : solid (mineral) conductivity
%   sigFl  : fluid conductivity (must exceed sigSol)
%   kSol   : solid bulk modulus
%   kFl    : fluid bulk modulus
%   muSol  : solid shear modulus
%
% Outputs (shape of sig)
%   k_Up,  k_Low  : upper and lower bounds on the bulk modulus
%   mu_Up, mu_Low : upper and lower bounds on the shear modulus (mu_Low = 0)
%   phi_Low, phi_Up : bounds on porosity implied by the conductivity
%
% Any cross-property DEM prediction must lie within these bounds; see the
% tests for an automated check.
%
% Reference:
% Carcione, J. M., Ursin, B., & Nordskag, J. I. (2007). Cross-property
% relations between electrical conductivity and the seismic velocity of
% rocks. Geophysics, 72(5), E193-E204, DOI:10.1190/1.2762224
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Revised for v2.0 (validation, porosity outputs).

cpdem_CheckPositive('sig', sig, 'sigSol', sigSol, 'sigFl', sigFl, ...
    'kSol', kSol, 'kFl', kFl, 'muSol', muSol);
if ~all([isscalar(sigSol), isscalar(sigFl), isscalar(kSol), isscalar(kFl), isscalar(muSol)])
    error('CPDEM:badInput', 'Constituent properties must be scalars.');
end
if ~(sigFl > sigSol)
    error('CPDEM:badInput', ...
        'cpdem_ElecElasHS assumes a fluid more conductive than the solid (sigFl > sigSol).');
end
if any(sig(:) < sigSol | sig(:) > sigFl)
    warning('CPDEM:outOfRange', ...
        'Some conductivities lie outside [sigSol, sigFl]; the porosity bounds are not meaningful there.');
end

% Porosity bounds from the Hashin-Shtrikman conductivity bounds
phi_Low = (sigSol - sig)./(sigSol - sigFl).*(3.*sigFl)./(sig + 2.*sigFl);
phi_Up  = (sigSol - sig)./(sigSol - sigFl).*(sigFl + 2.*sigSol)./(sig + 2.*sigSol);
iden    = ones(size(sig));

% Bulk modulus bounds
x     = 4./3.*muSol;
k_Up  = iden./((iden - phi_Low)./(kSol + x) + phi_Low./(kFl + x)) - x;
k_Low = iden./(phi_Up./kFl + (iden - phi_Up)./kSol);

% Shear modulus bounds
x      = muSol./6 .* ((9.*kSol + 8.*muSol)./(kSol + 2.*muSol));
mu_Up  = iden./((iden - phi_Low)./(muSol + x) + phi_Low./x) - x;
mu_Low = zeros(size(sig)); % Inviscid fluid phase

end
