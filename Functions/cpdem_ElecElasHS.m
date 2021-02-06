function [k_Up,k_Low,mu_Up,mu_Low] = cpdem_ElecElasHS(sig,sigSol,sigFl,kSol,kFl,muSol)
% Calculate electrical-elastic Hashin Shtrikman bounds
%
% Calculate electrical-elastic Hashin Shtrikman bounds using the method of
% Carcione et al. (2007). Note the method here
% sets the lower shear modulus bound is zero as the fluid phase is assumed
% to be inviscid (Carcione et al. 2007).
% 
% Reference: 
% Carcione, J. M., Ursin, B., & Nordskag, J. I. (2007). Cross-property
% relations between electrical conductivity and the seismic velocity of
% rocks. Geophysics, 72(5), E193-E204.
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

phiUp  = (sigSol-sig)./(sigSol-sigFl).*(3.*sigFl)./(sig+2.*sigFl);
phiLow = (sigSol-sig)./(sigSol-sigFl).*(sigFl+2.*sigSol)./(sig+2.*sigSol);
iden   = ones(size(sig));

%Bulk modulus bounds
x     = 4./3.*muSol;
k_Up  = iden./((iden-phiUp)./(kSol+x)+phiUp./(kFl+x))-x;
k_Low = iden./(phiLow./kFl+(iden-phiLow)./kSol);

%Shear modulus bounds
x      = muSol./6 .* ((9.*kSol + 8.*muSol)./(kSol + 2.*muSol));
mu_Up  = iden./((iden - phiUp)./(muSol + x) + phiUp./x) - x;
mu_Low = zeros(size(sig)); %Set to zero assuming inviscid fluid phase

end

