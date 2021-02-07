function epar = cpdem_SolveEPAR(k0,mu0,sig0,k1,mu1,sig1,k2,mu2,sig2,modFlag,normFlag,paramFlag,lb,ub)
% Solve for EPAR using elastic moduli and electrical resistivity
%
% Solve equations 18 & 19, or 20 & 21 of Cilli and Chapman (2021) for EPAR
% with all other parameters and measurements known.
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

ObjFn = @(x) cpdem_SolveEPARCost(x,k0,mu0,sig0,k1,mu1,sig1,k2,mu2,sig2,modFlag,normFlag,paramFlag);

x = fminbnd(ObjFn,lb,ub); 

epar = x;

end
