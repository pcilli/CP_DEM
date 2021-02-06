function [k0, mu0] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0,epar)
% Forward model elastic moduli from electrical resistivity
%
% Numerically evaluate equations 18 & 19, or 20 & 21 of Cilli and
% Chapman (2021).
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

global inVar;

inVar = zeros(4,1);
inVar(1) = k2; 
inVar(2) = mu2;
inVar(3) = sig2;
inVar(4) = epar;

%Optimisation routine only tested on MATLAB R2019B.
[tout,yout] = ode45('cpdem_Deriv',[sig1; sig0],[k1; mu1],1e-5);

k0  = yout(length(tout),1);
mu0 = yout(length(tout),2);

end