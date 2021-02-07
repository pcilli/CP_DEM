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

%Optimisation routine only tested on MATLAB R2019b.
[T,Y] = ode45(@(t,y)cpdem_Deriv(t,y,k2,mu2,sig2,epar),[sig1; sig0],[k1; mu1],1e-5);

k0  = Y(length(T),1);
mu0 = Y(length(T),2);

end
