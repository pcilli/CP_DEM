function yprime=cpdem_Deriv(t,y,k2,mu2,sig2,epar)
% Derivative used in cpdem_Forward
%
% The derivative terms in equations 18 & 19, or 20 & 21 of Cilli and
% Chapman (2021).
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
%
% Berryman, J. G. (1980). Long-wavelength propagation in composite elastic
% media II. Ellipsoidal inclusions. The Journal of the Acoustical Society
% of America, 68(6), 1820-1831, DOI:10.1121/1.385172
%
% Berryman, J. G. (1995). Mixture theories for rock properties. Rock
% physics and phase relations: A handbook of physical constants, 3,
% 205-228. DOI:10.1029/RF003
% 
% Written by Phil Cilli, January 2021 as a part of CP_DEM Version 1.0

k0  = y(1);
mu0 = y(2);

p  = P(k0,mu0,k2,mu2,epar); %From Berryman (1980)
q  = Q(k0,mu0,k2,mu2,epar); %From Berryman (1980)
mb = mbar(t,sig2,epar);     %From Cilli and Chapman (2021)

yprime(1,1) = (k2-k0)*p/((sig2-t)*mb);
yprime(2,1) = (mu2-mu0)*q/((sig2-t)*mb);

% To do the same calculation using R from Berryman (1995) instead of m-bar:
% 
% r = R(t,sig2,epar); %From Berryman (1995)
% yprime(1) = 1/(3*t)*(k2-k0)*p/((sig2-t)*r);
% yprime(2) = 1/(3*t)*(mu2-mu0)*q/((sig2-t)*r);

end
