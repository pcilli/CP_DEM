function dy = cpdem_Deriv(u, y, k1, k2, mu2, sig1, sig2, epar)
% Derivative of the cross-property DEM system used in cpdem_Forward
%
% Evaluates the right-hand side of the cross-property DEM equations of
% Cilli and Chapman (2021) (their equations 18 & 19, written with m-bar;
% equivalently 20 & 21 written with Berryman's R),
%
%   dK /d sig = (K2  - K ) P / ((sig2 - sig) m-bar),
%   dmu/d sig = (mu2 - mu) Q / ((sig2 - sig) m-bar),
%
% together with the electrical DEM for the (dummy) porosity,
%
%   dphi/d sig = (1 - phi) / ((sig2 - sig) m-bar),
%
% in the scaled variables used by the solver:
%
%   independent variable  u = ln(sig/sig1)
%   state                 y = [K/k1; mu/k1; phi]
%
% Integrating in ln(sig) removes the 1/sig behaviour of the derivatives at
% low conductivity (where m-bar is proportional to sig), and scaling the
% moduli by the background bulk modulus makes the absolute tolerance
% meaningful. Both were sources of the few-percent solver noise in v1.0.
%
% EPAR may be a scalar or a function handle of porosity (cpdem_EvalEPAR).
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Berryman, J. G. (1980). Long-wavelength propagation in composite elastic
% media II. Ellipsoidal inclusions. The Journal of the Acoustical Society
% of America, 68(6), 1820-1831, DOI:10.1121/1.385172
%
% Berryman, J. G. (1995). Mixture theories for rock properties. Rock
% physics and phase relations: A handbook of physical constants, 3,
% 205-228. DOI:10.1029/RF003p0205
%
% Written by Phil Cilli, January 2021 as a part of CP_DEM Version 1.0.
% Rewritten for v2.0.

sig = sig1*exp(u);
K   = y(1)*k1;
mu  = y(2)*k1;
phi = y(3);

alpha = cpdem_EvalEPAR(epar, phi);

[p, q] = cpdem_PQ(K, mu, k2, mu2, alpha);  % From Berryman (1980)
mb     = mbar(sig, sig2, alpha);           % From Cilli and Chapman (2021)

% g = sig d(.)/d sig per unit (1 - phi): the electrical DEM "speed".
g = sig/((sig2 - sig)*mb);

dy = [ (k2 - K)*p*g/k1;
       (mu2 - mu)*q*g/k1;
       (1 - phi)*g ];

end
