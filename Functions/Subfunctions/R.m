function r = R(sig0,sig2,epar)
% Calculate R
%
% Calculate R from Berryman (1995) used in equations 20 & 21 of Cilli and
% Chapman (2021):
%
%   R = (1/9) sum_j 1/(sig0 + L_j (sig2 - sig0)) = m-bar/(3 sig0).
%
% The electrical DEM written with R,
%   (1 - phi) d sig/d phi = 3 sig (sig2 - sig) R,
% is identical to the m-bar form (see mbar); both are provided so that the
% two forms of the cross-property equations in Cilli and Chapman (2021)
% can be evaluated.
%
% sig0 may be an array; sig2 and EPAR are scalars.
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Osborn, J. A. (1945). Demagnetizing factors of the general ellipsoid.
% Physical Review, 67(11-12), 351, DOI:10.1103/PhysRev.67.351
%
% Berryman, J. G. (1995). Mixture theories for rock properties. Rock
% physics and phase relations: A handbook of physical constants, 3,
% 205-228. DOI:10.1029/RF003p0205
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Revised for v2.0 (sphere-safe L, vectorised).

r = mbar(sig0, sig2, epar)./(3*sig0);

end
