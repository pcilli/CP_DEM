function p = P(k0,mu0,k2,mu2,epar)
% Calculate P
%
% Calculate P from Berryman (1980) [DOI:10.1121/1.385172] used in equations
% 18 & 19, or 20 & 21 of Cilli and Chapman (2021).
%
% Note we have defined EPAR as the ratio of the unique ellipsoidal axis
% length to the degenerate axis length, so it spans 0 to infinity for both
% types of spheroid. However, Berryman (1980) defines it as a/c, where the
% prolate ellipsoid has axis lengths a > b = c and oblate ellipsoid has a =
% b < c, so it spans 0 to 1 for both prolate and oblate. These equations
% are taken from the appendix of Berryman (1980), and so a and c are
% calculated from our definition of EPAR accordingly.
%
% Since v2.0 this is a thin wrapper around cpdem_PQ (which evaluates P and
% Q together and handles the sphere limit EPAR == 1 exactly); it is kept
% for backward compatibility.
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
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Revised for v2.0.

[p, ~] = cpdem_PQ(k0, mu0, k2, mu2, epar);

end
