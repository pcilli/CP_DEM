function mb = mbar(sig0,sig2,epar)
% Calculate m-bar
%
% Calculate m-bar in equations 18 & 19 of Cilli and Chapman (2021):
%
%   m-bar = (1/3) sig0 * sum_j 1/(sig0 + L_j (sig2 - sig0)),   j = 1,2,3
%
% where L_1 = L is the depolarisation factor along the unique axis of the
% spheroidal inclusion (Osborn 1945) and L_2 = L_3 = (1 - L)/2. m-bar is
% the geometric factor of the electrical DEM,
%
%   (1 - phi) d sig/d phi = (sig2 - sig) m-bar(sig, sig2, EPAR),
%
% and generalises Archie's cementation exponent: for insulating inclusions
% (sig2 = 0) it is constant, m-bar = (1/3) sum_j 1/(1 - L_j), which is 3/2
% for spheres (Sen et al. 1981; Mendelson & Cohen 1982), and the DEM then
% integrates to sig = sig0 (1 - phi)^m-bar. See also cpdem_Shape, R.
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
% Sen, P. N., Scala, C., and Cohen, M. H. (1981). A self-similar model for
% sedimentary rocks with application to the dielectric constant of fused
% glass beads. Geophysics, 46(5), 781-795, DOI:10.1190/1.1441215
%
% Mendelson, K. S., and Cohen, M. H. (1982). The effect of grain anisotropy
% on the electrical properties of sedimentary rocks. Geophysics, 47(2),
% 257-263, DOI:10.1190/1.1441332
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0. Revised for v2.0 (sphere-safe L, vectorised).

[~, ~, L] = cpdem_Shape(epar);

% m-bar from Cilli and Chapman (2021). The second term collects the two
% degenerate axes, each with depolarisation factor (1 - L)/2.
mb = (1/3)*sig0.*(1./(sig0 - L*(sig0 - sig2)) + ...
     4./(sig0 + sig2 + L*(sig0 - sig2)));

end
