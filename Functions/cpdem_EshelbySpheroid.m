function S = cpdem_EshelbySpheroid(nu, alpha)
% Eshelby tensor of a spheroid in an isotropic host (6x6, engineering Voigt)
%
%   S = cpdem_EshelbySpheroid(nu, alpha)
%
% Returns the Eshelby tensor of a spheroidal inclusion with semi-axes
% (a, a, c), aspect ratio alpha = c/a and symmetry axis x3, embedded in an
% isotropic host of Poisson's ratio nu, as the 6x6 matrix of the Voigt
% convention used throughout the anisotropic functions (order 11, 22, 33,
% 23, 13, 12; engineering shear strains, so S(4,4) = 2 S2323 etc.). This is
% the "trusted implementation of the isotropic-host spheroidal Eshelby
% tensor" that the closure of Cilli, "Anisotropic cross-property
% differential effective medium modelling" (section 2.3.4) calls for.
%
% The closed forms are those of Eshelby (1957) as tabulated by Mura (1987)
% and Tandon & Weng (1984), written with g = theta(alpha), Berryman's
% (1980) shape function. They cancel catastrophically near the sphere, so
% for |alpha - 1| < 0.02 an eighth-order Taylor series about the sphere
% (derived symbolically; truncation error below 1e-15) is used, and
% alpha = 1 returns the classical sphere tensor
%   S1111 = (7 - 5 nu)/(15 (1 - nu)),  S1122 = (5 nu - 1)/(15 (1 - nu)),
%   S1212 = (4 - 5 nu)/(15 (1 - nu)),  i.e. S(4,4) = 2 (4 - 5 nu)/(15 (1 - nu)),
% the consistency condition of eq. 25 of the paper.
%
% Verified (test_anisotropic_paper): the orientation average of the dilute
% strain-concentration tensor built from S reproduces Berryman's (1980) P
% and Q (cpdem_PQ) to machine precision for oblate, spherical and prolate
% inclusions, fluid and solid.
%
% References:
% Eshelby, J. D. (1957). Proc. R. Soc. Lond. A, 241, 376-396.
% Mura, T. (1987). Micromechanics of Defects in Solids, 2nd ed., Nijhoff.
% Tandon, G. P., and Weng, G. J. (1984). Polymer Composites, 5, 327-333.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

if ~isscalar(alpha) || ~isfinite(alpha) || alpha <= 0
    error('CPDEM:badInput', 'alpha must be a positive finite scalar.');
end
if ~isscalar(nu) || nu <= -1 || nu >= 0.5
    error('CPDEM:badInput', 'nu must be a scalar Poisson ratio in (-1, 0.5).');
end

e = alpha - 1;
if abs(e) < 0.02
    % Series about the sphere: each component times (1 - nu) is a
    % polynomial in e whose coefficients are linear in nu, (a nu + b)/den.
    T = struct( ...
        'S1111', [-5 7 15; -14 16 105; 3 -4 35; -176 268 3465; 260 -445 9009; -80 152 5005; 952 -1988 109395; -9728 22144 2078505; 4032 -9936 1616615], ...
        'S3333', [-5 7 15; 28 -32 105; -18 16 105; 352 -224 3465; -520 200 9009; 480 -64 15015; -1904 -224 109395; 19456 7168 2078505; -8064 -4992 1616615], ...
        'S1122', [5 -1 15; 14 -4 105; -9 2 105; 176 -28 3465; -260 25 9009; 240 -8 15015; -952 -28 109395; 9728 896 2078505; -4032 -624 1616615], ...
        'S1133', [5 -1 15; 14 2 105; -9 1 105; 176 -64 3465; -260 160 9009; 240 -208 15015; -952 1064 109395; 9728 -13312 2078505; -4032 6528 1616615], ...
        'S3311', [5 -1 15; -28 2 105; 18 1 105; -352 -64 3465; 520 160 9009; -480 -208 15015; 1904 1064 109395; -19456 -13312 2078505; 8064 6528 1616615], ...
        'S1212', [-5 4 15; -14 10 105; 9 -7 105; -176 148 3465; 260 -235 9009; -240 232 15015; 952 -980 109395; -9728 10624 2078505; 4032 -4656 1616615], ...
        'S1313', [-5 4 15; 7 -5 105; -9 11 210; 88 -152 3465; -130 290 9009; 120 -328 15015; -476 1540 109395; 4864 -18176 2078505; -2016 8544 1616615]);
    S1111 = horner(T.S1111, nu, e); S3333 = horner(T.S3333, nu, e);
    S1122 = horner(T.S1122, nu, e); S1133 = horner(T.S1133, nu, e);
    S3311 = horner(T.S3311, nu, e); S1212 = horner(T.S1212, nu, e);
    S1313 = horner(T.S1313, nu, e);
else
    q = alpha;
    g = cpdem_Shape(alpha);          % Berryman's theta == Mura's g
    d = q^2 - 1;
    S1111 = (3/(8*(1 - nu)))*q^2/d + (1/(4*(1 - nu)))*(1 - 2*nu - 9/(4*d))*g;
    S3333 = (1/(2*(1 - nu)))*(1 - 2*nu + (3*q^2 - 1)/d - (1 - 2*nu + 3*q^2/d)*g);
    S1122 = (1/(4*(1 - nu)))*(q^2/(2*d) - (1 - 2*nu + 3/(4*d))*g);
    S1133 = -(1/(2*(1 - nu)))*q^2/d + (1/(4*(1 - nu)))*(3*q^2/d - (1 - 2*nu))*g;
    S3311 = -(1/(2*(1 - nu)))*(1 - 2*nu + 1/d) + (1/(2*(1 - nu)))*(1 - 2*nu + 3/(2*d))*g;
    S1212 = (1/(4*(1 - nu)))*(q^2/(2*d) + (1 - 2*nu - 3/(4*d))*g);
    S1313 = (1/(4*(1 - nu)))*(1 - 2*nu - (q^2 + 1)/d - 0.5*(1 - 2*nu - 3*(q^2 + 1)/d)*g);
end

S = zeros(6);
S(1,1) = S1111; S(2,2) = S1111; S(3,3) = S3333;
S(1,2) = S1122; S(2,1) = S1122;
S(1,3) = S1133; S(2,3) = S1133;
S(3,1) = S3311; S(3,2) = S3311;
S(4,4) = 2*S1313; S(5,5) = 2*S1313;   % gamma_23, gamma_13
S(6,6) = 2*S1212;                     % gamma_12

end

function v = horner(tab, nu, e)
% Evaluate sum_i (a_i nu + b_i)/den_i e^i, divided by (1 - nu)
c = (tab(:,1)*nu + tab(:,2))./tab(:,3);
v = c(end);
for ii = numel(c)-1:-1:1
    v = v*e + c(ii);
end
v = v/(1 - nu);
end
