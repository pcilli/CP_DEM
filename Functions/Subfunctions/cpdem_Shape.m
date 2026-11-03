function [theta, f, L] = cpdem_Shape(epar)
% Spheroid shape functions theta, f (Berryman 1980) and depolarisation L
%
%   [theta, f, L] = cpdem_Shape(epar)
%
% Returns the geometric functions of a spheroidal inclusion with aspect
% ratio EPAR that are needed by the elastic (Berryman 1980) and electrical
% (Osborn 1945; Berryman 1995; Cilli & Chapman 2021) inclusion models:
%
%   theta, f : Berryman's (1980) theta and f, used in P and Q
%   L        : depolarisation factor along the unique (symmetry) axis of the
%              spheroid (Osborn 1945). The two degenerate axes have
%              depolarisation factor (1 - L)/2. Note L == 1 - theta.
%
% Definition of EPAR (Cilli & Chapman 2021): the ratio of the unique
% ellipsoidal axis length to the degenerate axis length, so EPAR < 1 is an
% oblate spheroid ("disk"/"penny"), EPAR == 1 is a sphere and EPAR > 1 is a
% prolate spheroid ("needle"). Berryman (1980) and Osborn (1945) use other
% conventions (see the Theory notes); the conversions are made here.
%
% The closed-form expressions suffer catastrophic cancellation near the
% sphere (EPAR -> 1), where they are undefined (0/0) and where f in
% particular loses accuracy (about 1e-8 relative at |EPAR - 1| = 1e-3 in
% double precision). For |EPAR - 1| < 0.1 a 13th-order Taylor series about
% the sphere is therefore used, with e = EPAR - 1,
%
%   theta = 2/3 + 4e/15 - 6e^2/35 + 32e^3/315 - 40e^4/693 + ...
%   f     = -2/5 - 12e/35 + 2e^2/15 - 32e^3/1155 - 40e^4/3003 + ...
%   L     = 1 - theta,
%
% whose truncation error is below 1e-16 over that range; outside it the
% closed forms are accurate to better than 1e-13.
%
% References:
% Berryman, J. G. (1980). Long-wavelength propagation in composite elastic
% media II. Ellipsoidal inclusions. The Journal of the Acoustical Society
% of America, 68(6), 1820-1831, DOI:10.1121/1.385172
%
% Osborn, J. A. (1945). Demagnetizing factors of the general ellipsoid.
% Physical Review, 67(11-12), 351, DOI:10.1103/PhysRev.67.351
%
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if ~isscalar(epar) || ~isreal(epar) || ~isfinite(epar) || epar <= 0
    error('CPDEM:badEPAR', ...
        'EPAR must be a positive, finite, real scalar (got %s).', mat2str(epar));
end

e = epar - 1;

if abs(e) < 0.1
    % Taylor series about the sphere (see help text above), evaluated by
    % Horner's rule. Exact rational coefficients:
    % theta: 2/3, 4/15, -6/35, 32/315, -40/693, 32/1001, -112/6435,
    %        1024/109395, -1152/230945, 2560/969969, -2816/2028117,
    %        12288/16900975, -13312/35102025, 28672/145422675
    % f:     -2/5, -12/35, 2/15, -32/1155, -40/3003, 32/1365, -784/36465,
    %        1024/62985, -18048/1616615, 7680/1062347, -76032/16900975,
    %        413696/152108775, -545792/339319575, 1404928/1502700975
    ct = [2/3, 4/15, -6/35, 32/315, -40/693, 32/1001, -112/6435, ...
          1024/109395, -1152/230945, 2560/969969, -2816/2028117, ...
          12288/16900975, -13312/35102025, 28672/145422675];
    cf = [-2/5, -12/35, 2/15, -32/1155, -40/3003, 32/1365, -784/36465, ...
          1024/62985, -18048/1616615, 7680/1062347, -76032/16900975, ...
          413696/152108775, -545792/339319575, 1404928/1502700975];
    theta = ct(end);
    f     = cf(end);
    for ii = numel(ct)-1:-1:1
        theta = theta*e + ct(ii);
        f     = f*e + cf(ii);
    end
    L = 1 - theta;

elseif epar < 1
    % Oblate spheroid: Berryman's a = 1 (degenerate axes), c = EPAR.
    a = 1;
    c = epar;
    theta = (a^2*c)/(a^2 - c^2)^(3/2)*(acos(c/a) - (c/a)*sqrt(1 - c^2/a^2));
    f     = c^2/(a^2 - c^2)*(3*theta - 2);

    % Osborn (1945) oblate spheroid with M = (long axis)/(short axis) > 1.
    M = 1/epar;
    L = M^2/(M^2 - 1)*(1 - asin(sqrt(M^2 - 1)/M)/sqrt(M^2 - 1));

else
    % Prolate spheroid: Berryman's a = EPAR (unique axis), c = 1.
    a = epar;
    c = 1;
    theta = (a*c^2)/(a^2 - c^2)^(3/2)*((a/c)*sqrt(a^2/c^2 - 1) - acosh(a/c));
    f     = a^2/(a^2 - c^2)*(2 - 3*theta);

    % Osborn (1945) prolate spheroid with M = (long axis)/(short axis) > 1.
    M = epar;
    s = sqrt(M^2 - 1);
    L = 1/(M^2 - 1)*(M/(2*s)*log((M + s)/(M - s)) - 1);
end

end
