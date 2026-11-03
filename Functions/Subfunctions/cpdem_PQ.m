function [p, q] = cpdem_PQ(k0, mu0, k2, mu2, epar)
% Berryman's (1980) P and Q polarisation factors for a spheroidal inclusion
%
%   [p, q] = cpdem_PQ(k0, mu0, k2, mu2, epar)
%
% P and Q describe the bulk and shear strain concentration in a spheroidal
% inclusion (moduli k2, mu2, aspect ratio EPAR) embedded in a host with
% moduli k0, mu0. They appear in the elastic DEM of Berryman (1992) and in
% the cross-property DEM equations of Cilli and Chapman (2021). The
% expressions are those of the appendix of Berryman (1980), also tabulated
% in Mavko et al. (2009).
%
% Inputs may be scalars; k0 and mu0 may also be equally sized arrays (for
% vectorised evaluation along a DEM path), in which case P and Q have the
% same size. EPAR follows the Cilli and Chapman (2021) definition (see
% cpdem_Shape). The sphere limit (EPAR == 1) reduces to
%   P = (k0 + 4/3 mu0)/(k2 + 4/3 mu0),
%   Q = (mu0 + zeta)/(mu2 + zeta),  zeta = mu0/6 (9 k0 + 8 mu0)/(k0 + 2 mu0).
%
% References:
% Berryman, J. G. (1980). Long-wavelength propagation in composite elastic
% media II. Ellipsoidal inclusions. The Journal of the Acoustical Society
% of America, 68(6), 1820-1831, DOI:10.1121/1.385172
%
% Berryman, J. G. (1992). Single-scattering approximations for coefficients
% in Biot's equations of poroelasticity. The Journal of the Acoustical
% Society of America, 91(2), 551-571, DOI:10.1121/1.402518
%
% Mavko, G., Mukerji, T., and Dvorkin, J. (2009). The Rock Physics
% Handbook, 2nd ed. Cambridge University Press.
%
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if ~isscalar(epar) || epar <= 0 || isinf(epar)
    error('CPDEM:badEPAR', ...
        'The elastic P and Q need 0 < EPAR < Inf (the disk and needle limits are singular).');
end
[T, f] = cpdem_Shape(epar);

nu0 = (3*k0 - 2*mu0)./(2*(3*k0 + mu0));   % host Poisson's ratio

A = mu2./mu0 - 1;
B = (1/3)*(k2./k0 - mu2./mu0);
RR = (1 - 2*nu0)./(2*(1 - nu0));         % Berryman's R (= 3 mu0/(3 k0 + 4 mu0))

F1 = 1 + A.*((3/2)*(f + T) - RR.*((3/2)*f + (5/2)*T - 4/3));

F2 = 1 + A.*(1 + (3/2)*(f + T) - (RR/2).*(3*f + 5*T)) + B.*(3 - 4*RR) + ...
    (A/2).*(A + 3*B).*(3 - 4*RR).*(f + T - RR.*(f - T + 2*T^2));

F3 = 1 + A.*(1 - (f + (3/2)*T) + RR.*(f + T));

F4 = 1 + (A/4).*(f + 3*T - RR.*(f - T));

F5 = A.*(-f + RR.*(f + T - 4/3)) + B.*T.*(3 - 4*RR);

F6 = 1 + A.*(1 + f - RR.*(f + T)) + B.*(1 - T).*(3 - 4*RR);

F7 = 2 + (A/4).*(3*f + 9*T - RR.*(3*f + 5*T)) + B.*T.*(3 - 4*RR);

F8 = A.*(1 - 2*RR + (f/2).*(RR - 1) + (T/2).*(5*RR - 3)) + B.*(1 - T).*(3 - 4*RR);

F9 = A.*((RR - 1).*f - RR.*T) + B.*T.*(3 - 4*RR);

p = F1./F2;
q = (1/5)*((2./F3) + (1./F4) + (F4.*F5 + F6.*F7 - F8.*F9)./(F2.*F4));

end
