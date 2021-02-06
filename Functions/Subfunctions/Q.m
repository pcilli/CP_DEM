function q = Q(k0,mu0,k2,mu2,epar)
% Calculate Q
%
% Calculate Q from Berryman (1982) [DOI:10.1121/1.385172] used in equations
% 18 & 19, or 20 & 21 of Cilli and Chapman (2021).
% 
% Note we have defined EPAR as the ratio of the unique ellipsoidal axis
% length to the degenerate axis length, so it spans 0 to infinity for both
% types of spheroid. However, Berryman (1982) defines it as a/c, where the
% prolate ellipsoid has axis lengths a > b = c and oblate ellipsoid has a =
% b < c, so it spans 0 to 1 for both prolate and oblate. These equations
% are taken from the appendix of Berryman (1982), and so a and c are
% calculated from our definition of EPAR accordingly.
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
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

if epar <= 1 % Oblate spheroid
    a = 1;
    c = epar;
    T = (a^2*c)/(a^2-c^2)^(3/2)*(acos(c/a)-(c/a)*sqrt(1-(c^2/a^2)));
    f = c^2/(a^2-c^2)*(3*T-2);

else % Prolate spheroid
    a = epar;
    c = 1;
    T = (a*c^2)/((a^2-c^2)^(3/2))*((a/c)*sqrt((a^2/c^2)-1)-acosh(a/c));
    f = a^2/(a^2-c^2)*(2-3*T);
end

nu0 = (3*k0-2*mu0)/(2*(3*k0+mu0));

A  = mu2/mu0-1;

B  = (1/3)*(k2/k0-mu2/mu0);

%Use double R here to avoid mix-ups with toolbox subfunction function "R"
%and its output "r"
RR  = (1-2*nu0)/(2*(1-nu0));

F1 = 1+A*((3/2)*(f+T)-RR*((3/2)*f+(5/2)*T-(4/3)));

F2 = 1+A*(1+(3/2)*(f+T)-(RR/2)*(3*f+5*T))+B*(3-4*RR)+...
    (A/2)*(A+3*B)*(3-4*RR)*(f+T-RR*(f-T+2*T^2));

F3 = 1+A*(1-(f+(3/2)*T)+RR*(f+T));

F4 = 1+(A/4)*(f+3*T-RR*(f-T));

F5 = A*(-f+RR*(f+T-4/3))+B*T*(3-4*RR);

F6 = 1+A*(1+f-RR*(f+T))+B*(1-T)*(3-4*RR);

F7 = 2+(A/4)*(3*f+9*T-RR*(3*f+5*T))+B*T*(3-4*RR);

F8 = A*(1-2*RR+(f/2)*(RR-1)+(T/2)*(5*RR-3))+B*(1-T)*(3-4*RR);

F9 = A*((RR-1)*f-RR*T)+B*T*(3-4*RR);

q = (1/5)*((2/F3)+(1/F4)+(F4*F5+F6*F7-F8*F9)/(F2*F4));

end

