function r = R(sig0,sig2,epar)
% Calculate R
%
% Calculate R from Berryman (1995) used in equations 20 & 21 of Cilli and
% Chapman (2021).
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Osborn, J. A. (1945). Demagnetizing factors of the general ellipsoid.
% Physical review, 67(11-12), 351, DOI:10.1103/PhysRev.67.351
%
% Berryman, J. G. (1995). Mixture theories for rock properties. Rock
% physics and phase relations: A handbook of physical constants, 3,
% 205-228. DOI:10.1029/RF003
%
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

% M and L are from Osborn (1945)
if epar < 1 %Oblate spheroid
    M = 1/epar;
    
    L = M^2/(M^2-1)*...
        (1-1/((M^2-1)^(1/2)) *...
        asin((M^2-1)^(1/2)/M));
    
else %Prolate spheroid
    M = epar;
    
    L = 1/(M^2-1)*...
        (M/(2*(M^2-1)^(1/2)) *...
        log((M+(M^2-1)^(1/2))/(M-(M^2-1)^(1/2)))-1);
end

% m-bar is from Cilli and Chapman (2021)
mb = ((1/3)*sig0*(1/(sig0-L*(sig0-sig2))+...
     4/(sig0+sig2+L*(sig0-sig2))));
 
% R from Berryman (1995)
r =  mb/(3*sig0);

end