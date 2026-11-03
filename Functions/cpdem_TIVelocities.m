function v = cpdem_TIVelocities(C11, C13, C33, C44, C66, rho)
% Phase velocities and Thomsen parameters of a transversely isotropic medium
%
%   v = cpdem_TIVelocities(C11, C13, C33, C44, C66, rho)
%
% Inputs are the five independent stiffnesses (Pa; arrays of a common size
% are allowed) and the density (kg/m^3). The symmetry axis is x3. Output
% struct fields:
%   vp0    : P velocity along the symmetry axis,  sqrt(C33/rho)
%   vp90   : P velocity normal to it,             sqrt(C11/rho)
%   vsv    : SV (and SH along the axis) velocity, sqrt(C44/rho)
%   vsh90  : SH velocity normal to the axis,      sqrt(C66/rho)
%   epsilon, gamma, delta : Thomsen (1986) anisotropy parameters
%
% Reference: Thomsen, L. (1986). Weak elastic anisotropy. Geophysics,
% 51(10), 1954-1966.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

v.vp0   = sqrt(C33./rho);
v.vp90  = sqrt(C11./rho);
v.vsv   = sqrt(C44./rho);
v.vsh90 = sqrt(C66./rho);
v.epsilon = (C11 - C33)./(2*C33);
v.gamma   = (C66 - C44)./(2*C44);
v.delta   = ((C13 + C44).^2 - (C33 - C44).^2)./(2*C33.*(C33 - C44));

end
