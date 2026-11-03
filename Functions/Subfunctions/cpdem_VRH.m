function [K, mu, nu, KV, KR, muV, muR] = cpdem_VRH(C)
% Voigt-Reuss-Hill isotropic moduli of a transversely isotropic stiffness
%
%   [K, mu, nu] = cpdem_VRH(C)
%   [K, mu, nu, KV, KR, muV, muR] = cpdem_VRH(C)
%
% Implements eqs 27-32 of Cilli, "Anisotropic cross-property differential
% effective medium modelling": the Voigt and Reuss bulk and shear moduli of
% the 6x6 Voigt stiffness C (engineering shear convention), their Hill
% averages and the isotropic Poisson ratio nu used by the isotropic-host
% Eshelby closure. For an isotropic C the exact moduli are returned.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

s = inv(C);
C11 = C(1,1); C12 = C(1,2); C13 = C(1,3); C33 = C(3,3); C44 = C(4,4); C66 = C(6,6);
KV  = (2*C11 + C33 + 2*C12 + 4*C13)/9;
muV = (2*C11 + C33 - C12 - 2*C13 + 6*C44 + 3*C66)/15;
KR  = 1/(2*s(1,1) + 2*s(1,2) + 4*s(1,3) + s(3,3));
muR = 15/(8*s(1,1) + 4*s(3,3) - 4*s(1,2) - 8*s(1,3) + 6*s(4,4) + 3*s(6,6));
K  = (KV + KR)/2;
mu = (muV + muR)/2;
nu = (3*K - 2*mu)/(2*(3*K + mu));

end
