function out = cpdem_AnisoElasDEM(phi, k1, mu1, k2, alpha, S, varargin)
% Transversely isotropic elastic DEM of the anisotropic cross-property paper
%
%   out = cpdem_AnisoElasDEM(phi, k1, mu1, k2, alpha)
%   out = cpdem_AnisoElasDEM(phi, k1, mu1, k2, alpha, S, opts)
%
% Integrates the five-equation elastic DEM of Cilli, "Anisotropic
% cross-property differential effective medium modelling" (eqs 20-39) for
% aligned fluid-filled spheroids (bulk modulus k2, zero shear) of aspect
% ratio alpha in an isotropic matrix (k1, mu1), with the isotropic-host
% Eshelby closure (section 2.3.4): at each step the current transversely
% isotropic stiffness is reduced to isotropic moduli by the Voigt-Reuss-Hill
% average (cpdem_VRH), the spheroidal Eshelby tensor is evaluated for that
% host (cpdem_EshelbySpheroid), and
%
%   (1 - phi) dC/dphi = (C2 - C) [I + S_cl C^-1 (C2 - C)]^-1
%
% is integrated for C11, C13, C33, C44, C66, with C12 = C11 - 2 C66.
% S in [0, 1] is the alignment order parameter of eq. 57 (S = 1, default:
% aligned spheroids; S = 0: the spherical-reference branch).
%
% Output struct with fields C11, C12, C13, C33, C44, C66 (Pa, shape of
% phi) and the VRH moduli K, mu along the path. For spheres (alpha = 1)
% the result equals the isotropic DEM of Berryman (cpdem_ElasDEM) exactly.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

if nargin < 6 || isempty(S), S = 1; end
opts = cpdem_AnisoOptions(varargin);
cpdem_CheckPositive('k1', k1, 'mu1', mu1, 'k2', k2, 'phi>=0', phi);
cpdem_AnisoCheck(alpha, S, 'isotropic-host');
if any(phi(:) >= 1), error('CPDEM:badInput', 'Porosity must be < 1.'); end

prm = struct('sig2', 1, 'k2', k2, 'alpha', alpha, 'S', S, 'L', 1/3, 'closure', 'isotropic-host');
lam = k1 - 2*mu1/3;
c0 = [lam + 2*mu1; lam; lam + 2*mu1; mu1; mu1]/k1;
fun = @(p, y) elasOnly(p, y, prm, k1);
Y = cpdem_Integrate(fun, phi(:), c0, opts)*k1;

shape = size(phi);
out.C11 = reshape(Y(:,1), shape); out.C13 = reshape(Y(:,2), shape);
out.C33 = reshape(Y(:,3), shape); out.C44 = reshape(Y(:,4), shape);
out.C66 = reshape(Y(:,5), shape); out.C12 = out.C11 - 2*out.C66;
out.K = zeros(shape); out.mu = zeros(shape);
for ii = 1:numel(phi)
    [out.K(ii), out.mu(ii)] = cpdem_VRH(cpdem_TIStiffness(Y(ii, :)));
end

end

function dy = elasOnly(phi, y, prm, k1)
r = cpdem_AnisoRates(1, 1, y*k1, prm);
dy = r(3:7)/(k1*(1 - phi));
end
