function dy = cpdem_AnisoDeriv(phi, y, prm)
% d/dphi of the anisotropic DEM state [sigPar; sigPerp; C11; C13; C33; C44; C66]
%
%   dy = cpdem_AnisoDeriv(phi, y, prm)
%
% Thin wrapper around cpdem_AnisoRates (which returns (1 - phi) dy/dphi);
% see that function for the equations.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

dy = cpdem_AnisoRates(y(1), y(2), y(3:7), prm)/(1 - phi);

end
