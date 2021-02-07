function sigInv = cpdem_InvertConduct(k0,mu0,k1,mu1,sig1,k2,mu2,sig2,epar,lb,ub,modFlag,normFlag,paramFlag)
% Invert for electrical conductivity from elastic measurements
%
% Numerically invert equations 18 & 19, or 20 & 21 of Cilli and
% Chapman (2021) for conductivity.
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

if paramFlag == 0
    for ii = 1 : length(k0)
        ObjFn = @(x) cpdem_InvertConductCost(x,k0(ii),mu0(ii),k1,mu1,sig1,k2,mu2,sig2,epar,modFlag,normFlag);
        x = fminbnd(ObjFn,lb,ub);
        sigInv(ii)  = x;      
    end
    
elseif paramFlag == 1
    for ii = 1 : length(k0)
        ObjFn = @(x) cpdem_InvertConductCost(x,k0(ii),mu0(ii),k1(ii),mu1(ii),sig1(ii),k2(ii),mu2(ii),sig2(ii),epar(ii),modFlag,normFlag);
        x = fminbnd(ObjFn,lb,ub);
        sigInv(ii)  = x;
    end
    
else
    disp('Invalid paramFlag. Set to 0 or 1')
end

end
