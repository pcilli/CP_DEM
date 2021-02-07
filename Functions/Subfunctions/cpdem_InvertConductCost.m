function cost = cpdem_InvertConductCost(x0,k0,mu0,k1,mu1,sig1,k2,mu2,sig2,epar,modFlag,normFlag)
% Cost function for cpdem_InvertConduct
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

clearvars -except k0 mu0 k1 mu1 sig1 k2 mu2 sig2 epar modFlag normFlag paramFlag x0;

sig0 = x0;

[kModel,muModel] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0,epar);

if normFlag == 1 %l1 norm
    if modFlag == 0
        cost = sum(abs(k0  - kModel));
    elseif modFlag == 1
        cost = sum(abs(mu0 - muModel));
    else
        disp('Invalid modFlag. Set to 0 or 1')
    end
    
elseif normFlag == 2 %l2 norm
    if modFlag == 0
        cost = norm(k0  - kModel);
    elseif modFlag == 1
        cost = norm(mu0 - muModel);
    else
        disp('Invalid modFlag. Set to 0 or 1')
    end
    
else
    disp('Invalid normFlag. Set to 1 or 2')
end

end
