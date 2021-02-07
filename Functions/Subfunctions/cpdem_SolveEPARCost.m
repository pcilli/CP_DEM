function cost = cpdem_SolveEPARCost(x0,k0,mu0,sig0,k1,mu1,sig1,k2,mu2,sig2,modFlag,normFlag,paramFlag)
% Cost function for cpdem_SolveEPAR
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Written by Phil Cilli, January 2021 as a part of Cross-Property DEM
% Toolbox Version 1.0

clearvars -except k0 mu0 sig0 k1 mu1 sig1 k2 mu2 sig2 modFlag normFlag paramFlag x0;

epar = x0;

kModel  = zeros(size(sig0));
muModel = zeros(size(sig0));

if paramFlag == 0 % If all parameters are scalar:
    for ii = 1 : length(sig0)
        [kTmp,muTmp] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0(ii),epar);
        kModel(ii)  = kTmp;
        muModel(ii) = muTmp;
    end
    
elseif paramFlag == 1 % If all parameters are sample dependent:
    for ii = 1 : length(sig0)
        [kTmp,muTmp] = cpdem_Forward(k1(ii),mu1(ii),sig1(ii),k2(ii),mu2(ii),sig2(ii),sig0(ii),epar);
        kModel(ii)  = kTmp;
        muModel(ii) = muTmp;
    end
    
else
    disp('Invalid paramFlag. Set to 0 or 1')
end


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
