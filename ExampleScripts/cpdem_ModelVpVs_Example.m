%% Example: Model Vp/Vs ratios
%
% An example script demonstrating how to estimate Vp/Vs ratios from
% electrical measurements using function "cpdem_Forward".
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, DOI:10.1093/gji/ggab046
% 
% Han, T., Best, A.I., Sothcott, J. and MacGregor, L.M. (2011), Joint
% elastic?electrical properties of reservoir sandstones and their
% relationships with petrophysical parameters. Geophysical Prospecting, 59:
% 518-535. DOI:10.1111/j.1365-2478.2010.00940.x
%
% Gardner, G. H. F., Gardner, L. W., & Gregory, A. R. (1974). Formation
% velocity and density—The diagnostic basics for stratigraphic traps.
% Geophysics, 39(6), 770-780.
% 
% Written by Phil Cilli, January 2021 as a part of CP_DEM Version 1.0

close all
clear all

%Set random number generator seed for general reproduceability
rng(1,'twister');

%% Preamble
% Here we need to model both effective bulk modulus using K-EPAR then model
% shear modulus using mu-EPAR, so both EPARs need to be known before
% modelling (Cilli and Chapman (2021) show they seem to be different in
% general for real rocks as real rocks do not honour the assumptions of
% theoretical inclusion modelling).

% See script "cpdem_ModelModuli_Template" for more details on how to use
% function "cpdem_Forward".

%% Model moduli
% Make sure both K-EPAR and mu-EPAR are set for this script to run. Cilli
% and Chapman (2021) calibrate their model to the sandstone data of Han et
% al. (2011) at 60 MPa differential pressure, finding the following EPARs:
%
% Mixed sandstones:
% K-EPAR  = 16.4
% mu-EPAR = 12.8
%
% Clean sandstones:
% K-EPAR  = 16.6
% mu-EPAR = 12.3
%
% Clay-bearing sandstones:
% K-EPAR  = 16.3
% mu-EPAR = 12.9
%
% The function cpdem_Forward is set up here for a vector of measurements
% and scalar parameters. To use sample-dependent parameters in a vector
% form, put index "(ii)" after any vector parameters listed in the
% functions argument. Make sure there is one parameter per measurement in
% this case, so the parameter and measurement (sig0) vectors have the same
% length.

kModel  = zeros(size(sig0));
muModel = zeros(size(sig0));

for ii = 1 : numel(sig0)
    
    %Forward model bulk modulus with K-EPAR
    [kTmp,muTmp] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0(ii),kepar);
    kModel(ii)  = kTmp;
    
    %Forward model shear modulus with mu-EPAR
    [kTmp,muTmp] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0(ii),muepar);
    muModel(ii)  = muTmp;
    
end

%% Convert to velocities and calculate vp/vs
% Convert moduli to velocities. Here we use the Gardner et al. (1974)
% relation which is only appropriate for sandstones and is useful when
% density is not available. If reliable density measurements are available,
% it is probably better to use them and the equations of linear elasticity
% to estimate velocities.

%Gardner relation parameters
alpha = 0.31*1000; %scale by 1000 for a density in kg/m^3
beta  = 0.25;

%Gardner relations
vp = ((kModel+4./3.*muModel)./alpha).^(1/(2+beta));
vs = sqrt(muModel./(alpha.* vp.^beta));

VpVsRatio  = vp./vs;
