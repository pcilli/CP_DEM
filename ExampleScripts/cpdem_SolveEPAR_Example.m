%% Example: Solve for cross-property EPAR
%
% An example script demonstrating how to solve for cross-property EPAR
% using function "cpdem_SolveEPAR".
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
% Written by Phil Cilli, January 2021 as a part of CP_DEM Version 1.0

close all
clear all

%Set random number generator seed for general reproduceability
rng(1,'twister'); 

%% Preamble
% To solve for cross-property DEM EPAR, you will need to set the following
% variables. They can either be scalar values (that is, for a single
% measurement), or vectors, for a set of measurements. This has only been
% tested with moduli units Pa and conductivity units 1/(ohm m).
%
% k1:   Measured (effective) bulk modulus 
% mu1:  Measured (effective) shear modulus
% sig0: Measured (effective) conductivity 
%
% k1:   Background (mineral) bulk modulus 
% mu1:  Background (mineral) shear modulus
% sig1: Background (mineral) conductivity 
%
% k2:   Inclusion (pore fluid) bulk modulus 
% mu2:  Inclusion (pore fluid) shear modulus 
% sig2: Inclusion (pore fluid) conductivity 
%
% Note, a different EPAR - bulk or shear modulus EPAR - must be solved for
% in general to predict either the bulk modulus or shear modulus with
% accuracy, as discussed by Cilli and Chapman (2021). For context, Cilli
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

%% Set options
% modFlag = 0 -> solve for the K-EPAR for forward modelling K
% modFlag = 1 -> solve for the mu-EPAR for forward modelling mu
modFlag = 0;

%normFlag = 1 -> minimise l1 norm of residuals
%normFlag = 2 -> minimise l2 norm of residuals
normFlag = 2;

% paramFlag = 0 -> All parameters are scalars (the same for a whole
% dataset), but measurements k0, mu0, sig0 can still be vectors, all of
% equal length.
% paramFlag = 1 -> All parameters are sample dependent, so
% all k's, mu's, and sig's going into the solver are vectors of the same
% length.
paramFlag = 0;
%% Solve for EPAR
%Solver's search bounds for EPAR. Cilli and Chapman (2021) found the model
%tends to work better for EPAR > 1 (prolate spheroids) but there's no
%reason why the search cant be extended to include oblate spheroids (0 <=
%EPAR < 1) too. The optimisation seems to be convex over all prolate
%spheroids but may have multiple local solutions over all spheroids.
lb = 2;
ub = 50;

%Optimisation routine only tested on MATLAB R2019B.
epar = cpdem_SolveEPAR(k0,mu0,sig0,k1,mu1,sig1,k2,mu2,sig2,modFlag,normFlag,paramFlag,lb,ub);

