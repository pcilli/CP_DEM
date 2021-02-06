%% Example: Model elastic moduli
%
% An example script demonstrating how to forward model elastic moduli from
% resistivity measurements using function "cpdem_Forward".
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
% To estimate effective bulk or shear modulus from conductivity
% measurements, you will need to set the following variables. They can
% either be scalar values (that is, for a single measurement), or
% vectors, for a set of measurements. This has only been tested with
% moduli units Pa and conductivity units 1/(ohm m).
%
% k1:   Background (mineral) bulk modulus 
% mu1:  Background (mineral) shear modulus
% sig1: Background (mineral) conductivity 
%
% k2:   Inclusion (pore fluid) bulk modulus 
% mu2:  Inclusion (pore fluid) shear modulus 
% sig2: Inclusion (pore fluid) conductivity 
%
% sig0: Measured bulk (effective) conductivity 
% epar: Equivalent pore aspect ratio (EPAR)
%
% Note, a different EPAR - bulk or shear modulus EPAR - must be used in
% general to predict either the bulk modulus or shear modulus with
% accuracy, as discussed by Cilli and Chapman (2021).
%
% Cilli and Chapman (2021) calibrate their model to the sandstone data of
% Han et al. (2011) at 60 MPa differential pressure, finding the following
% EPARs:
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
% If you have data with both electrical and elastic measurements available,
% you can invert for the EPARs of individual samples or the best-fitting
% EPARs for a collection of samples using function cpdem_SolveEPAR; see
% script cpdem_SolveEPAR_Example for details on how to do this. EPAR is
% assumed to be known a priori when forward-modelling elastic moduli from
% resistivity measurements in this script.

%% Forward model multiple samples with constant model parameters
% Use the following block of code for a scalar or vector of resistivity
% measurements. If using k-EPAR, only the kModel output is physically
% meaningful. If using mu-EPAR, only the muModel output is physically
% meaningful. The other modulus will probably be inaccurate and unreliable.
% The solver outputs both as they are solved for simultaneously, being the
% two outputs of two coupled differential equations.
% 
% The function cpdem_Forward is set up here for a vector of measurements
% and scalar parameters. To use sample-dependent parameters in a vector
% form, put index "(ii)" after any vector parameters listed in the
% functions argument. Make sure there is one parameter per measurement in
% this case, so the parameter and measurement (sig0) vectors have the same
% length.

kModel  = zeros(size(sig0));
muModel = zeros(size(sig0));

for ii = 1:numel(sig0)
[kTmp,muTmp] = cpdem_Forward(k1,mu1,sig1,k2,mu2,sig2,sig0(ii),epar);
    
    kModel(ii)  = kTmp;
    muModel(ii) = muTmp;
end
