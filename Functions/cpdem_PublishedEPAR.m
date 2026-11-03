function tab = cpdem_PublishedEPAR(which)
% Published cross-property EPAR calibrations of Cilli and Chapman (2021)
%
%   tab = cpdem_PublishedEPAR()          struct array of all calibrations
%   s   = cpdem_PublishedEPAR('clean')   one calibration by name
%
% Cilli and Chapman (2021) calibrated the cross-property DEM to the joint
% elastic-electrical laboratory measurements of Han et al. (2011) on
% brine-saturated reservoir sandstones at 60 MPa differential pressure,
% and found the following bulk modulus (alpha_K) and shear modulus
% (alpha_mu) equivalent pore aspect ratios:
%
%   Data set                         alpha_K   alpha_mu
%   all samples                      16.4      12.8
%   clean (clay < 3.5 %)             16.6      12.3
%   clay-bearing (clay > 3.5 %)      16.3      12.9
%
% All are prolate ("needle-like"). They are a useful starting point when
% no joint data are available to calibrate EPAR for a particular rock
% (see cpdem_CalibrateEPAR), but note that EPAR is an effective, not a
% literal, pore shape, that it depends on the constituent properties
% assumed in the calibration, and that it varies (weakly) with pressure
% (Cilli & Chapman, EGU General Assembly 2022, EGU22-9251).
%
% Fields: name, description, eparK, eparMu, source.
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Han, T., Best, A.I., Sothcott, J. and MacGregor, L.M. (2011), Joint
% elastic-electrical properties of reservoir sandstones and their
% relationships with petrophysical parameters. Geophysical Prospecting,
% 59, 518-535, DOI:10.1111/j.1365-2478.2010.00940.x
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

src = 'Cilli & Chapman (2021); Han et al. (2011) sandstones at 60 MPa';
tab = struct( ...
    'name',        {'all', 'clean', 'clay-bearing'}, ...
    'description', {'All sandstone samples', 'Clean sandstones (clay < 3.5 %)', ...
                    'Clay-bearing sandstones (clay > 3.5 %)'}, ...
    'eparK',       {16.4, 16.6, 16.3}, ...
    'eparMu',      {12.8, 12.3, 12.9}, ...
    'source',      {src, src, src});

if nargin > 0
    hit = strcmpi({tab.name}, which);
    if ~any(hit)
        error('CPDEM:badInput', 'Unknown calibration ''%s''. Choose ''all'', ''clean'' or ''clay-bearing''.', which);
    end
    tab = tab(hit);
end

end
