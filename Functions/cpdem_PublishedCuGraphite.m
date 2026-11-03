function c = cpdem_PublishedCuGraphite()
% Constituent properties and calibration of Cilli & Chapman (arXiv:2102.12434)
%
%   c = cpdem_PublishedCuGraphite()
%
% The thermal-electrical cross-property DEM of Cilli and Chapman (2021,
% arXiv:2102.12434) was tested on the copper-graphite composites of
% Mazloum et al. (2016): five samples of copper containing randomly
% oriented, oblate graphite flakes of measured aspect ratio 0.1 at volume
% fractions 0, 0.1, 0.3, 0.4 and 0.5. The constituent conductivities
% (Section 3 of the paper, taken from Kovacik & Emmer 2011 and Kovacik &
% Bielek 1996) and the parameters solved for in Section 4 are returned in
% the units quoted by the paper: W/(m K) for thermal conductivity and
% "1/(ohm m)" for electrical conductivity as printed there. The
% cross-property relation is invariant to a common scale factor of the
% electrical conductivities, so the unit prefix does not affect the
% inferred aspect ratio or the modelled thermal conductivities.
%
% Fields
%   kappa1, sig1        : copper (isotropic host)      348.6, 58.8e-8
%   kappaFlakeAxis      : graphite along the flake's unique (short) axis   274
%   kappaFlakePlane     : graphite in the plane of the two long axes        10
%   sigFlakeAxis        : likewise, electrical                            0.59e-8
%   sigFlakePlane       :                                                 2.26e-8
%   alphaMeasured       : measured flake aspect ratio                      0.1
%   phiSamples          : graphite volume fractions of the five samples
%   phiCalibration      : fractions used in the calibration (<= 0.4)
%   kappa2, sig2, alpha : effective inclusion conductivities and aspect
%                         ratio solved by the paper: 24.1, 0.61e-8, 0.097
%   alphaCurves         : aspect ratios of the modelled trends in Figure 1
%
% The measured conductivities of the samples themselves are in Mazloum et
% al. (2016) and are not reproduced here.
%
% References:
% Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal
% conductivity through cross-property inclusion modelling. arXiv:2102.12434
%
% Mazloum, A., Kovacik, J., Emmer, S., and Sevostianov, I. (2016).
% Copper-graphite composites: thermal expansion, thermal and electrical
% conductivities, and cross-property connections. Journal of Materials
% Science, 51(17), 7977-7990.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

c = struct( ...
    'kappa1',          348.6, ...
    'sig1',            58.8e-8, ...
    'kappaFlakeAxis',  274, ...
    'kappaFlakePlane', 10, ...
    'sigFlakeAxis',    0.59e-8, ...
    'sigFlakePlane',   2.26e-8, ...
    'alphaMeasured',   0.1, ...
    'phiSamples',      [0 0.1 0.3 0.4 0.5], ...
    'phiCalibration',  [0 0.1 0.3 0.4], ...
    'kappa2',          24.1, ...
    'sig2',            0.61e-8, ...
    'alpha',           0.097, ...
    'alphaCurves',     [0 1e-2 1e-1 1 Inf], ...
    'source',          'Cilli & Chapman (2021), arXiv:2102.12434, Sections 3-4');

end
