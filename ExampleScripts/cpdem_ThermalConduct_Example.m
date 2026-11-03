%% Example: Thermal-electrical cross-property DEM (Cilli & Chapman, arXiv:2102.12434)
%
% Reproduces the modelling of the paper "Linking electrical and thermal
% conductivity through cross-property inclusion modelling":
%
%   * equation 5, thermal conductivity from electrical conductivity, and
%     its reciprocal equation 6, electrical from thermal
%     (cpdem_ForwardConduct with the properties swapped);
%   * the trends of Figure 1 for inclusion aspect ratios 0 (disks), 1e-2,
%     1e-1, 1 (spheres) and Inf (needles), which all converge at the
%     end members and lie on or within the thermal-electrical
%     Hashin-Shtrikman bounds (cpdem_CondCondHS, Carcione stratagem);
%   * the three-parameter calibration of Section 3: aspect ratio and
%     effective (isotropic) inclusion conductivities, the latter bounded by
%     the graphite flake's axis and in-plane values (cpdem_CalibrateConduct).
%
% The constituent properties and solved parameters are those of the paper
% for the copper-graphite composites of Mazloum et al. (2016)
% (cpdem_PublishedCuGraphite). The measured sample conductivities are in
% Mazloum et al. and are not redistributed, so the calibration below is
% demonstrated on SYNTHETIC measurements generated with the paper's solved
% parameters at the paper's volume fractions.
%
% Reference:
% Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal
% conductivity through cross-property inclusion modelling. arXiv:2102.12434
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

%% The paper's parameters
c = cpdem_PublishedCuGraphite();
fprintf('Copper: kappa1 = %.1f W/(m K), sigma1 = %.3g (paper units)\n', c.kappa1, c.sig1);
fprintf('Graphite effective inclusion values solved by the paper: kappa2 = %.1f, sigma2 = %.3g, alpha = %.3f\n', ...
    c.kappa2, c.sig2, c.alpha);

%% Figure 1a: thermal from electrical (equation 5) for several aspect ratios
sig = linspace(c.sig2*1.001, c.sig1, 200);
kap = zeros(numel(c.alphaCurves), numel(sig));
for ii = 1:numel(c.alphaCurves)
    kap(ii, :) = cpdem_ForwardConduct(c.kappa1, c.sig1, c.kappa2, c.sig2, sig, c.alphaCurves(ii));
end
[kUp, kLow] = cpdem_CondCondHS(sig, c.sig1, c.sig2, c.kappa1, c.kappa2);

%% Figure 1b: electrical from thermal (equation 6), the mirror image
kapAxis = linspace(c.kappa2*1.001, c.kappa1, 200);
sigFromKap = zeros(numel(c.alphaCurves), numel(kapAxis));
for ii = 1:numel(c.alphaCurves)
    sigFromKap(ii, :) = cpdem_ForwardConduct(c.sig1, c.kappa1, c.sig2, c.kappa2, kapAxis, c.alphaCurves(ii));
end
[sUp, sLow] = cpdem_CondCondHS(kapAxis, c.kappa1, c.kappa2, c.sig1, c.sig2);

% The two directions are inverses of each other
chk = cpdem_ForwardConduct(c.sig1, c.kappa1, c.sig2, c.kappa2, kap(3, 2:end), c.alphaCurves(3));
fprintf('Equations 5 and 6 are reciprocal: max relative round-trip error %.1e\n', max(abs(chk - sig(2:end))./sig(2:end)));

%% Section 3: three-parameter calibration on synthetic data
% "Measurements" at the paper's calibration fractions (0.1, 0.3, 0.4),
% generated with the solved parameters. The paper solved for alpha, kappa2
% and sigma2 simultaneously by minimising the l2 misfit in thermal
% conductivity, with the inclusion values bounded by the flake's principal
% values and alpha searched over [0, Inf]. cpdem_CalibrateConduct profiles
% the misfit over aspect ratio and recovers the generating parameters.
phi = c.phiCalibration(c.phiCalibration > 0);
sigMeas = cpdem_ElecDEM(phi, c.sig1, c.sig2, c.alpha);
kapMeas = cpdem_ElecDEM(phi, c.kappa1, c.kappa2, c.alpha);
tic
cal = cpdem_CalibrateConduct(kapMeas, sigMeas, c.kappa1, c.sig1, ...
    [c.kappaFlakePlane c.kappaFlakeAxis], [c.sigFlakeAxis c.sigFlakePlane], 'Refine', 1);
t = toc;
fprintf('\nCalibration on noise-free synthetic data (%.1f s):\n', t);
fprintf('  alpha  = %.4f   (paper: %.3f; measured flakes %.1f)\n', cal.epar, c.alpha, c.alphaMeasured);
fprintf('  kappa2 = %.2f    (paper: %.1f; bounds [%g %g])\n', cal.p2, c.kappa2, c.kappaFlakePlane, c.kappaFlakeAxis);
fprintf('  sigma2 = %.3g (paper: %.3g; bounds [%g %g])\n', cal.sig2, c.sig2, c.sigFlakeAxis, c.sigFlakePlane);
fprintf('  Profile of the misfit over aspect ratio (inclusion values optimised at each):\n');
fprintf('    alpha : %s\n', mat2str(cal.profile.alpha, 3));
fprintf('    misfit: %s\n', mat2str(cal.profile.misfit, 3));

%% ... and the effect of measurement noise
% With three informative samples and three parameters, small noise makes
% the profile nearly flat over a range of aspect ratios: the aspect ratio
% trades off against the effective inclusion values and the three are then
% poorly resolved. This is worth bearing in mind when interpreting
% calibrations of real data sets of this size.
sigNoisy = sigMeas.*(1 + 0.005*[1 -1 0.5]);
kapNoisy = kapMeas.*(1 + 0.005*[-0.5 1 -1]);
calN = cpdem_CalibrateConduct(kapNoisy, sigNoisy, c.kappa1, c.sig1, ...
    [c.kappaFlakePlane c.kappaFlakeAxis], [c.sigFlakeAxis c.sigFlakePlane], ...
    'GridSize', [9 3 3], 'Refine', 1, 'TolX', 1e-3, 'Limits', false);
fprintf('\nWith 0.5 %% noise: alpha = %.3f, kappa2 = %.1f, sigma2 = %.3g, R^2 = %.4f\n', ...
    calN.epar, calN.p2, calN.sig2, calN.r2);
fprintf('  profile misfit over alpha %s:\n  %s\n', mat2str(calN.profile.alpha, 2), mat2str(calN.profile.misfit, 3));

%% Plot (Figure 1 of the paper)
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Thermal-electrical cross-property DEM');
    labels = {'\alpha = 0', '\alpha = 10^{-2}', '\alpha = 10^{-1}', '\alpha = 1', '\alpha = \infty'};
    subplot(1,2,1); hold on
    plot(sig, kap', 'LineWidth', 1.5); plot(sig, kUp, 'k--', sig, kLow, 'k--');
    plot(c.sig1, c.kappa1, 'ks', c.sig2, c.kappa2, 'kd', 'MarkerFaceColor', 'k');
    plot(sigMeas, kapMeas, 'ko', 'MarkerSize', 8);
    grid on; xlabel('Electrical conductivity (paper units)'); ylabel('Thermal conductivity (W m^{-1} K^{-1})');
    legend([labels {'HS bounds', '100 % copper', '100 % graphite (effective)', 'synthetic data'}], 'Location', 'SouthEast');
    title('a) Thermal from electrical (eq. 5)');
    subplot(1,2,2); hold on
    plot(kapAxis, sigFromKap', 'LineWidth', 1.5); plot(kapAxis, sUp, 'k--', kapAxis, sLow, 'k--');
    plot(c.kappa1, c.sig1, 'ks', c.kappa2, c.sig2, 'kd', 'MarkerFaceColor', 'k');
    grid on; xlabel('Thermal conductivity (W m^{-1} K^{-1})'); ylabel('Electrical conductivity (paper units)');
    title('b) Electrical from thermal (eq. 6)');
    cpdem_FinishFigure(fig, 'ThermalConduct');
end
