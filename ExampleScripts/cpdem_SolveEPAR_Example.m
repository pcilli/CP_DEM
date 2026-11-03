%% Example: Calibrate the cross-property EPARs to a joint data set
%
% Solves for the bulk modulus EPAR (alpha_K) and the shear modulus EPAR
% (alpha_mu) that best fit a set of joint elastic-electrical measurements,
% with cpdem_CalibrateEPAR, then calibrates every sample individually. The
% data here are SYNTHETIC (cpdem_SyntheticData: generated with the
% published EPARs and 2 % noise) so that the script runs out of the box;
% replace them with your own k0, mu0 and sig0.
%
% The v1.0 interface, cpdem_SolveEPAR, is also shown.
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

%% Data
p = cpdem_ExampleParameters();
d = cpdem_SyntheticData(p);      % d.k0, d.mu0 (Pa), d.sig0 (S/m)

%% Calibrate both EPARs to the whole data set
% The search is over opts.Bounds (default [1 100], i.e. prolate pores, in
% log space). Cilli and Chapman (2021) found prolate EPARs fit sandstones
% better than oblate ones, and the misfit is well behaved over prolate
% EPARs; extend the bounds below 1 only with care, as several local minima
% may then exist.
tic
cal = cpdem_CalibrateEPAR(d.k0, d.mu0, d.sig0, p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2);
t = toc;
fprintf('Calibrated to %d synthetic samples in %.2f s:\n', numel(d.sig0), t);
fprintf('  alpha_K  = %6.2f  (generated with %.1f)  RMS %.2f GPa, MAPE %.2f %%, R^2 %.4f\n', ...
    cal.eparK, p.eparK, cal.rmsK/1e9, cal.mapeK, cal.r2K);
fprintf('  alpha_mu = %6.2f  (generated with %.1f)  RMS %.2f GPa, MAPE %.2f %%, R^2 %.4f\n', ...
    cal.eparMu, p.eparMu, cal.rmsMu/1e9, cal.mapeMu, cal.r2Mu);

%% Options: l1 norm, relative residuals, custom bounds, weights
cal1 = cpdem_CalibrateEPAR(d.k0, d.mu0, d.sig0, p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, ...
    cpdem_Options('Norm', 1, 'Relative', true, 'Bounds', [2 50]));
fprintf('  l1 / relative residuals: alpha_K = %.2f, alpha_mu = %.2f\n', cal1.eparK, cal1.eparMu);

%% Per-sample EPARs
% Each sample gets its own EPAR; the scatter reflects the measurement noise.
% (A subset of samples is used here to keep the example quick: each sample
% needs its own bounded search.)
sub = 1:3:numel(d.sig0);
cps = cpdem_CalibrateEPAR(d.k0(sub), d.mu0(sub), d.sig0(sub), p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, ...
    struct('PerSample', true, 'Bounds', [2 50]));
fprintf('  per-sample alpha_K : mean %.2f, std %.2f (%d samples)\n', mean(cps.eparK), std(cps.eparK), numel(sub));
fprintf('  per-sample alpha_mu: mean %.2f, std %.2f\n', mean(cps.eparMu), std(cps.eparMu));

%% The v1.0 interface still works
% modFlag 0/1 = bulk/shear EPAR, normFlag 1/2 = l1/l2, paramFlag 0/1 =
% scalar/per-sample parameters, then the search bounds.
eparK_v1 = cpdem_SolveEPAR(d.k0, d.mu0, d.sig0, p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, 0, 2, 0, 2, 50);
fprintf('  cpdem_SolveEPAR (v1.0 call): alpha_K = %.2f\n', eparK_v1);

%% Misfit as a function of EPAR
% Illustrates the shape of the objective over prolate aspect ratios.
epars = logspace(0, 2, 40);
misfitK = zeros(size(epars));
for ii = 1:numel(epars)
    kTry = cpdem_Forward(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,d.sig0,epars(ii));
    misfitK(ii) = sqrt(mean((d.k0 - kTry).^2))/1e9;
end

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'EPAR calibration');
    subplot(1,3,1)
    semilogx(d.res0, d.k0/1e9, 'o', d.res0, cal.kModel/1e9, '-', 'LineWidth', 1.5); grid on
    xlabel('Resistivity (\Omega m)'); ylabel('Bulk modulus (GPa)');
    legend({'synthetic data', sprintf('model, \\alpha_K = %.1f', cal.eparK)});
    subplot(1,3,2)
    semilogx(d.res0, d.mu0/1e9, 'o', d.res0, cal.muModel/1e9, '-', 'LineWidth', 1.5); grid on
    xlabel('Resistivity (\Omega m)'); ylabel('Shear modulus (GPa)');
    legend({'synthetic data', sprintf('model, \\alpha_\\mu = %.1f', cal.eparMu)});
    subplot(1,3,3)
    semilogx(epars, misfitK, 'LineWidth', 1.5); hold on
    plot(cal.eparK, cal.rmsK/1e9, 'rp', 'MarkerSize', 12); grid on
    xlabel('\alpha_K'); ylabel('RMS misfit in K (GPa)'); title('Objective function');
    cpdem_FinishFigure(fig, 'SolveEPAR');
end
