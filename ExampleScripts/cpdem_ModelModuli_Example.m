%% Example: Model elastic moduli from electrical conductivity
%
% Forward models bulk and shear moduli from resistivity with the
% cross-property DEM (cpdem_Forward), compares the three published EPAR
% calibrations, shows the porosity the model implies, checks the
% predictions against the electrical-elastic Hashin-Shtrikman bounds
% (cpdem_ElecElasHS) and demonstrates sample-dependent parameters.
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Carcione, J. M., Ursin, B., & Nordskag, J. I. (2007). Cross-property
% relations between electrical conductivity and the seismic velocity of
% rocks. Geophysics, 72(5), E193-E204.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

%% Set up
% Run cpdem_Setup once (or add the Functions folder to the path).
p = cpdem_ExampleParameters();

% A sweep of resistivities (ohm m) and the equivalent conductivities (S/m)
res0 = logspace(log10(1), log10(200), 60);
sig0 = 1./res0;

%% Forward model with each published calibration
% Remember: with the bulk modulus EPAR only the bulk modulus is physically
% meaningful, and with the shear modulus EPAR only the shear modulus is.
% The whole sweep is computed in one ODE integration per call.
tab = cpdem_PublishedEPAR();
kModel  = zeros(numel(tab), numel(sig0));
muModel = zeros(numel(tab), numel(sig0));
phiK    = zeros(numel(tab), numel(sig0));
for ii = 1:numel(tab)
    [kModel(ii,:), ~, phiK(ii,:)] = cpdem_Forward(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,tab(ii).eparK);
    [~, muModel(ii,:)]            = cpdem_Forward(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,tab(ii).eparMu);
end

fprintf('Bulk and shear moduli (GPa) at selected resistivities, "%s" calibration:\n', tab(1).description);
fprintf('%10s %10s %10s %10s\n', 'rho (ohm m)', 'K', 'mu', 'phi_K');
for jj = [1 15 30 45 60]
    fprintf('%10.1f %10.2f %10.2f %10.3f\n', res0(jj), kModel(1,jj)/1e9, muModel(1,jj)/1e9, phiK(1,jj));
end

%% Hashin-Shtrikman bounds
% Every cross-property DEM prediction lies within the bounds of Carcione et
% al. (2007), which only use the constituent properties.
[kUp, kLow, muUp, muLow] = cpdem_ElecElasHS(sig0, p.sig1, p.sig2, p.k1, p.k2, p.mu1);
assert(all(kModel(:) <= max(kUp)) && all(muModel(:) <= max(muUp)), 'bounds violated');

%% Sample-dependent parameters
% Any of k1, mu1, sig1, k2, mu2, sig2 (and EPAR) may be given per sample,
% e.g. a brine conductivity that differs between samples.
sig2PerSample = p.sig2*(1 + 0.2*sin(linspace(0, pi, numel(sig0))));
kVar = cpdem_Forward(p.k1,p.mu1,p.sig1,p.k2,p.mu2,sig2PerSample,sig0,p.eparK);
fprintf('\nMax change in K from a +-20%% brine conductivity variation: %.1f %%\n', ...
    100*max(abs(kVar - kModel(1,:))./kModel(1,:)));

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Moduli from resistivity');
    subplot(1,2,1); hold on
    for ii = 1:numel(tab)
        semilogx(res0, kModel(ii,:)/1e9, 'LineWidth', 1.5);
    end
    semilogx(res0, kUp/1e9, 'k--', res0, kLow/1e9, 'k--');
    set(gca, 'XScale', 'log'); grid on
    xlabel('Resistivity (\Omega m)'); ylabel('Bulk modulus (GPa)');
    legend([{tab.description} {'HS bounds'}], 'Location', 'SouthEast');
    title('Bulk modulus with \alpha_K');
    subplot(1,2,2); hold on
    for ii = 1:numel(tab)
        semilogx(res0, muModel(ii,:)/1e9, 'LineWidth', 1.5);
    end
    semilogx(res0, muUp/1e9, 'k--', res0, muLow/1e9, 'k--');
    set(gca, 'XScale', 'log'); grid on
    xlabel('Resistivity (\Omega m)'); ylabel('Shear modulus (GPa)');
    title('Shear modulus with \alpha_\mu');
    cpdem_FinishFigure(fig, 'ModelModuli');
end
