%% Example: Model Vp, Vs and Vp/Vs from electrical conductivity
%
% Estimates P- and S-wave velocities and their ratio from resistivity with
% cpdem_ModelVpVs, which models the bulk modulus with alpha_K and the shear
% modulus with alpha_mu (both are needed, as Cilli and Chapman (2021) find
% they differ for real rocks) and converts to velocities with a density.
% Three density options are compared: the porosity implied by the model
% itself, measured densities, and the Gardner et al. (1974) relation used
% in v1.0 of the toolbox.
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Gardner, G. H. F., Gardner, L. W., & Gregory, A. R. (1974). Formation
% velocity and density - The diagnostic basics for stratigraphic traps.
% Geophysics, 39(6), 770-780.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

%% Set up
p = cpdem_ExampleParameters();
res0 = logspace(log10(1), log10(200), 60);
sig0 = 1./res0;

%% Velocities with a density from the model's implied porosity
% rho = (1 - phi) rho1 + phi rho2 with phi the porosity the electrical DEM
% implies for each conductivity (mean of the alpha_K and alpha_mu paths).
dens = struct('method', 'implied', 'rho1', p.rho1, 'rho2', p.rho2);
imp = cpdem_ModelVpVs(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,p.eparK,p.eparMu,dens);

%% ... with the Gardner relation (sandstones only; v1.0 behaviour)
gar = cpdem_ModelVpVs(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,p.eparK,p.eparMu);

%% ... and with "measured" densities (here taken from the implied porosity)
mea = cpdem_ModelVpVs(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,p.eparK,p.eparMu, ...
                      struct('method', 'measured', 'rho', imp.rho));

fprintf('%12s %8s %8s %8s %8s %8s\n', 'rho (ohm m)', 'Vp', 'Vs', 'Vp/Vs', 'phi', 'rho_b');
for jj = [1 15 30 45 60]
    fprintf('%12.1f %8.0f %8.0f %8.3f %8.3f %8.0f\n', res0(jj), imp.vp(jj), imp.vs(jj), ...
        imp.vpvs(jj), (imp.phiK(jj) + imp.phiMu(jj))/2, imp.rho(jj));
end
fprintf('\nVp/Vs does not depend on density; Vp differs by up to %.1f %% between the\n', ...
    100*max(abs(gar.vp - imp.vp)./imp.vp));
fprintf('implied-porosity and Gardner densities over this sweep.\n');

%% Compare the three published calibrations
tab = cpdem_PublishedEPAR();
vpvs = zeros(numel(tab), numel(sig0));
for ii = 1:numel(tab)
    o = cpdem_ModelVpVs(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,tab(ii).eparK,tab(ii).eparMu,dens);
    vpvs(ii,:) = o.vpvs;
end

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Velocities from resistivity');
    subplot(1,2,1)
    semilogx(res0, imp.vp, 'LineWidth', 1.5); hold on
    semilogx(res0, gar.vp, '--', 'LineWidth', 1.5);
    semilogx(res0, imp.vs, 'LineWidth', 1.5);
    semilogx(res0, gar.vs, '--', 'LineWidth', 1.5);
    grid on; xlabel('Resistivity (\Omega m)'); ylabel('Velocity (m/s)');
    legend({'Vp, implied \rho', 'Vp, Gardner \rho', 'Vs, implied \rho', 'Vs, Gardner \rho'}, 'Location', 'SouthEast');
    subplot(1,2,2)
    semilogx(res0, vpvs', 'LineWidth', 1.5);
    grid on; xlabel('Resistivity (\Omega m)'); ylabel('Vp/Vs');
    legend({tab.description}, 'Location', 'NorthEast');
    cpdem_FinishFigure(fig, 'ModelVpVs');
end
