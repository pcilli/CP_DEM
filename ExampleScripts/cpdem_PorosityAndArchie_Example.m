%% Example: The porosity behind the model, Archie's law and EPAR(phi)
%
% Porosity is a dummy variable of the cross-property DEM (Cilli & Chapman
% 2021): the model links moduli to conductivity through the two single-
% property DEMs without ever needing it. This script looks "inside" the
% model:
%
%   * the electrical DEM (cpdem_ElecDEM) and elastic DEM (cpdem_ElasDEM)
%     as functions of porosity, and the porosity the model implies for a
%     measured conductivity (cpdem_ImpliedPorosity);
%   * m-bar as a generalisation of Archie's cementation exponent: the
%     formation factor F = sig2/sig follows a power law in porosity whose
%     exponent depends on pore shape (Sen et al. 1981; Mendelson & Cohen
%     1982) and, in general, on porosity itself;
%   * a porosity-dependent aspect ratio, as in the power-law relation of
%     Cilli and Chapman (2020), passed as a function handle.
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Cilli, P. A., and Chapman, M. (2020). The power-law relation between
% inclusion aspect ratio and porosity: Implications for electrical and
% elastic modeling. Journal of Geophysical Research: Solid Earth, 125,
% e2019JB019187, DOI:10.1029/2019JB019187
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

%% Set up
p = cpdem_ExampleParameters();
phi = linspace(0.005, 0.4, 80);

%% Single-property DEM curves for a few pore shapes
epars = [0.1 1 16.4];
labels = {'oblate, \alpha = 0.1', 'sphere, \alpha = 1', 'prolate, \alpha = 16.4'};
sig = zeros(numel(epars), numel(phi));
k   = zeros(numel(epars), numel(phi));
mu  = zeros(numel(epars), numel(phi));
for ii = 1:numel(epars)
    sig(ii,:) = cpdem_ElecDEM(phi, p.sig1, p.sig2, epars(ii));
    [k(ii,:), mu(ii,:)] = cpdem_ElasDEM(phi, p.k1, p.mu1, p.k2, p.mu2, epars(ii));
end

%% Implied porosity for a set of conductivities
sig0 = [0.02 0.05 0.1 0.2 0.4];
phiImplied = cpdem_ImpliedPorosity(sig0, p.sig1, p.sig2, p.eparK);
[~, ~, phiFwd] = cpdem_Forward(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,p.eparK);
fprintf('Implied porosity (alpha_K = %.1f):\n', p.eparK);
fprintf('%8s %8s\n', 'sig0', 'phi');
fprintf('%8.3f %8.4f\n', [sig0; phiImplied]);
fprintf('(cpdem_Forward returns the same porosity as its third output: max diff %.1e)\n', ...
    max(abs(phiImplied - phiFwd)));

%% Archie's law and m-bar
% For the electrical DEM, F = sig2/sig. In Archie's law F = phi^-m. The
% apparent cementation exponent m_app = -log(F)/log(phi) varies with
% porosity for a brine-saturated rock, and tends to m-bar's insulating
% limit (1/3) sum_j 1/(1 - L_j) only when the contrast is extreme.
F = p.sig2./sig;
mApp = -log(F)./log(phi(ones(numel(epars),1),:));
mArchie = zeros(size(epars));
for ii = 1:numel(epars)
    mArchie(ii) = mbar(1, 0, epars(ii));
end
fprintf('\nInsulating-limit cementation exponents (1/3) sum 1/(1 - L_j):\n');
for ii = 1:numel(epars)
    fprintf('  %-22s m = %.3f\n', labels{ii}, mArchie(ii));
end

%% Porosity-dependent aspect ratio
% A floored power law alpha(phi) = max(a phi^b, floor): every DEM path
% starts at phi = 0, so the law must stay positive there.
law = @(ph) max(8*ph.^0.5, 0.02);
sigLaw = cpdem_ElecDEM(phi, p.sig1, p.sig2, law);
[kLaw, muLaw] = cpdem_ElasDEM(phi, p.k1, p.mu1, p.k2, p.mu2, law);
[kX, ~, phiX] = cpdem_Forward(p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2,sig0,law);
fprintf('\nWith alpha(phi) = max(8 phi^0.5, 0.02), K at sig0 = %s S/m is %s GPa\n', ...
    mat2str(sig0), mat2str(kX/1e9, 4));

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Inside the model');
    subplot(2,2,1)
    semilogy(phi, sig', 'LineWidth', 1.5); hold on; semilogy(phi, sigLaw, 'k--', 'LineWidth', 1.5)
    grid on; xlabel('Porosity'); ylabel('Conductivity (S/m)'); legend([labels {'\alpha(\phi)'}], 'Location', 'SouthEast');
    subplot(2,2,2)
    plot(phi, k'/1e9, 'LineWidth', 1.5); hold on; plot(phi, kLaw/1e9, 'k--', 'LineWidth', 1.5)
    grid on; xlabel('Porosity'); ylabel('Bulk modulus (GPa)');
    subplot(2,2,3)
    loglog(phi, F', 'LineWidth', 1.5); grid on
    xlabel('Porosity'); ylabel('Formation factor \sigma_2/\sigma');
    subplot(2,2,4)
    plot(phi, mApp', 'LineWidth', 1.5); hold on
    for ii = 1:numel(epars), plot(phi([1 end]), mArchie(ii)*[1 1], ':'); end
    grid on; xlabel('Porosity'); ylabel('Apparent cementation exponent'); ylim([0 4]);
    cpdem_FinishFigure(fig, 'PorosityAndArchie');
end
