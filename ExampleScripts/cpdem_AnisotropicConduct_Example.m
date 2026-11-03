%% Example: Anisotropic inclusions and aligned inclusions (conductivity)
%
% Two extensions of the thermal-electrical cross-property model of Cilli
% and Chapman (2021, arXiv:2102.12434), which assumes randomly oriented
% inclusions:
%
%   1. an inclusion material that is itself anisotropic, such as graphite
%      flakes in a copper matrix, whose in-plane and through-plane
%      conductivities differ by orders of magnitude (cpdem_CondIncrement,
%      accepted by cpdem_ElecDEM, cpdem_ForwardConduct and
%      cpdem_ImpliedPorosity as sig2 = [axis plane]);
%
%   2. inclusions that are all aligned, so that the composite is
%      transversely isotropic (cpdem_AlignedCondDEM) and the
%      cross-property relation links the components of the two property
%      tensors (cpdem_AlignedForwardConduct).
%
% The constituent values below are REPRESENTATIVE of a copper-graphite
% composite and are not those of any particular data set.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

%% Constituents (SI)
lamCu  = 390;      % thermal conductivity of copper, W/(m K)
sigCu  = 5.9e7;    % electrical conductivity of copper, S/m
% Graphite flake: [through-plane (flake axis)  in-plane]
lamGr  = [10 400];     % W/(m K)
sigGr  = [1e3 1e5];    % S/m
alpha  = 0.1;          % flake aspect ratio (thickness / diameter)
phi    = linspace(0, 0.5, 26);
phi(1) = 1e-6;

%% 1. Randomly oriented flakes: isotropic composite, anisotropic inclusions
lamIso  = cpdem_ElecDEM(phi, lamCu, lamGr, alpha);          % anisotropic flakes
lamLow  = cpdem_ElecDEM(phi, lamCu, lamGr(1), alpha);       % if flakes had the low value everywhere
lamHigh = cpdem_ElecDEM(phi, lamCu, lamGr(2), alpha);       % ... or the high value
sigIso  = cpdem_ElecDEM(phi, sigCu, sigGr, alpha);

fprintf('Randomly oriented flakes, alpha = %.2f\n', alpha);
fprintf('%6s %12s %12s %12s %14s\n', 'phi', 'lam(aniso)', 'lam(low)', 'lam(high)', 'sigma');
for jj = [2 6 11 16 21 26]
    fprintf('%6.2f %12.1f %12.1f %12.1f %14.3e\n', phi(jj), lamIso(jj), lamLow(jj), lamHigh(jj), sigIso(jj));
end

% Cross-property: thermal from electrical with anisotropic flakes, no
% volume fraction needed
lamFromSig = cpdem_ForwardConduct(lamCu, sigCu, lamGr, sigGr, sigIso(2:end), alpha);
fprintf('\nThermal from electrical (randomly oriented): max relative error vs the\n');
fprintf('volume-fraction DEM %.1e\n', max(abs(lamFromSig - lamIso(2:end))./lamIso(2:end)));

%% 2. Aligned flakes: transversely isotropic composite
[lamPar, lamPerp] = cpdem_AlignedCondDEM(phi, lamCu, lamGr, alpha);
[sigPar, sigPerp] = cpdem_AlignedCondDEM(phi, sigCu, sigGr, alpha);
fprintf('\nAligned flakes (axis = flake normal):\n');
fprintf('%6s %10s %10s %12s %12s\n', 'phi', 'lam_par', 'lam_perp', 'sig_par', 'sig_perp');
for jj = [2 6 11 16 21 26]
    fprintf('%6.2f %10.1f %10.1f %12.3e %12.3e\n', phi(jj), lamPar(jj), lamPerp(jj), sigPar(jj), sigPerp(jj));
end

% Cross-property: measure sigma across the flakes (par) -> both thermal
% components and the in-plane electrical conductivity
o = cpdem_AlignedForwardConduct(lamCu, sigCu, lamGr, sigGr, sigPar(2:end), alpha, 'par');
fprintf('\nAligned cross-property from sigma_par: max relative error in lam_par %.1e, lam_perp %.1e, sigma_perp %.1e\n', ...
    max(abs(o.pPar - lamPar(2:end))./lamPar(2:end)), max(abs(o.pPerp - lamPerp(2:end))./lamPerp(2:end)), ...
    max(abs(o.sigPerp - sigPerp(2:end))./sigPerp(2:end)));

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Anisotropic conductivity');
    subplot(1,2,1)
    plot(phi, lamIso, 'k-', phi, lamPar, 'b-', phi, lamPerp, 'r-', phi, lamLow, 'k:', phi, lamHigh, 'k--', 'LineWidth', 1.5)
    grid on; xlabel('Flake volume fraction'); ylabel('Thermal conductivity (W m^{-1} K^{-1})');
    legend({'random', 'aligned, along axis', 'aligned, in plane', 'random, low value', 'random, high value'}, 'Location', 'SouthWest');
    subplot(1,2,2)
    semilogy(phi, sigIso, 'k-', phi, sigPar, 'b-', phi, sigPerp, 'r-', 'LineWidth', 1.5)
    grid on; xlabel('Flake volume fraction'); ylabel('Electrical conductivity (S/m)');
    legend({'random', 'aligned, along axis', 'aligned, in plane'}, 'Location', 'SouthWest');
    cpdem_FinishFigure(fig, 'AnisotropicConduct');
end
