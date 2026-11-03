%% Example: Thermal conductivity from electrical conductivity
%
% The chain-rule construction of the cross-property DEM applies to any two
% properties with a DEM: Cilli and Chapman (2021, arXiv:2102.12434) link
% electrical and thermal conductivity, which obey the same Laplace-type
% equation and therefore share the geometric factor m-bar. The
% cross-property relation is again free of porosity, depending only on the
% end-member properties and the pore aspect ratio (cpdem_ForwardConduct).
%
% References:
% Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal
% conductivity through cross-property inclusion modelling.
% arXiv:2102.12434.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

%% Set up
p = cpdem_ExampleParameters();
lam1 = 7.7;   % thermal conductivity of quartz, W/(m K)
lam2 = 0.6;   % thermal conductivity of water, W/(m K)
res0 = logspace(log10(1), log10(200), 60);
sig0 = 1./res0;

%% Thermal conductivity for several pore shapes
epars = [0.1 1 16.4];
lam = zeros(numel(epars), numel(sig0));
for ii = 1:numel(epars)
    lam(ii,:) = cpdem_ForwardConduct(lam1, p.sig1, lam2, p.sig2, sig0, epars(ii));
end
fprintf('Thermal conductivity (W/(m K)) at selected resistivities:\n');
fprintf('%12s %10s %10s %10s\n', 'rho (ohm m)', 'a=0.1', 'a=1', 'a=16.4');
for jj = [1 15 30 45 60]
    fprintf('%12.1f %10.3f %10.3f %10.3f\n', res0(jj), lam(:,jj));
end

%% The relation is invertible: electrical from thermal
sigBack = cpdem_ForwardConduct(p.sig1, lam1, p.sig2, lam2, lam(3,:), 16.4);
fprintf('\nRound trip thermal -> electrical: max relative error %.1e\n', max(abs(sigBack - sig0)./sig0));

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Thermal from electrical');
    semilogx(res0, lam', 'LineWidth', 1.5); grid on
    xlabel('Resistivity (\Omega m)'); ylabel('Thermal conductivity (W m^{-1} K^{-1})');
    legend({'\alpha = 0.1', '\alpha = 1', '\alpha = 16.4'}, 'Location', 'SouthEast');
    cpdem_FinishFigure(fig, 'ThermalConduct');
end
