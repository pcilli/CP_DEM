%% Example: Invert for electrical conductivity from elastic moduli
%
% Inverts measured bulk or shear moduli for conductivity with the
% cross-property DEM. Since v2.0 this is done directly, by integrating the
% cross-property equations with the modulus as the independent variable
% (cpdem_Inverse), which is exact to solver tolerance and fast. The v1.0
% bounded-search interface (cpdem_InvertConduct) is retained.
%
% As discussed by Cilli and Chapman (2021), a different EPAR is needed for
% each modulus: sigma_K is obtained from the bulk modulus with alpha_K, and
% sigma_mu from the shear modulus with alpha_mu. They found sigma_mu to be
% the better estimate of the true conductivity for the clay-bearing
% sandstones of Han et al. (2011), where the shear modulus-conductivity
% cross-plot scatters less than the bulk modulus one; for the few clean
% sandstones sigma_K was slightly better.
%
% The data are SYNTHETIC (see cpdem_SyntheticData).
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
d = cpdem_SyntheticData(p);

%% Direct inversion
tic
[sigK, phiK]   = cpdem_Inverse(d.k0,  'K',  p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, p.eparK);
[sigMu, phiMu] = cpdem_Inverse(d.mu0, 'mu', p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, p.eparMu);
t = toc;

errK  = 100*abs(sigK  - d.sig0)./d.sig0;
errMu = 100*abs(sigMu - d.sig0)./d.sig0;
fprintf('Inverted %d samples twice in %.2f s.\n', numel(d.sig0), t);
fprintf('Mean absolute error in conductivity from %.0f %% noise in the moduli:\n', 100*d.noiseLevel);
fprintf('  sigma_K : %.1f %%   sigma_mu: %.1f %%\n', mean(errK), mean(errMu));
fprintf('(Conductivity is sensitive to the moduli: a small modulus error maps to a\n');
fprintf(' larger conductivity error, more so at low conductivity where the moduli\n');
fprintf(' change slowly with log resistivity.)\n');

%% Measurements the model cannot reach
% A modulus at or above the mineral modulus, or at or below the fluid
% modulus, has no conductivity in the model: cpdem_Inverse returns NaN
% with a warning, whereas the v1.0 wrapper returns the nearest bound.
ws = warning('off', 'CPDEM:outOfRange');
sBad = cpdem_Inverse([p.k1*1.02, d.k0(5)], 'K', p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, p.eparK);
warning(ws);
fprintf('\nUnreachable modulus -> %s; reachable one -> %.4f S/m\n', mat2str(sBad(1)), sBad(2));

%% v1.0 interface
% (k0, mu0, parameters, EPAR, lower/upper conductivity bounds, modFlag 0/1
% for K/mu, normFlag, paramFlag)
sigK_v1 = cpdem_InvertConduct(d.k0, d.mu0, p.k1,p.mu1,p.sig1,p.k2,p.mu2,p.sig2, p.eparK, p.sig1, p.sig2, 0, 2, 0);
fprintf('cpdem_InvertConduct agrees with cpdem_Inverse to %.1e (relative).\n', max(abs(sigK_v1 - sigK)./sigK));

%% Plot
if cpdem_PlotsEnabled()
    fig = figure('Name', 'Conductivity inversion');
    loglog(d.res0, 1./sigK, 'o', d.res0, 1./sigMu, 's', d.res0, d.res0, 'k-'); grid on
    xlabel('True resistivity (\Omega m)'); ylabel('Inverted resistivity (\Omega m)');
    legend({'from K with \alpha_K', 'from \mu with \alpha_\mu', '1:1'}, 'Location', 'NorthWest');
    cpdem_FinishFigure(fig, 'SolveConduct');
end
