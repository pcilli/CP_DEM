%% Example: Anisotropic (aligned-spheroid) cross-property DEM
%
% Reproduces the numerical analysis of Cilli, "Anisotropic cross-property
% differential effective medium modelling" (Section 3): the closure-based
% transversely isotropic (TI) cross-property DEM for aligned spheroidal
% inclusions, with the trajectories of the five TI stiffnesses against
% the axial and transverse conductivities for needles (alpha = 10),
% spheres (alpha = 1) and discs (alpha = 0.1), conductive (sigma2 = 10 S/m)
% and pore-like (sigma2 = 0.05 S/m) inclusions, and alignment order
% parameters S in {0, 0.3, 0.6, 0.85, 1} (Figures 1-6 of the paper). The
% end points at phi = 0.3 are compared with the values quoted in the paper,
% the inverse problem is demonstrated, and the paper's electrical closure
% is compared with the exact coaxial solution.
%
% Functions: cpdem_AnisoElecDEM, cpdem_AnisoElasDEM, cpdem_AnisoForward,
% cpdem_AnisoInverse, cpdem_TIVelocities.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

%% Numerical setup of the paper
sig1 = 1;  K1 = 40e9;  MU1 = 30e9;  K2 = 2.2e9;
phi = linspace(0, 0.3, 31); phi(1) = 1e-9;
shapes = struct('name', {'needles', 'spheres', 'discs'}, 'alpha', {10, 1, 0.1});
contrasts = [10 0.05];
Svals = [0 0.3 0.6 0.85 1];

%% Trajectories (volume-fraction DEMs) and end-point table
fprintf('%-10s %-6s %8s %8s %8s %8s %8s %8s\n', 'inclusion', 'sig2', 'sig_par', 'sig_perp', 'C11', 'C33', 'C44', 'C66');
traj = struct();
for ic = 1:numel(contrasts)
    sig2 = contrasts(ic);
    for is = 1:numel(shapes)
        a = shapes(is).alpha;
        for iS = 1:numel(Svals)
            S = Svals(iS);
            [sp, sq] = cpdem_AnisoElecDEM(phi, sig1, sig2, a, S);
            el = cpdem_AnisoElasDEM(phi, K1, MU1, K2, a, S);
            traj(ic, is, iS).sp = sp; traj(ic, is, iS).sq = sq; traj(ic, is, iS).el = el;
        end
        t = traj(ic, is, end);
        fprintf('%-10s %-6g %8.3f %8.3f %8.2f %8.2f %8.2f %8.2f\n', shapes(is).name, sig2, ...
            t.sp(end), t.sq(end), t.el.C11(end)/1e9, t.el.C33(end)/1e9, t.el.C44(end)/1e9, t.el.C66(end)/1e9);
    end
end
fprintf(['\nPaper (S = 1, phi = 0.3): needles 3.54/1.74, C11 32.39, C33 54.15, C44 14.90, C66 10.06;\n' ...
         'spheres 2.03, 41.11, 14.74; discs 1.43/3.17, C11 52.54, C33 9.88, C44 2.70, C66 20.05;\n' ...
         'pore-like: needles 0.71/0.54, spheres 0.62, discs 0.26/0.70 (GPa, S/m).\n']);

%% The porosity-free cross-property system: stiffness from a measured conductivity
sigMeas = [1.5 2.0 2.5 3.0];                   % axial conductivities, conductive needles
f = cpdem_AnisoForward(K1, MU1, sig1, K2, 10, sigMeas, 10, 1, 'par');
fprintf('\nTI stiffness from sigma_par (needles, S = 1):\n');
fprintf('%8s %9s %8s %8s %8s %8s %8s\n', 'sig_par', 'sig_perp', 'C11', 'C13', 'C33', 'C44', 'C66');
for ii = 1:numel(sigMeas)
    fprintf('%8.2f %9.3f %8.2f %8.2f %8.2f %8.2f %8.2f\n', sigMeas(ii), f.sigPerp(ii), ...
        f.C11(ii)/1e9, f.C13(ii)/1e9, f.C33(ii)/1e9, f.C44(ii)/1e9, f.C66(ii)/1e9);
end
v = cpdem_TIVelocities(f.C11, f.C13, f.C33, f.C44, f.C66, 2500);
fprintf('Thomsen epsilon %s, gamma %s\n', mat2str(v.epsilon, 3), mat2str(v.gamma, 3));

%% Inverse: conductivity tensor from a measured C33 (Section 2.7)
inv1 = cpdem_AnisoInverse(f.C33, 'C33', K1, MU1, sig1, K2, 10, 10, 1);
fprintf('\nRound trip through C33: max relative error in sigma_par %.1e\n', max(abs(inv1.sigPar - sigMeas)./sigMeas));

%% The paper's electrical closure versus the exact coaxial TI solution
fprintf('\nElectrical closure (paper) vs exact affine-stretch solution at phi = 0.3, S = 1:\n');
for sig2 = contrasts
    for a = [10 0.1]
        [cp, cq] = cpdem_AnisoElecDEM(0.3, sig1, sig2, a);
        [ep, eq] = cpdem_AnisoElecDEM(0.3, sig1, sig2, a, 1, 'Closure', 'exact');
        fprintf('  sig2 = %5g alpha = %4g: closure %.4f / %.4f   exact %.4f / %.4f\n', sig2, a, cp, cq, ep, eq);
    end
end

%% Plot: Figures 1 and 3 of the paper (conductive needles and discs)
if cpdem_PlotsEnabled()
    comps = {'C11', 'C13', 'C33', 'C44', 'C66'};
    for is = [1 3]
        fig = figure('Name', sprintf('Anisotropic cross-property DEM: conductive %s', shapes(is).name));
        for ic2 = 1:5
            for side = 1:2
                subplot(5, 2, 2*(ic2 - 1) + side); hold on
                for iS = 1:numel(Svals)
                    t = traj(1, is, iS);
                    if side == 1, x = t.sq; else, x = t.sp; end
                    plot(x, t.el.(comps{ic2})/1e9, 'LineWidth', 1.2);
                end
                grid on
                if side == 1, xlabel('\sigma_\perp (S/m)'); else, xlabel('\sigma_{||} (S/m)'); end
                ylabel([comps{ic2} ' (GPa)']);
                if ic2 == 1 && side == 1
                    legend(arrayfun(@(s) sprintf('S = %.2f', s), Svals, 'UniformOutput', false), 'Location', 'Best');
                end
            end
        end
        cpdem_FinishFigure(fig, sprintf('AnisotropicCrossProperty_%s', shapes(is).name));
    end
end
