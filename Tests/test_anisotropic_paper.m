function test_anisotropic_paper()
% The closure-based anisotropic cross-property DEM of Cilli, "Anisotropic
% cross-property differential effective medium modelling": depolarisation
% factors (eq. 3), the spheroidal Eshelby tensor (checked through
% Berryman's P and Q), the VRH reduction (eqs 27-32), the isotropic and
% spherical limits, the dry-inclusion closed form (eq. 70), the Section 3
% end points, the cross-property chain rule and the inverse problem.

% --- eq. 3 equals Osborn's depolarisation factor used elsewhere
for a = [0.05 0.1 0.5 0.9 1.1 2 10 50]
    if a > 1
        e = sqrt(1 - 1/a^2); Lp = (1 - e^2)/(2*e^3)*(log((1 + e)/(1 - e)) - 2*e);
    else
        ep = sqrt(1/a^2 - 1); Lp = (1 + ep^2)/ep^3*(ep - atan(ep));
    end
    [~, ~, L] = cpdem_Shape(a);
    cpdem_AssertClose(L, Lp, 1e-11, 0, sprintf('eq. 3 at alpha %g', a));
end

% --- Eshelby tensor: sphere (eq. 25), series/closed-form continuity, and
% Berryman's P, Q from the orientation average of the dilute concentration
nu = 0.25;
S = cpdem_EshelbySpheroid(nu, 1);
a = (7 - 5*nu)/(15*(1 - nu)); b = (5*nu - 1)/(15*(1 - nu)); c = 2*(4 - 5*nu)/(15*(1 - nu));
cpdem_AssertClose([S(1,1) S(1,2) S(1,3) S(3,1) S(3,3) S(4,4) S(6,6)], [a b b b a c c], 1e-14, 0, 'sphere Eshelby tensor');
for e = [-0.02 0.02]
    S1 = cpdem_EshelbySpheroid(nu, 1 + e*(1 - 1e-9)); S2 = cpdem_EshelbySpheroid(nu, 1 + e*(1 + 1e-9));
    cpdem_AssertClose(S1(:), S2(:), 1e-9, 1e-12, sprintf('Eshelby series/closed form at e = %g', e));
end
k1 = 37e9; mu1 = 44e9;
nu1 = (3*k1 - 2*mu1)/(2*(3*k1 + mu1));
C = cpdem_TIStiffness([k1 + 4/3*mu1, k1 - 2/3*mu1, k1 + 4/3*mu1, mu1, mu1]);
for incl = {[2.25e9 0], [10e9 5e9], [80e9 60e9]}
    k2 = incl{1}(1); mu2 = incl{1}(2);
    C2 = cpdem_TIStiffness([k2 + 4/3*mu2, k2 - 2/3*mu2, k2 + 4/3*mu2, mu2, mu2]);
    for alpha = [0.05 0.2 0.9 1 1.2 5 30]
        S = cpdem_EshelbySpheroid(nu1, alpha);
        A = (eye(6) + S*(C\(C2 - C)))\eye(6);
        Aiijj = sum(sum(A(1:3, 1:3)));
        Aijij = trace(A(1:3, 1:3)) + A(4,4) + A(5,5) + A(6,6);
        [P, Q] = cpdem_PQ(k1, mu1, k2, mu2, alpha);
        cpdem_AssertClose([Aiijj/3, (Aijij - Aiijj/3)/5], [P Q], 1e-10, 0, ...
            sprintf('Berryman P,Q from Eshelby, k2 %g mu2 %g alpha %g', k2, mu2, alpha));
    end
end

% --- VRH reduction returns exact moduli for an isotropic stiffness and
% matches the general Voigt/Reuss definitions for a TI one
[K, mu, nuv] = cpdem_VRH(C);
cpdem_AssertClose([K mu nuv], [k1 mu1 nu1], 1e-12, 0, 'VRH of isotropic C');
Ct = cpdem_TIStiffness([60e9 15e9 45e9 12e9 20e9]);
[~, ~, ~, KV, KR, muV, muR] = cpdem_VRH(Ct);
s = inv(Ct);
KVg = sum(sum(Ct(1:3,1:3)))/9;
muVg = (trace(Ct(1:3,1:3)) - (Ct(1,2) + Ct(1,3) + Ct(2,3)) + 3*(Ct(4,4) + Ct(5,5) + Ct(6,6)))/15;
KRg = 1/sum(sum(s(1:3,1:3)));
muRg = 15/(4*trace(s(1:3,1:3)) - 4*(s(1,2) + s(1,3) + s(2,3)) + 3*(s(4,4) + s(5,5) + s(6,6)));
cpdem_AssertClose([KV KR muV muR], [KVg KRg muVg muRg], 1e-12, 0, 'VRH TI formulas');

% --- Paper's numerical setup (Section 3)
sig1 = 1; K1 = 40e9; MU1 = 30e9; K2 = 2.2e9;
phi = [0.1 0.2 0.3];

% Spheres: the model collapses to the isotropic DEM (Bruggeman + Berryman)
[sp, sq] = cpdem_AnisoElecDEM(phi, sig1, 10, 1);
cpdem_AssertClose([sp sq], [cpdem_ElecDEM(phi, sig1, 10, 1) cpdem_ElecDEM(phi, sig1, 10, 1)], 1e-9, 0, 'spheres: electrical = Bruggeman');
el = cpdem_AnisoElasDEM(phi, K1, MU1, K2, 1);
[kk, mm] = cpdem_ElasDEM(phi, K1, MU1, K2, 0, 1);
cpdem_AssertClose([el.K el.mu], [kk mm], 1e-9, 0, 'spheres: elastic = Berryman DEM');
cpdem_AssertClose([el.C11 el.C33], [kk + 4/3*mm, kk + 4/3*mm], 1e-9, 0, 'spheres: C11 = C33 = K + 4/3 mu');
cpdem_AssertClose([el.C44 el.C66], [mm mm], 1e-9, 0, 'spheres: C44 = C66 = mu');
cpdem_AssertClose(el.C13, kk - 2/3*mm, 1e-9, 0, 'spheres: C13 = lambda');

% S = 0: the electrical branch is the randomly oriented isotropic DEM
[sp, sq] = cpdem_AnisoElecDEM(phi, sig1, 10, 0.1, 0);
cpdem_AssertClose([sp sq], [cpdem_ElecDEM(phi, sig1, 10, 0.1) cpdem_ElecDEM(phi, sig1, 10, 0.1)], 1e-9, 0, 'S = 0: randomly oriented electrical DEM');
el0 = cpdem_AnisoElasDEM(phi, K1, MU1, K2, 0.1, 0);
cpdem_AssertClose([el0.C11 el0.C44], [el.C11 el.C44], 1e-9, 0, 'S = 0: spherical-reference elastic branch');

% Dry inclusions, eq. 70
[~, ~, L] = cpdem_Shape(0.1);
[sp, sq] = cpdem_AnisoElecDEM(phi, sig1, 1e-300, 0.1);
cpdem_AssertClose([sp sq], [sig1*(1 - phi).^(1/(1 - L)), sig1*(1 - phi).^(2/(1 + L))], 1e-6, 0, 'dry inclusions, eq. 70');

% Section 3 end points at phi = 0.3, S = 1 (values as quoted in the paper)
paper = struct( ...
    'needles', struct('sig2', 10, 'alpha', 10, 'sp', 3.54, 'sq', 1.74, 'C11', 32.39, 'C33', 54.15, 'C44', 14.90, 'C66', 10.06), ...
    'spheres', struct('sig2', 10, 'alpha', 1,  'sp', 2.03, 'sq', 2.03, 'C11', 41.11, 'C33', 41.11, 'C44', 14.74, 'C66', 14.74), ...
    'discs',   struct('sig2', 10, 'alpha', 0.1,'sp', 1.43, 'sq', 3.17, 'C11', 52.54, 'C33', 9.88,  'C44', 2.70,  'C66', 20.05));
names = fieldnames(paper);
for ii = 1:numel(names)
    p = paper.(names{ii});
    [sp, sq] = cpdem_AnisoElecDEM(0.3, sig1, p.sig2, p.alpha);
    el = cpdem_AnisoElasDEM(0.3, K1, MU1, K2, p.alpha);
    cpdem_AssertClose([sp sq], [p.sp p.sq], 0, 0.006, sprintf('Section 3 conductivities, %s', names{ii}));
    cpdem_AssertClose([el.C11 el.C33 el.C44 el.C66]/1e9, [p.C11 p.C33 p.C44 p.C66], 0, 0.006, sprintf('Section 3 stiffnesses, %s', names{ii}));
end
pore = struct('needles', [10 0.71 0.54], 'discs', [0.1 0.26 0.70], 'spheres', [1 0.62 0.62]);
names = fieldnames(pore);
for ii = 1:numel(names)
    v = pore.(names{ii});
    [sp, sq] = cpdem_AnisoElecDEM(0.3, sig1, 0.05, v(1));
    cpdem_AssertClose([sp sq], v(2:3), 0, 0.006, sprintf('Section 3 pore-like conductivities, %s', names{ii}));
end

% --- Cross-property forward (from sig_par and from sig_perp) reproduces the
% phi-parametrised DEMs, for both closures and partial alignment
for S = [1 0.6]
    for alpha = [10 0.1]
        [sp, sq] = cpdem_AnisoElecDEM(phi, sig1, 10, alpha, S);
        el = cpdem_AnisoElasDEM(phi, K1, MU1, K2, alpha, S);
        f = cpdem_AnisoForward(K1, MU1, sig1, K2, 10, sp, alpha, S, 'par');
        cpdem_AssertClose([f.sigPerp f.phi], [sq phi], 1e-7, 1e-10, sprintf('forward from par: sigPerp, phi (alpha %g, S %g)', alpha, S));
        cpdem_AssertClose([f.C11 f.C13 f.C33 f.C44 f.C66], [el.C11 el.C13 el.C33 el.C44 el.C66], 1e-7, 0, sprintf('forward from par: stiffnesses (alpha %g, S %g)', alpha, S));
        g = cpdem_AnisoForward(K1, MU1, sig1, K2, 10, sq, alpha, S, 'perp');
        cpdem_AssertClose([g.sigPar g.C33 g.C66], [sp el.C33 el.C66], 1e-7, 0, sprintf('forward from perp (alpha %g, S %g)', alpha, S));
    end
end
% exact electrical closure equals cpdem_AlignedCondDEM
[sp, sq] = cpdem_AnisoElecDEM(phi, sig1, 0.05, 0.1, 1, struct('Closure', 'exact'));
[ap, aq] = cpdem_AlignedCondDEM(phi, sig1, 0.05, 0.1);
cpdem_AssertClose([sp sq], [ap aq], 1e-9, 0, 'exact closure = cpdem_AlignedCondDEM');
fx = cpdem_AnisoForward(K1, MU1, sig1, K2, 0.05, sp, 0.1, 1, 'par', 'Closure', 'exact');
cpdem_AssertClose(fx.sigPerp, aq, 1e-7, 0, 'forward with exact closure');

% --- Inverse: stiffness -> conductivities, round trip
[sp, sq] = cpdem_AnisoElecDEM(phi, sig1, 10, 0.1);
el = cpdem_AnisoElasDEM(phi, K1, MU1, K2, 0.1);
for w = {'C33', 'C11', 'C44', 'C66'}
    inv1 = cpdem_AnisoInverse(el.(w{1}), w{1}, K1, MU1, sig1, K2, 10, 0.1);
    cpdem_AssertClose([inv1.sigPar inv1.sigPerp inv1.phi], [sp sq phi], 1e-6, 1e-9, sprintf('inverse from %s', w{1}));
    cpdem_AssertClose(inv1.C13, el.C13, 1e-6, 0, sprintf('inverse from %s: C13', w{1}));
end
ws = warning('off', 'CPDEM:outOfRange');
bad = cpdem_AnisoInverse([K1*3, el.C33(1)], 'C33', K1, MU1, sig1, K2, 10, 0.1);
warning(ws);
assert(isnan(bad.sigPar(1)) && abs(bad.sigPar(2) - sp(1)) < 1e-6*sp(1), 'unreachable stiffness -> NaN');

% --- TI velocities and Thomsen parameters: isotropic case has zero anisotropy
v = cpdem_TIVelocities(el.C11(1), el.C13(1), el.C33(1), el.C44(1), el.C66(1), 2500);
assert(v.vp90 > v.vp0 && v.epsilon > 0, 'aligned discs: faster across the axis');
v0 = cpdem_TIVelocities(80e9, 20e9, 80e9, 30e9, 30e9, 2500);
cpdem_AssertClose([v0.epsilon v0.gamma v0.delta], [0 0 0], 0, 1e-14, 'isotropic Thomsen parameters');

end
