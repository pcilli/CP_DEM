function r = cpdem_AnisoRates(sp, sq, c5, prm)
% Rates (1 - phi) d/dphi of the anisotropic DEM state [sigPar sigPerp C11 C13 C33 C44 C66]
%
%   r = cpdem_AnisoRates(sigPar, sigPerp, c5, prm)
%
% Implements the closure-based model of Cilli, "Anisotropic cross-property
% differential effective medium modelling":
%
%   electrical (eqs 7-12 with the partial-alignment surrogate eqs 54-56)
%     r_p = (sig2 - sig_p) mbar_p^(S),  p = par, perp,
%     mbar_par  = sig_par /(sig_par  + L (sig2 - sig_par)),
%     mbar_perp = sig_perp/(sig_perp + (1 - L)/2 (sig2 - sig_perp)),
%     mbar_rand = (2 mbar_perp + mbar_par)/3,
%     mbar_p^(S) = mbar_rand + S (mbar_p - mbar_rand);
%   or, with closure 'exact', the exact coaxial transversely isotropic
%   increments of cpdem_AlignedIncrement (S = 1 only);
%
%   elastic (eqs 20-39 with eq. 57)
%     M = (C2 - C) A,  A = [I + S_cl C^-1 (C2 - C)]^-1,
%     S_cl = (1 - S) S_sphere(nu) + S S_spheroid(nu, alpha),
%     nu from cpdem_VRH(C),  C2 = fluid (K2 in the 3x3 normal block),
%   returning the entries M11, M13, M33, M44, M66.
%
% prm fields: sig2, k2, alpha, S, L (axial depolarisation factor for alpha),
% closure ('isotropic-host' | 'exact'). Inputs in SI.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

% --- electrical part
switch prm.closure
    case 'isotropic-host'
        L = prm.L;
        mPar  = sp/(sp + L*(prm.sig2 - sp));
        mPerp = sq/(sq + (1 - L)/2*(prm.sig2 - sq));
        mRand = (2*mPerp + mPar)/3;
        dsp = (prm.sig2 - sp)*(mRand + prm.S*(mPar - mRand));
        dsq = (prm.sig2 - sq)*(mRand + prm.S*(mPerp - mRand));
    case 'exact'
        [dsp, dsq] = cpdem_AlignedIncrement(sp, sq, prm.sig2, prm.alpha);
    otherwise
        error('CPDEM:badOptions', 'closure must be ''isotropic-host'' or ''exact''.');
end

% --- elastic part
C = cpdem_TIStiffness(c5);
[~, ~, nu] = cpdem_VRH(C);
if prm.S == 1
    Scl = cpdem_EshelbySpheroid(nu, prm.alpha);
elseif prm.S == 0
    Scl = cpdem_EshelbySpheroid(nu, 1);
else
    Scl = (1 - prm.S)*cpdem_EshelbySpheroid(nu, 1) + prm.S*cpdem_EshelbySpheroid(nu, prm.alpha);
end
C2 = zeros(6); C2(1:3, 1:3) = prm.k2;
dC = C2 - C;
A = (eye(6) + Scl*(C\dC))\eye(6);
M = dC*A;

r = [dsp; dsq; M(1,1); M(1,3); M(3,3); M(4,4); M(6,6)];

end
