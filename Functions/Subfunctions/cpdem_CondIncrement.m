function G = cpdem_CondIncrement(sig, sig2, epar)
% DEM increment of a conductivity-type property for randomly oriented spheroids
%
%   G = cpdem_CondIncrement(sig, sig2, epar)
%
% Returns G = (1 - phi) d sig/d phi for randomly oriented spheroidal
% inclusions of aspect ratio EPAR in an isotropic host of conductivity sig:
%
%   isotropic inclusions (scalar sig2):
%       G = (sig2 - sig) m-bar(sig, sig2, EPAR)          (Cilli & Chapman 2021)
%
%   anisotropic inclusions (sig2 = [sig2_axis sig2_plane], the principal
%   conductivities of the inclusion material along the unique axis of the
%   spheroid and in the plane of its two equal axes, e.g. the through-plane
%   and in-plane conductivities of a graphite flake):
%       G = (1/3) sig [ (sig2_axis  - sig)/(sig + L (sig2_axis  - sig))
%                   + 2 (sig2_plane - sig)/(sig + L' (sig2_plane - sig)) ],
%       L' = (1 - L)/2,
%
% which is the orientation average of the dilute field-concentration
% factors of an ellipsoidal inclusion whose principal conductivity axes
% coincide with its geometric axes (Landau & Lifshitz 1960, section 8;
% Sen et al. 1981). With sig2_axis = sig2_plane it reduces to the
% isotropic form above. The anisotropic option is an extension of the
% toolbox for inclusions such as graphite, whose intrinsic anisotropy
% Cilli and Chapman (2021, arXiv:2102.12434) discuss.
%
% sig may be an array; sig2 is a scalar or a 2-vector; EPAR is a scalar.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if isscalar(sig2)
    G = (sig2 - sig).*mbar(sig, sig2, epar);
elseif numel(sig2) == 2
    [~, ~, L] = cpdem_Shape(epar);
    Lp = (1 - L)/2;
    sa = sig2(1);
    sp = sig2(2);
    G = (1/3)*sig.*((sa - sig)./(sig + L*(sa - sig)) + 2*(sp - sig)./(sig + Lp*(sp - sig)));
else
    error('CPDEM:badInput', ...
        'Inclusion conductivity must be a scalar or [sig_axis sig_plane].');
end

end
