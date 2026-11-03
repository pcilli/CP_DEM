function [Gpar, Gperp] = cpdem_AlignedIncrement(sigPar, sigPerp, sig2, epar)
% DEM increments for aligned spheroids in a transversely isotropic host
%
%   [Gpar, Gperp] = cpdem_AlignedIncrement(sigPar, sigPerp, sig2, epar)
%
% Returns (1 - phi) d sig_par/d phi and (1 - phi) d sig_perp/d phi for
% spheroidal inclusions of aspect ratio EPAR whose unique axes are all
% aligned with the symmetry axis of a transversely isotropic host with
% conductivity sig_par along that axis and sig_perp in the plane normal
% to it. sig2 is the inclusion conductivity: a scalar, or
% [sig2_axis sig2_plane] for an inclusion material that is itself
% anisotropic with its principal axes along the spheroid axes.
%
% Derivation (affine transformation). Stretching each coordinate by
% 1/sqrt(sig_i) turns the host into an isotropic unit conductor, the
% field ratios inside and outside the inclusion are unchanged component by
% component, the volume fraction is unchanged, the inclusion conductivity
% becomes sig2_i/sig_i, and the spheroid remains a spheroid with
%
%   EPAR' = EPAR sqrt(sig_perp/sig_par).
%
% The dilute field-concentration factor along axis i is therefore
% 1/(1 + L_i(EPAR') (sig2_i/sig_i - 1)), with L_axis = L(EPAR') and
% L_plane = (1 - L(EPAR'))/2 (Osborn 1945), and the DEM increments are
%
%   G_i = sig_i (sig2_i - sig_i) / (sig_i + L_i(EPAR') (sig2_i - sig_i)).
%
% Checks (see test_anisotropic_conductivity): aligned disks give the
% series (harmonic) average along the axis and the parallel (arithmetic)
% average in the plane; aligned needles give the arithmetic average along
% the axis and the two-dimensional Bruggeman model in the plane; spheres
% give the isotropic Bruggeman model in both.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if isscalar(sig2)
    sa = sig2; sp = sig2;
elseif numel(sig2) == 2
    sa = sig2(1); sp = sig2(2);
else
    error('CPDEM:badInput', 'Inclusion conductivity must be a scalar or [sig_axis sig_plane].');
end

aeff = epar*sqrt(sigPerp/sigPar);
[~, ~, L] = cpdem_Shape(aeff);
Lp = (1 - L)/2;

Gpar  = sigPar.*(sa - sigPar)./(sigPar + L*(sa - sigPar));
Gperp = sigPerp.*(sp - sigPerp)./(sigPerp + Lp*(sp - sigPerp));

end
