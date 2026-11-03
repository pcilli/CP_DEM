function [pUp, pLow, phiMin, phiMax, pOn1, pOn2] = cpdem_CondCondHS(sig, sig1, sig2, p1, p2)
% Cross-property Hashin-Shtrikman bounds between two conductivities
%
%   [pUp, pLow] = cpdem_CondCondHS(sig, sig1, sig2, p1, p2)
%   [pUp, pLow, phiMin, phiMax, pOn1, pOn2] = cpdem_CondCondHS(...)
%
% Bounds on a composite's property p (e.g. thermal conductivity) given a
% measurement of another Laplace-type property sig (e.g. electrical
% conductivity), following the stratagem of Carcione et al. (2007) as used
% by Cilli and Chapman (arXiv:2102.12434, Section 3): the upper and lower
% two-phase Hashin-Shtrikman bounds on sig are each rearranged for the
% inclusion volume fraction, and those fractions are substituted into the
% upper and lower Hashin-Shtrikman bounds on p respectively. Either phase
% may be the more conductive one, and the bounds for sig given p follow by
% swapping the arguments.
%
% Each Hashin-Shtrikman bound is attained by the coated-sphere assemblage
% built on one phase, and the same assemblage attains the same-type bound
% of every Laplace property, which is why the bounds are paired type with
% type at a common volume fraction: pOn1 is the curve built on phase 1 and
% pOn2 the curve built on phase 2. The differential scheme with disk-
% shaped inclusions reproduces the bound built on the inclusion phase for
% each property (Norris 1985), so the EPAR = 0 trend of
% cpdem_ForwardConduct lies exactly on pOn2, as noted in Section 4 of the
% paper. Which of the two curves is the upper envelope in the
% cross-property plane depends on the phase ordering: when the inclusion
% phase is lower in both properties (copper-graphite), the single-property
% lower bound built on it is the UPPER envelope of p versus sig, because a
% composite on that bound needs less inclusion to reach a given sig. pUp
% and pLow are therefore returned as the larger and smaller of pOn1 and
% pOn2.
%
% Inputs
%   sig      : measured values of the first property (array), between
%              sig1 and sig2
%   sig1, p1 : values of phase 1 (fraction 1 - phi, the host of the DEM)
%   sig2, p2 : values of phase 2 (fraction phi, the inclusions)
%
% Outputs (shape of sig): upper and lower envelopes of p, the volume
% fractions obtained from the two conductivity bounds as an interval
% [phiMin, phiMax], and the two paired curves pOn1 and pOn2.
%
% For two phases with properties a (fraction 1 - phi) and b (fraction phi)
% the Hashin-Shtrikman bounds are a + phi/(1/(b - a) + (1 - phi)/(3a))
% (built on phase a) and b + (1 - phi)/(1/(a - b) + phi/(3b)) (built on
% phase b); the one built on the more conductive phase is the upper bound.
%
% References:
% Carcione, J. M., Ursin, B., & Nordskag, J. I. (2007). Cross-property
% relations between electrical conductivity and the seismic velocity of
% rocks. Geophysics, 72(5), E193-E204.
%
% Cilli, P. A., and Chapman, M. (2021). Linking electrical and thermal
% conductivity through cross-property inclusion modelling. arXiv:2102.12434
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

cpdem_CheckPositive('sig', sig, 'sig1', sig1, 'sig2', sig2, 'p1', p1, 'p2', p2);
if ~all([isscalar(sig1), isscalar(sig2), isscalar(p1), isscalar(p2)])
    error('CPDEM:badInput', 'Phase properties must be scalars.');
end
lo = min(sig1, sig2); hi = max(sig1, sig2);
if any(sig(:) < lo | sig(:) > hi)
    warning('CPDEM:outOfRange', ...
        'Some values lie outside the interval between the phase properties; bounds are not meaningful there.');
end

% Volume fractions at which each conductivity bound equals the measurement.
% The bound built on phase 1 is the upper bound when phase 1 is the more
% conductive phase, and vice versa.
phi1 = (sig - sig1).*(2*sig1 + sig2)./((sig2 - sig1).*(sig + 2*sig1));      % bound built on phase 1
phi2 = 1 - (sig - sig2).*(2*sig2 + sig1)./((sig1 - sig2).*(sig + 2*sig2));  % bound built on phase 2
phi1 = min(max(phi1, 0), 1);
phi2 = min(max(phi2, 0), 1);
if sig1 >= sig2
    phiFromUpper = phi1; phiFromLower = phi2;
else
    phiFromUpper = phi2; phiFromLower = phi1;
end
phiMin = min(phiFromUpper, phiFromLower);
phiMax = max(phiFromUpper, phiFromLower);

% Same-type bounds on p at those volume fractions: the curve built on
% phase 1 uses the fraction from the sig bound built on phase 1, etc.
pOn1 = hsBuiltOn1(phi1, p1, p2);
pOn2 = hsBuiltOn2(phi2, p1, p2);
pUp  = max(pOn1, pOn2);
pLow = min(pOn1, pOn2);

end

function b = hsBuiltOn1(phi, a1, a2)
% Hashin-Shtrikman bound with phase 1 (fraction 1 - phi) as the coating
b = a1 + phi./(1./(a2 - a1) + (1 - phi)./(3*a1));
end

function b = hsBuiltOn2(phi, a1, a2)
% Hashin-Shtrikman bound with phase 2 (fraction phi) as the coating
b = a2 + (1 - phi)./(1./(a1 - a2) + phi./(3*a2));
end
