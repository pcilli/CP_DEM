function cpdem_AnisoCheck(alpha, S, closure)
% Validate the shape and alignment inputs of the anisotropic functions
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || alpha <= 0
    error('CPDEM:badInput', 'alpha must be a positive finite scalar aspect ratio.');
end
if ~isscalar(S) || ~isreal(S) || S < 0 || S > 1
    error('CPDEM:badInput', 'The alignment order parameter S must lie in [0, 1].');
end
if strcmp(closure, 'exact') && S ~= 1
    error('CPDEM:badInput', 'The exact closure is only defined for perfect alignment (S = 1).');
end

end
