function cpdem_CheckPositive(varargin)
% Validate that named parameters are real, finite and positive
%
%   cpdem_CheckPositive('name1', value1, 'name2', value2, ...)
%
% Values may be scalars or arrays. Use 'name>=0' style by passing a name
% ending in '>=0' to allow zeros (e.g. an inviscid fluid shear modulus).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

for ii = 1:2:numel(varargin)
    name = varargin{ii};
    v = varargin{ii+1};
    allowZero = numel(name) > 3 && strcmp(name(end-2:end), '>=0');
    if allowZero
        name = name(1:end-3);
    end
    if ~isnumeric(v) || ~isreal(v) || any(~isfinite(v(:)))
        error('CPDEM:badInput', '%s must be real, numeric and finite.', name);
    end
    if allowZero
        if any(v(:) < 0)
            error('CPDEM:badInput', '%s must be non-negative.', name);
        end
    elseif any(v(:) <= 0)
        error('CPDEM:badInput', ...
            '%s must be strictly positive (a perfectly insulating host or a zero modulus is not supported by DEM).', name);
    end
end

end
