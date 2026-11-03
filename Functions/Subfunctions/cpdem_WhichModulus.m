function useK = cpdem_WhichModulus(which)
% Interpret a modulus selector: 'K'/'k'/0 -> true, 'mu'/'G'/1 -> false
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if ischar(which)
    switch lower(strtrim(which))
        case {'k', 'bulk'}
            useK = true;
        case {'mu', 'g', 'shear'}
            useK = false;
        otherwise
            error('CPDEM:badInput', 'Modulus selector must be ''K'' or ''mu'' (got ''%s'').', which);
    end
elseif isnumeric(which) && isscalar(which) && any(which == [0 1])
    useK = (which == 0);   % v1.0 modFlag convention
elseif islogical(which) && isscalar(which)
    useK = which;
else
    error('CPDEM:badInput', 'Modulus selector must be ''K'', ''mu'', 0 or 1.');
end

end
