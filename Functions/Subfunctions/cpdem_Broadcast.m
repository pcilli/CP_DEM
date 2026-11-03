function varargout = cpdem_Broadcast(n, varargin)
% Expand scalar parameters to n-by-1 columns and check vector sizes
%
%   [p1, p2, ...] = cpdem_Broadcast(n, p1, p2, ...)
%
% Each input must be a real scalar or an array with n elements. Scalars are
% repeated n times; arrays are reshaped to n-by-1 columns.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

varargout = cell(size(varargin));
for ii = 1:numel(varargin)
    p = varargin{ii};
    if ~isnumeric(p) || ~isreal(p)
        error('CPDEM:badInput', 'Parameter %d must be real and numeric.', ii);
    end
    if isscalar(p)
        varargout{ii} = repmat(p, n, 1);
    elseif numel(p) == n
        varargout{ii} = p(:);
    else
        error('CPDEM:badInput', ...
            'Parameter %d must be a scalar or have %d elements (one per measurement); it has %d.', ...
            ii, n, numel(p));
    end
end

end
