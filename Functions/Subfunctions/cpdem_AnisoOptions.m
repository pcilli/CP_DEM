function [opts, closure] = cpdem_AnisoOptions(varargin)
% Options for the anisotropic functions: cpdem_Options plus 'Closure'
%
%   [opts, closure] = cpdem_AnisoOptions(opts_struct_or_pairs)
%
% Separates the extra option Closure ('isotropic-host' default, or 'exact')
% from the standard cpdem_Options fields.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

closure = 'isotropic-host';
args = varargin;
if numel(args) == 1 && iscell(args{1}), args = args{1}; end
if numel(args) == 1 && isempty(args{1}), args = {}; end
if numel(args) >= 1 && isempty(args{1}) && ~isstruct(args{1}), args(1) = []; end
keep = true(size(args));
if ~isempty(args) && isstruct(args{1})
    s = args{1};
    if isfield(s, 'Closure')
        closure = s.Closure;
        s = rmfield(s, 'Closure');
    end
    args{1} = s;
end
for ii = 1:numel(args)
    if ischar(args{ii}) && strcmpi(args{ii}, 'Closure') && ii < numel(args)
        closure = args{ii+1};
        keep([ii ii+1]) = false;
    end
end
opts = cpdem_Options(args{keep});
closure = lower(closure);
if ~any(strcmp(closure, {'isotropic-host', 'exact'}))
    error('CPDEM:badOptions', 'Closure must be ''isotropic-host'' or ''exact''.');
end

end
