function ok = cpdem_InRange(x, a, b)
% True where x lies between a (inclusive) and b (exclusive), either order
%
%   ok = cpdem_InRange(x, a, b)
%
% a is the background (start) value, which may be attained exactly; b is
% the inclusion (end) value, which a DEM path only reaches at 100 %
% porosity and is therefore excluded. Inputs are equally sized columns.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

lo = min(a, b);
hi = max(a, b);
ok = (x >= lo) & (x <= hi) & (x ~= b);

end
