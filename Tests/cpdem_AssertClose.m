function cpdem_AssertClose(actual, expected, relTol, absTol, what)
% Assert that two arrays agree to within relative or absolute tolerance
%
%   cpdem_AssertClose(actual, expected, relTol)
%   cpdem_AssertClose(actual, expected, relTol, absTol, what)
%
% Passes where |actual - expected| <= absTol + relTol*|expected| for every
% element; otherwise throws an error naming the worst element. NaNs must
% match NaNs.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 4 || isempty(absTol), absTol = 0; end
if nargin < 5, what = 'values'; end

if ~isequal(size(actual), size(expected))
    error('CPDEM:assert', '%s: size mismatch (%s vs %s).', what, ...
        mat2str(size(actual)), mat2str(size(expected)));
end
a = actual(:); e = expected(:);
nanMismatch = xor(isnan(a), isnan(e));
if any(nanMismatch)
    error('CPDEM:assert', '%s: NaN pattern mismatch at element %d.', what, find(nanMismatch, 1));
end
ok = isnan(a) | (abs(a - e) <= absTol + relTol*abs(e));
if ~all(ok)
    [~, worst] = max(abs(a - e) - (absTol + relTol*abs(e)));
    error('CPDEM:assert', ...
        '%s: element %d is %.10g, expected %.10g (rel err %.3g, relTol %.1g, absTol %.1g).', ...
        what, worst, a(worst), e(worst), abs(a(worst) - e(worst))/max(abs(e(worst)), realmin), relTol, absTol);
end

end
