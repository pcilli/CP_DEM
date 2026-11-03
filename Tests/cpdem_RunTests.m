function results = cpdem_RunTests(pattern)
% Run the Cross-Property DEM Toolbox test suite (MATLAB and GNU Octave)
%
%   cpdem_RunTests            run every test_*.m in this folder
%   cpdem_RunTests('inverse') run only tests whose name contains 'inverse'
%   results = cpdem_RunTests  also return a struct array of outcomes
%
% Each test is a function file test_<name>.m that throws an error on
% failure (use assert or cpdem_AssertClose). The runner prints a summary
% and, when called without an output argument, raises an error if any test
% failed so that it can be used from a CI job:
%
%   octave --no-gui -q --eval "cd Tests; cpdem_RunTests"
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 1
    pattern = '';
end

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(fullfile(root, 'Functions')));
addpath(here);

files = dir(fullfile(here, 'test_*.m'));
names = sort({files.name});
if ~isempty(pattern)
    names = names(~cellfun(@isempty, strfind(names, pattern)));
end

results = struct('name', {}, 'passed', {}, 'message', {}, 'seconds', {});
fprintf('Cross-Property DEM Toolbox v%s test suite (%s)\n', cpdem_Version(), version());
fprintf('%s\n', repmat('-', 1, 60));
for ii = 1:numel(names)
    name = names{ii}(1:end-2);
    t0 = tic;
    try
        feval(name);
        ok = true;
        msg = '';
    catch err
        ok = false;
        msg = err.message;
    end
    dt = toc(t0);
    results(end+1) = struct('name', name, 'passed', ok, 'message', msg, 'seconds', dt); %#ok<AGROW>
    if ok
        fprintf('PASS  %-40s (%.2f s)\n', name, dt);
    else
        fprintf('FAIL  %-40s (%.2f s)\n      %s\n', name, dt, strrep(msg, sprintf('\n'), sprintf('\n      ')));
    end
end

nPass = nnz([results.passed]);
nFail = numel(results) - nPass;
fprintf('%s\n%d passed, %d failed, %.1f s total\n', repmat('-', 1, 60), ...
    nPass, nFail, sum([results.seconds]));

if nargout == 0
    clear results
    if nFail > 0
        error('CPDEM:testsFailed', '%d test(s) failed.', nFail);
    end
end

end
