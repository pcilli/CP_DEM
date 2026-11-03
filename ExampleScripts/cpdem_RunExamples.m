function cpdem_RunExamples(savefigs)
% Run every example script (optionally headless) and report failures
%
%   cpdem_RunExamples          run all examples; plots shown if available
%   cpdem_RunExamples(true)    also save the figures to ExampleScripts/output
%
% Set the environment variable CPDEM_NOPLOT=1 (or the global variable
% CPDEM_NOPLOT = true) to skip plotting, e.g. on a CI machine without a
% display. The function raises an error if any example fails.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

global CPDEM_SAVEFIGS CPDEM_NOPLOT
if nargin < 1
    savefigs = false;
end
CPDEM_SAVEFIGS = savefigs;

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(fileparts(here), 'Functions')));
addpath(here);

% On Octave without any graphics toolkit (e.g. a CI machine), fall back
% to no plotting rather than failing.
if exist('OCTAVE_VERSION', 'builtin') && isempty(available_graphics_toolkits())
    CPDEM_NOPLOT = true;
end

if ~cpdem_PlotsEnabled()
    fprintf('(plotting disabled)\n');
end

files = dir(fullfile(here, 'cpdem_*_Example.m'));
names = sort({files.name});
nFail = 0;
for ii = 1:numel(names)
    name = names{ii}(1:end-2);
    fprintf('\n==== %s ====\n', name);
    t0 = tic;
    try
        evalin('base', ['clearvars -except CPDEM_SAVEFIGS CPDEM_NOPLOT; ' name ';']);
        fprintf('---- %s finished in %.1f s\n', name, toc(t0));
    catch err
        nFail = nFail + 1;
        fprintf('---- %s FAILED: %s\n', name, err.message);
    end
end
if nFail > 0
    error('CPDEM:examplesFailed', '%d example script(s) failed.', nFail);
end

end
