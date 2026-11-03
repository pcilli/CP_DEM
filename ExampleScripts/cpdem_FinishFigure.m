function cpdem_FinishFigure(fig, name)
% Save a figure of an example script as PNG when requested
%
%   cpdem_FinishFigure(fig, 'name')
%
% If the global variable CPDEM_SAVEFIGS is true, the figure is written to
% ExampleScripts/output/<name>.png (the folder is created if needed).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

global CPDEM_SAVEFIGS
if isempty(CPDEM_SAVEFIGS) || ~CPDEM_SAVEFIGS || ~cpdem_PlotsEnabled()
    return
end
outdir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
print(fig, '-dpng', '-r110', fullfile(outdir, [name '.png']));

end
