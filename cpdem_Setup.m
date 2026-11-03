function cpdem_Setup()
% Add the Cross-Property DEM Toolbox to the MATLAB/Octave path
%
%   cpdem_Setup
%
% Adds the Functions folder (and its subfolders) of the toolbox to the
% search path for the current session. Call this once, or add it to your
% startup.m. To make the change permanent in MATLAB use savepath.
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'Functions')));
addpath(fullfile(root, 'Tests'));
addpath(fullfile(root, 'ExampleScripts'));
fprintf('Cross-Property DEM Toolbox v%s added to the path.\n', cpdem_Version());

end
