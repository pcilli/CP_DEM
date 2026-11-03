function tf = cpdem_PlotsEnabled()
% True unless plotting has been disabled for headless runs
%
% Plotting is disabled when the global variable CPDEM_NOPLOT is true or
% the environment variable CPDEM_NOPLOT is set to 1 (as cpdem_RunExamples
% does on a machine without graphics).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

global CPDEM_NOPLOT
tf = true;
if ~isempty(CPDEM_NOPLOT) && CPDEM_NOPLOT
    tf = false;
end
if strcmp(getenv('CPDEM_NOPLOT'), '1')
    tf = false;
end

end
