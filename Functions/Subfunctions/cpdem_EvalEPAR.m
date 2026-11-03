function alpha = cpdem_EvalEPAR(epar, phi)
% Evaluate the (possibly porosity-dependent) equivalent pore aspect ratio
%
%   alpha = cpdem_EvalEPAR(epar, phi)
%
% EPAR may be
%   * a positive scalar: a constant aspect ratio (Cilli & Chapman 2021);
%   * a function handle alpha = epar(phi) giving the aspect ratio as a
%     function of the inclusion volume fraction along the DEM path (the
%     porosity when pores are the inclusions, as in the cross-property DEM;
%     the grain fraction in a grain-inclusion model), e.g. the power law of
%     Cilli & Chapman (2020),
%     epar = @(phi) max(a*phi.^b, floor). Because porosity is carried as
%     a state variable along every DEM path in this toolbox, such a
%     relation can be used in all forward and inverse calculations. The
%     law must return a strictly positive value at every porosity on the
%     path, including phi = 0 where every path starts: a pure power law is
%     zero there, so give it a small floor as above.
%
% References:
% Cilli, P. A., and Chapman, M. (2020). The power-law relation between
% inclusion aspect ratio and porosity: Implications for electrical and
% elastic modeling. Journal of Geophysical Research: Solid Earth, 125,
% e2019JB019187, DOI:10.1029/2019JB019187
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if isa(epar, 'function_handle')
    % ODE solvers evaluate trial states that may lie marginally outside
    % [0, 1); evaluate the law at the nearest physical porosity.
    alpha = epar(min(max(phi, 0), 1));
    if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || alpha <= 0
        error('CPDEM:badEPAR', ...
            ['EPAR function handle must return a positive finite scalar; it returned %s at phi = %g. ', ...
             'A power law a*phi^b is zero at phi = 0, where every DEM path starts: use e.g. @(phi) max(a*phi.^b, 1e-3).'], ...
            mat2str(alpha), phi);
    end
elseif isnumeric(epar) && isscalar(epar)
    alpha = epar;
else
    error('CPDEM:badEPAR', ...
        'EPAR must be a positive scalar or a function handle of porosity.');
end

end
