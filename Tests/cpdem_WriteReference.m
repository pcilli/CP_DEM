function cpdem_WriteReference(outfile)
% Write reference values for the Python port's cross-check (python/test_cpdem.py)
%
%   cpdem_WriteReference            writes python/reference_values.json
%   cpdem_WriteReference(outfile)   writes to the given file
%
% Requires jsonencode (MATLAB R2016b+, GNU Octave 7+).
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(fileparts(here), 'Functions')));
if nargin < 1
    outfile = fullfile(fileparts(here), 'python', 'reference_values.json');
end

p = cpdem_TestParams();
ref.generator = sprintf('CP_DEM v%s, %s, %s', cpdem_Version(), version(), datestr(now, 'yyyy-mm-dd'));
ref.params = struct('k1', p.k1, 'mu1', p.mu1, 'sig1', p.sig1, 'k2', p.k2, ...
                    'mu2', p.mu2, 'sig2', p.sig2, 'sig0', p.sig0);
ref.thermal = struct('lam1', 7.7, 'lam2', 0.6);
epars = [0.1 0.95 1 1.05 12.8 16.4 100];
cases = cell(1, numel(epars));
for ii = 1:numel(epars)
    [k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, epars(ii));
    lam = cpdem_ForwardConduct(7.7, p.sig1, 0.6, p.sig2, p.sig0, epars(ii));
    cases{ii} = struct('epar', epars(ii), 'k0', k0, 'mu0', mu0, 'phi0', phi0, 'lam', lam);
end
ref.cases = [cases{:}];

txt = jsonencode(ref);
fid = fopen(outfile, 'w');
fprintf(fid, '%s\n', txt);
fclose(fid);
fprintf('Wrote %s\n', outfile);

end
