function cpdem_WriteAnisoReference(outfile)
% Write reference values of the anisotropic functions for the Python port
%
%   cpdem_WriteAnisoReference            writes python/reference_aniso.json
%
% Part of the Cross-Property DEM Toolbox (CP_DEM), anisotropic extension.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(fileparts(here), 'Functions')));
if nargin < 1
    outfile = fullfile(fileparts(here), 'python', 'reference_aniso.json');
end
ref.generator = sprintf('CP_DEM anisotropic extension, %s, %s', version(), datestr(now, 'yyyy-mm-dd'));
ref.params = struct('sig1', 1, 'k1', 40e9, 'mu1', 30e9, 'k2', 2.2e9, 'phi', [0.05 0.15 0.3]);
cases = {};
for sig2 = [10 0.05]
    for alpha = [10 1 0.1]
        for S = [1 0.6 0]
            [sp, sq] = cpdem_AnisoElecDEM(ref.params.phi, 1, sig2, alpha, S);
            el = cpdem_AnisoElasDEM(ref.params.phi, 40e9, 30e9, 2.2e9, alpha, S);
            f = cpdem_AnisoForward(40e9, 30e9, 1, 2.2e9, sig2, sp, alpha, S, 'par');
            cases{end+1} = struct('sig2', sig2, 'alpha', alpha, 'S', S, 'sigPar', sp, 'sigPerp', sq, ...
                'C11', el.C11, 'C13', el.C13, 'C33', el.C33, 'C44', el.C44, 'C66', el.C66, ...
                'fwdC33', f.C33, 'fwdSigPerp', f.sigPerp); %#ok<AGROW>
        end
    end
end
ref.cases = [cases{:}];
S = cpdem_EshelbySpheroid(0.3, 0.2); ref.eshelby_nu03_a02 = S(:)';
S = cpdem_EshelbySpheroid(0.3, 1.01); ref.eshelby_nu03_a101 = S(:)';
fid = fopen(outfile, 'w'); fprintf(fid, '%s\n', jsonencode(ref)); fclose(fid);
fprintf('Wrote %s\n', outfile);

end
