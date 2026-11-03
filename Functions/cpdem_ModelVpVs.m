function out = cpdem_ModelVpVs(k1,mu1,sig1,k2,mu2,sig2,sig0,kepar,muepar,dens,opts)
% Model P- and S-wave velocities and Vp/Vs from electrical conductivity
%
%   out = cpdem_ModelVpVs(k1,mu1,sig1,k2,mu2,sig2,sig0,kepar,muepar)
%   out = cpdem_ModelVpVs(k1,mu1,sig1,k2,mu2,sig2,sig0,kepar,muepar,dens)
%   out = cpdem_ModelVpVs(..., dens, opts)
%
% Follows Cilli and Chapman (2021): the bulk modulus is modelled with the
% bulk modulus EPAR (alpha_K) and the shear modulus with the shear modulus
% EPAR (alpha_mu), each with cpdem_Forward, and the moduli are converted to
% velocities with a bulk density. In the paper, Vp/Vs trends for mixed wet
% sandstones are predicted in the absence of porosity and density
% measurements by combining the modelled moduli with the Gardner et al.
% (1974) velocity-density relation for sandstones; that is the default
% here ('gardner'). Measured densities ('measured') are the natural
% alternative when available. The 'implied' option, which builds a density
% from the porosity the electrical DEM passes through, is an extension of
% this toolbox and not part of the paper. Vp/Vs itself needs no density.
%
% Inputs
%   k1,mu1,sig1,k2,mu2,sig2,sig0 : as in cpdem_Forward (sig0 may be an array)
%   kepar, muepar : bulk and shear modulus EPARs (scalar, per-sample array,
%                   or function handle of porosity)
%   dens  : how to obtain density (optional; default 'gardner'):
%           struct('method','measured','rho',rho)   measured densities (kg/m^3)
%           struct('method','implied','rho1',rho1,'rho2',rho2[,'path','mean'])
%                 rho = (1 - phi) rho1 + phi rho2 with the porosity implied
%                 by the electrical DEM along the alpha_K path ('K'), the
%                 alpha_mu path ('mu') or their mean ('mean', default)
%           struct('method','gardner'[,'a',310,'b',0.25])
%                 Gardner et al. (1974): rho = a Vp^b (SI: a = 310 kg m^-3
%                 (m/s)^-b, b = 0.25); only appropriate for sandstones, and
%                 solved jointly with Vp
%   opts  : optional cpdem_Options struct
%
% Output: struct with fields (each shaped like sig0)
%   k, mu    : modelled bulk (alpha_K) and shear (alpha_mu) moduli (Pa)
%   phiK, phiMu : porosity implied along each path
%   rho      : density used (kg/m^3)
%   vp, vs   : velocities (m/s)
%   vpvs     : Vp/Vs ratio
%   M        : P-wave modulus k + 4/3 mu (Pa)
%
% References:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Gardner, G. H. F., Gardner, L. W., and Gregory, A. R. (1974). Formation
% velocity and density - The diagnostic basics for stratigraphic traps.
% Geophysics, 39(6), 770-780, DOI:10.1190/1.1440465
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 10 || isempty(dens)
    dens = struct('method', 'gardner');
end
if nargin < 11 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end

[k, ~, phiK]   = cpdem_Forward(k1, mu1, sig1, k2, mu2, sig2, sig0, kepar, opts);
[~, mu, phiMu] = cpdem_Forward(k1, mu1, sig1, k2, mu2, sig2, sig0, muepar, opts);

M = k + 4/3*mu;

switch lower(dens.method)
    case 'measured'
        rho = dens.rho;
        if numel(rho) == 1
            rho = repmat(rho, size(sig0));
        elseif numel(rho) ~= numel(sig0)
            error('CPDEM:badInput', 'dens.rho must have one density per measurement.');
        end
        rho = reshape(rho, size(sig0));
        vp = sqrt(M./rho);
        vs = sqrt(mu./rho);

    case 'implied'
        if isfield(dens, 'path')
            path = dens.path;
        else
            path = 'mean';
        end
        switch lower(path)
            case 'k',    phi = phiK;
            case 'mu',   phi = phiMu;
            case 'mean', phi = (phiK + phiMu)/2;
            otherwise
                error('CPDEM:badInput', 'dens.path must be ''K'', ''mu'' or ''mean''.');
        end
        rho = (1 - phi)*dens.rho1 + phi*dens.rho2;
        vp = sqrt(M./rho);
        vs = sqrt(mu./rho);

    case 'gardner'
        if isfield(dens, 'a'), a = dens.a; else, a = 0.31*1000; end
        if isfield(dens, 'b'), b = dens.b; else, b = 0.25; end
        % rho = a vp^b and rho vp^2 = M  =>  vp = (M/a)^(1/(2+b))
        vp  = (M./a).^(1/(2 + b));
        rho = a*vp.^b;
        vs  = sqrt(mu./rho);

    otherwise
        error('CPDEM:badInput', 'dens.method must be ''measured'', ''implied'' or ''gardner''.');
end

out = struct('k', k, 'mu', mu, 'phiK', phiK, 'phiMu', phiMu, 'rho', rho, ...
             'vp', vp, 'vs', vs, 'vpvs', vp./vs, 'M', M);

end
