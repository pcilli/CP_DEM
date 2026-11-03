function [sig0, phi0, mOther] = cpdem_Inverse(m0, which, k1,mu1,sig1,k2,mu2,sig2,epar,opts)
% Invert for electrical conductivity from an elastic modulus (direct ODE)
%
%   sig0 = cpdem_Inverse(k0,  'K',  k1,mu1,sig1,k2,mu2,sig2,epar)
%   sig0 = cpdem_Inverse(mu0, 'mu', k1,mu1,sig1,k2,mu2,sig2,epar)
%   [sig0, phi0, mOther] = cpdem_Inverse(..., opts)
%
% The cross-property DEM of Cilli and Chapman (2021) is a one-to-one
% mapping between conductivity and each modulus, so it can be inverted
% without optimisation by integrating the same equations with the modulus
% as the independent variable, e.g. for the bulk modulus
%
%   d sig/dK = (sig2 - sig) m-bar / ((K2 - K) P),
%   d mu /dK = (mu2 - mu)   Q     / ((K2 - K) P),
%   d phi/dK = (1 - phi)          / ((K2 - K) P),
%
% from K = k1 (sig = sig1) to the measured modulus. This is exact to solver
% tolerance, typically ~100x faster than the bounded search of v1.0's
% cpdem_InvertConduct (which is retained as a wrapper), and gives all
% measurements in a single integration when the parameters are shared.
%
% Inputs
%   m0     : measured modulus (Pa); scalar or array of any shape
%   which  : 'K' to invert the bulk modulus (giving sigma_K, use alpha_K)
%            or 'mu' to invert the shear modulus (sigma_mu, use alpha_mu);
%            the v1.0 modFlag values 0 ('K') and 1 ('mu') are also accepted
%   k1,mu1,sig1, k2,mu2,sig2, epar, opts : as in cpdem_Forward
%
% Outputs (shape of m0)
%   sig0   : conductivity consistent with the measured modulus (S/m)
%   phi0   : porosity implied by the model path
%   mOther : the other modulus along the same path (mu for 'K', K for 'mu')
%
% Measurements outside the interval between the inclusion and background
% modulus cannot be matched by the model and return NaN with a warning.
%
% Reference:
% Cilli, P.A., and Chapman, M. (2021), Linking elastic and electrical
% properties of rocks using cross-property DEM. Geophysical Journal
% International, 225(3), 1812-1823, DOI:10.1093/gji/ggab046
%
% Part of the Cross-Property DEM Toolbox (CP_DEM) v2.0.

if nargin < 9
    error('CPDEM:badInput', 'cpdem_Inverse needs at least 9 inputs.');
end
if nargin < 10 || isempty(opts)
    opts = cpdem_Options();
else
    opts = cpdem_Options(opts);
end

useK = cpdem_WhichModulus(which);

cpdem_CheckPositive('k1', k1, 'mu1', mu1, 'sig1', sig1, 'k2', k2, ...
    'mu2>=0', mu2, 'sig2', sig2, 'm0>=0', m0);

shape = size(m0);
n = numel(m0);
m0 = m0(:);
[k1v, mu1v, sig1v, k2v, mu2v, sig2v] = cpdem_Broadcast(n, k1, mu1, sig1, k2, mu2, sig2);

if isnumeric(epar) && ~isscalar(epar)
    if numel(epar) ~= n
        error('CPDEM:badEPAR', ...
            'EPAR must be a scalar, a function handle, or have one value per measurement.');
    end
    eparv = epar(:);
    perSampleEPAR = true;
else
    perSampleEPAR = false;
end

if useK
    m1v = k1v; m2v = k2v;
else
    m1v = mu1v; m2v = mu2v;
end

bad = ~cpdem_InRange(m0, m1v, m2v);
if any(bad)
    warning('CPDEM:outOfRange', ...
        '%d of %d moduli lie outside the interval between the background and inclusion modulus; NaN returned for them.', ...
        nnz(bad), n);
end

sig0   = nan(n, 1);
phi0   = nan(n, 1);
mOther = nan(n, 1);

shared = all([isscalar(k1), isscalar(mu1), isscalar(sig1), isscalar(k2), ...
              isscalar(mu2), isscalar(sig2)]) && ~perSampleEPAR;

if shared
    ok = ~bad;
    if any(ok)
        [x, y0, fun] = setup(useK, k1, mu1, sig1, k2, mu2, sig2, epar);
        Y = cpdem_Integrate(fun, x(m0(ok)), y0, opts);
        sig0(ok)   = sig1*exp(Y(:, 1));
        mOther(ok) = Y(:, 2)*k1;
        phi0(ok)   = Y(:, 3);
    end
else
    for ii = find(~bad).'
        if perSampleEPAR
            e = eparv(ii);
        else
            e = epar;
        end
        [x, y0, fun] = setup(useK, k1v(ii), mu1v(ii), sig1v(ii), k2v(ii), mu2v(ii), sig2v(ii), e);
        Y = cpdem_Integrate(fun, x(m0(ii)), y0, opts);
        sig0(ii)   = sig1v(ii)*exp(Y(1));
        mOther(ii) = Y(2)*k1v(ii);
        phi0(ii)   = Y(3);
    end
end

sig0   = reshape(sig0, shape);
phi0   = reshape(phi0, shape);
mOther = reshape(mOther, shape);

end

function [x, y0, fun] = setup(useK, k1, mu1, sig1, k2, mu2, sig2, epar)
% Independent variable x = (m - m1)/k1 starts at 0; state y = [ln(sig/sig1);
% other modulus/k1; phi].
if useK
    x   = @(m) (m - k1)/k1;
    y0  = [0; mu1/k1; 0];
    fun = @(x, y) derivK(x, y, k1, k2, mu2, sig1, sig2, epar);
else
    x   = @(m) (m - mu1)/k1;
    y0  = [0; 1; 0];
    fun = @(x, y) derivMu(x, y, k1, mu1, k2, mu2, sig1, sig2, epar);
end
end

function dy = derivK(x, y, k1, k2, mu2, sig1, sig2, epar)
K   = k1 + x*k1;
sig = sig1*exp(y(1));
mu  = y(2)*k1;
phi = y(3);
alpha = cpdem_EvalEPAR(epar, phi);
[p, q] = cpdem_PQ(K, mu, k2, mu2, alpha);
mb = mbar(sig, sig2, alpha);
h = k1/((k2 - K)*p);                 % dphi/dx per unit (1 - phi)
dy = [ (sig2 - sig)*mb/sig*h;
       (mu2 - mu)*q/k1*h;
       (1 - phi)*h ];
end

function dy = derivMu(x, y, k1, mu1, k2, mu2, sig1, sig2, epar)
mu  = mu1 + x*k1;
sig = sig1*exp(y(1));
K   = y(2)*k1;
phi = y(3);
alpha = cpdem_EvalEPAR(epar, phi);
[p, q] = cpdem_PQ(K, mu, k2, mu2, alpha);
mb = mbar(sig, sig2, alpha);
h = k1/((mu2 - mu)*q);
dy = [ (sig2 - sig)*mb/sig*h;
       (k2 - K)*p/k1*h;
       (1 - phi)*h ];
end
