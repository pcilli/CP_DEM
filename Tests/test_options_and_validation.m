function test_options_and_validation()
% Options handling, input validation and out-of-range behaviour.

p = cpdem_TestParams();

% Options: defaults, struct and name/value overrides, unknown names
o = cpdem_Options();
assert(o.RelTol == 1e-10 && o.AbsTol == 1e-13 && strcmp(o.Solver, 'ode45') && isequal(o.Bounds, [1 100]), 'defaults');
o = cpdem_Options(struct('RelTol', 1e-6), 'norm', 1);
assert(o.RelTol == 1e-6 && o.Norm == 1 && o.AbsTol == 1e-13, 'overrides (case-insensitive)');
threw = false;
try
    cpdem_Options('NoSuchOption', 1);
catch
    threw = true;
end
assert(threw, 'unknown option must throw');

% Looser tolerances still agree with the defaults to the looser level
[k0, mu0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 16.4);
[kl, ml] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, p.sig0, 16.4, ...
    struct('RelTol', 1e-5, 'AbsTol', 1e-8));
cpdem_AssertClose([kl ml], [k0 mu0], 1e-4, 0, 'loose tolerance');

% Invalid inputs
bad = { {p.k1, p.mu1, 0, p.k2, p.mu2, p.sig2, 0.1, 16.4},      ... insulating host
        {p.k1, p.mu1, p.sig1, p.k2, -1, p.sig2, 0.1, 16.4},    ... negative mu2
        {p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 0.1, 0},    ... EPAR 0
        {p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, 0.1, [1 2]},... wrong EPAR length
        {[1 2 3]*p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, [0.1 0.2], 16.4} }; % wrong param length
for ii = 1:numel(bad)
    threw = false;
    try
        cpdem_Forward(bad{ii}{:});
    catch
        threw = true;
    end
    assert(threw, sprintf('invalid input set %d must throw', ii));
end

% Out of range conductivities give NaN (with a warning) and do not disturb
% the others
ws = warning('off', 'CPDEM:outOfRange');
[k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, [0.1 6 1e-5 0.2], 16.4);
warning(ws);
assert(all(isnan([k0(2:3) mu0(2:3) phi0(2:3)])), 'NaN out of range');
[kr, mr] = cpdem_Forward(p.k1, p.mu1, p.sig1, p.k2, p.mu2, p.sig2, [0.1 0.2], 16.4);
cpdem_AssertClose([k0([1 4]) mu0([1 4])], [kr mr], 1e-12, 0, 'in-range samples unaffected');

% Reverse contrast (conductive host, resistive inclusions) integrates too
[k0, mu0, phi0] = cpdem_Forward(p.k1, p.mu1, 5, p.k2, p.mu2, 1e-4, [1 0.1], 1);
assert(all(isfinite([k0 mu0 phi0])) && all(phi0 > 0) && all(k0 < p.k1), 'reverse contrast');

end
