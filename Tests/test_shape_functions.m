function test_shape_functions()
% theta, f and L: sphere limits, continuity of the series switch, identity
% L = 1 - theta, and symmetry of conventions.

% Sphere limits
[t, f, L] = cpdem_Shape(1);
cpdem_AssertClose([t f L], [2/3 -2/5 1/3], 1e-14, 0, 'sphere limits');

% Series (|e| < 0.1) and closed form (|e| >= 0.1) agree across the switch
for e = [-0.1 0.1]
    [t1, f1, L1] = cpdem_Shape(1 + e*(1 - 1e-9));   % series side
    [t2, f2, L2] = cpdem_Shape(1 + e*(1 + 1e-9));   % closed-form side
    cpdem_AssertClose([t1 f1 L1], [t2 f2 L2], 1e-9, 0, sprintf('continuity at e = %g', e));
end
% Leading terms of the series
for e = [-1e-3 1e-3]
    [t, f, L] = cpdem_Shape(1 + e);
    ts = 2/3 + 4*e/15 - 6*e^2/35 + 32*e^3/315 - 40*e^4/693;
    fs = -2/5 - 12*e/35 + 2*e^2/15 - 32*e^3/1155 - 40*e^4/3003;
    cpdem_AssertClose([t f L], [ts fs 1 - ts], 1e-12, 0, sprintf('series at e = %g', e));
end
% The series is smooth through the sphere: central difference of theta
% equals 4/15 to high accuracy
h = 1e-4;
[tp] = cpdem_Shape(1 + h); [tm] = cpdem_Shape(1 - h);
cpdem_AssertClose((tp - tm)/(2*h), 4/15, 1e-7, 0, 'dtheta/dEPAR at the sphere');

% L == 1 - theta over a wide range of aspect ratios
for a = [1e-3 0.01 0.1 0.5 0.9 1.1 2 10 16.4 100 1e3]
    [t, ~, L] = cpdem_Shape(a);
    cpdem_AssertClose(L, 1 - t, 1e-10, 0, sprintf('L = 1 - theta at %g', a));
end

% Known depolarisation limits: thin disk L -> 1, long needle L -> 0
[~, ~, Ld] = cpdem_Shape(1e-4);
[~, ~, Ln] = cpdem_Shape(1e4);
assert(Ld > 0.999 && Ln < 1e-6, 'disk/needle depolarisation limits');

% Invalid aspect ratios are rejected
for bad = {0, -1, Inf, NaN, [1 2]}
    threw = false;
    try
        cpdem_Shape(bad{1});
    catch
        threw = true;
    end
    assert(threw, 'cpdem_Shape must reject invalid EPAR');
end

end
