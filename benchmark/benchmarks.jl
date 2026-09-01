using BenchmarkTools, RationalFunctionApproximation, ComplexRegions
using Logging

# silence convergence warnings so benchmark output stays clean
global_logger(SimpleLogger(stderr, Logging.Error))

# AirspeedVelocity runs this one script against both the base revision and the PR head,
# so it has to work with either calling convention: the old API selected the interpolant
# type with a `method = Barycentric` keyword, the current one takes an instance
# (`Barycentric()`) as the last positional argument. Pick the form once, at load time, so
# no branch survives into the benchmarked call. Note the deprecation shim is deliberately
# *not* used here: it warns on every call, which would skew the timings it appears in.
if hasmethod(Barycentric, Tuple{})
    approx(f, domain, M; kw...) = approximate(f, domain, M(); kw...)
else
    approx(f, domain, M; kw...) = approximate(f, domain; method = M, kw...)
end

const SUITE = BenchmarkGroup()

# --- construction cost, by method and domain ---
SUITE["approximate"] = BenchmarkGroup()
for (name, method) in (("aaa", Barycentric), ("thiele", Thiele))
    g = SUITE["approximate"][name] = BenchmarkGroup()
    g["exp_interval"] = @benchmarkable approx(exp, $unit_interval, $method; allowed = true)
    g["tanh_steep"]   = @benchmarkable approx(x -> tanh(50x), $unit_interval, $method; allowed = true)
    g["abs_circle"]   = @benchmarkable approx(z -> abs(z - 1.0001im), $unit_circle, $method; allowed = true)
end

# --- construction cost on a discrete point set ---
# log-clustered points near 0, matching the discrete-domain test setup
let zc = 10.0 .^ range(-15, 0, 500)
    global const DISCRETE_PTS = [-reverse(zc); 0.0; zc]
end
SUITE["approximate_discrete"] = BenchmarkGroup()
for (name, method) in (("aaa", Barycentric), ("thiele", Thiele))
    g = SUITE["approximate_discrete"][name] = BenchmarkGroup()
    g["tanh_steep"]  = @benchmarkable approx(x -> tanh(100x), $DISCRETE_PTS, $method; allowed = true)
    g["abs_shift"]   = @benchmarkable approx(x -> abs(x + 0.5 + 0.01im), $DISCRETE_PTS, $method; allowed = true)
    g["sin_recip"]   = @benchmarkable approx(x -> sin(1 / (1.05 - x)), $DISCRETE_PTS, $method; allowed = true)
end

# --- evaluation cost on a fixed approximant ---
SUITE["evaluate"] = BenchmarkGroup()
let r = approx(x -> tanh(50x), unit_interval, Barycentric; allowed = true),
    z = collect(range(-1, 1, 1000))

    SUITE["evaluate"]["bary_vector"] = @benchmarkable $r.($z)
end
let r = approx(x -> tanh(50x), unit_interval, Thiele; allowed = true),
    z = collect(range(-1, 1, 1000))

    SUITE["evaluate"]["thiele_vector"] = @benchmarkable $r.($z)
end

# --- pole solve ---
SUITE["poles"] = BenchmarkGroup()
let r = approx(x -> 1 / sqrt(x^2 + 0.01), unit_interval, Thiele)
    SUITE["poles"]["thiele"] = @benchmarkable poles($r)
end
