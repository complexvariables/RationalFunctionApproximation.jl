using BenchmarkTools, RationalFunctionApproximation, ComplexRegions
using Logging

# silence convergence warnings so benchmark output stays clean
global_logger(SimpleLogger(stderr, Logging.Error))

const SUITE = BenchmarkGroup()

# --- construction cost, by method and domain ---
SUITE["approximate"] = BenchmarkGroup()
for (name, method) in (("aaa", Barycentric), ("thiele", Thiele))
    g = SUITE["approximate"][name] = BenchmarkGroup()
    g["exp_interval"] = @benchmarkable approximate(exp, $unit_interval; method = $method, allowed = true)
    g["tanh_steep"]   = @benchmarkable approximate(x -> tanh(50x), $unit_interval; method = $method, allowed = true)
    g["abs_circle"]   = @benchmarkable approximate(z -> abs(z - 1.0001im), $unit_circle; method = $method, allowed = true)
end

# --- construction cost on a discrete point set ---
# log-clustered points near 0, matching the discrete-domain test setup
let zc = 10.0 .^ range(-15, 0, 500)
    global const DISCRETE_PTS = [-reverse(zc); 0.0; zc]
end
SUITE["approximate_discrete"] = BenchmarkGroup()
for (name, method) in (("aaa", Barycentric), ("thiele", Thiele))
    g = SUITE["approximate_discrete"][name] = BenchmarkGroup()
    g["tanh_steep"]  = @benchmarkable approximate(x -> tanh(100x), $DISCRETE_PTS; method = $method, allowed = true)
    g["abs_shift"]   = @benchmarkable approximate(x -> abs(x + 0.5 + 0.01im), $DISCRETE_PTS; method = $method, allowed = true)
    g["sin_recip"]   = @benchmarkable approximate(x -> sin(1 / (1.05 - x)), $DISCRETE_PTS; method = $method, allowed = true)
end

# --- evaluation cost on a fixed approximant ---
SUITE["evaluate"] = BenchmarkGroup()
let r = approximate(x -> tanh(50x), unit_interval, method = Barycentric, allowed = true),
    z = collect(range(-1, 1, 1000))

    SUITE["evaluate"]["bary_vector"] = @benchmarkable $r.($z)
end
let r = approximate(x -> tanh(50x), unit_interval; method = Thiele, allowed = true),
    z = collect(range(-1, 1, 1000))

    SUITE["evaluate"]["thiele_vector"] = @benchmarkable $r.($z)
end

# --- pole solve ---
SUITE["poles"] = BenchmarkGroup()
let r = approximate(x -> 1 / sqrt(x^2 + 0.01), unit_interval; method = Thiele)
    SUITE["poles"]["thiele"] = @benchmarkable poles($r)
end
