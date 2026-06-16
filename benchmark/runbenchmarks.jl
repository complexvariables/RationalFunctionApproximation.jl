# Convenience runner:  julia --project=benchmark benchmark/runbenchmarks.jl
using BenchmarkTools
include("benchmarks.jl")

tune!(SUITE)
results = run(SUITE; verbose = true)
show(stdout, MIME"text/plain"(), results)
println()
