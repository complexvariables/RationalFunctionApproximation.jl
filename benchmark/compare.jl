# Save/compare runner for before/after performance checks.
#   Baseline:  julia --project=benchmark benchmark/compare.jl save /tmp/bench_base.json
#   After:     julia --project=benchmark benchmark/compare.jl judge /tmp/bench_base.json
using BenchmarkTools
include("benchmarks.jl")

mode = length(ARGS) >= 1 ? ARGS[1] : "save"
path = length(ARGS) >= 2 ? ARGS[2] : "/tmp/bench_base.json"

tune!(SUITE)
results = run(SUITE; verbose = true)

if mode == "save"
    BenchmarkTools.save(path, median(results))
    println("\nSaved median results to $path")
elseif mode == "judge"
    baseline = BenchmarkTools.load(path)[1]
    j = judge(median(results), baseline)
    println("\n=== judge (current vs baseline) ===")
    show(stdout, MIME"text/plain"(), j)
    println()
    println("\n=== regressions ===")
    show(stdout, MIME"text/plain"(), regressions(j))
    println()
    println("\n=== improvements ===")
    show(stdout, MIME"text/plain"(), improvements(j))
    println()
end
