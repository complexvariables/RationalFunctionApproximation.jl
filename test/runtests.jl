using RationalFunctionApproximation, Test, ComplexRegions, DoubleFloats, Logging
const RFA = RationalFunctionApproximation

pass(f, r, z; kw...) = isapprox(f.(z), r.(z), norm=u->maximum(abs, u); kw...)
logger = Logging.SimpleLogger(stderr, Logging.Error)
global_logger(logger)

@testset "AAA" verbose=true begin
    include("aaa.jl")
end

@testset "Discretzed path" verbose=true begin
    include("discretized-path.jl")
end

@testset "Real interval" verbose=true begin
    include("interval.jl")
end

@testset "Imaginary axis" verbose=true begin
    include("imagaxis.jl")
end

@testset "Circle and disk" verbose=true begin
    include("circle.jl")
end

@testset "Custom curve" verbose=true begin
    include("custom.jl")
end
