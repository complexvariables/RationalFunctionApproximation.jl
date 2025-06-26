using RationalFunctionApproximation, Test, ComplexRegions, DoubleFloats, Logging
const RFA = RationalFunctionApproximation

pass(f, r, z; kw...) = isapprox(f.(z), r.(z), norm=u->maximum(abs, u); kw...)
logger = Logging.SimpleLogger(stderr, Logging.Error)
global_logger(logger)

@testset "AAA" verbose=true begin
    include("aaa.jl")
end

@testset "Discrete domains" verbose=true begin
    include("discrete.jl")
end

@testset "Discretzed path" verbose=true begin
    include("discretized-path.jl")
end

@testset "Real interval" verbose=true begin
    include("interval.jl")
end

@testset "Operations on functions" verbose=true begin
    include("operations.jl")
end

@testset "Imaginary interval" verbose=true begin
    include("imaginterval.jl")
end

@testset "Circle and disk" verbose=true begin
    include("circle.jl")
end

@testset "Custom curve" verbose=true begin
    include("custom.jl")
end

@testset "Prescribed poles" verbose=true begin
    include("parfrac.jl")
end
