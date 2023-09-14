using RationalFunctionApproximation, Test, LinearAlgebra, ComplexRegions

# ii = 1im * range(-300,400,500)
# zz = cispi.(2 * (0:500) / 500)

pass(f, r, z; kw...) = isapprox(f.(z), r.(z), norm=u->norm(u,Inf); kw...)

@testset "AAA" verbose=true begin
    include("aaa.jl")
end

@testset "Interval" verbose=true begin
    include("interval.jl")
end

@testset "Circle and disk" verbose=true begin
    include("circle.jl")
end
