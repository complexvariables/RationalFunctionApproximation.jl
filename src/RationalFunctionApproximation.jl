module RationalFunctionApproximation

using LinearAlgebra, ComplexRegions
using PyFormattedStrings

export Barycentric, nodes, weights, degree, rewind,
    unit_interval, unit_circle, unit_disk, isclosed
include("types.jl")

export poles, residues, roots
include("methods.jl")

export aaa
include("aaa.jl")

export approximate, check
include("approximate.jl")

export minimax
include("lawson.jl")

include("operations.jl")

# These are overloaded by plotting extensions.
export convergenceplot, errorplot

function convergenceplot(r::Union{Barycentric,Approximation})
    @error "Load the Makie package to make a convergence plot"
end

function errorplot(r::Approximation; args...)
    @error "Load the Makie package to make an error plot"
end

end
