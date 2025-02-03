module RationalFunctionApproximation

using LinearAlgebra, Statistics, GenericLinearAlgebra, ComplexRegions
using PyFormattedStrings

export Barycentric, Thiele, nodes, weights, degree, rewind,
    unit_interval, unit_circle, unit_disk, isclosed
include("types.jl")

export poles, residues, roots
include("barycentric.jl")
include("thiele.jl")

export aaa
include("aaa.jl")

export approximate, check
include("approximate.jl")

export minimax
include("lawson.jl")

include("operations.jl")

# These are overloaded by plotting extensions.
export convergenceplot, errorplot, animate

function animate(::Any)
    error("Load the Makie or Plots package first.")
end

function convergenceplot(::Any)
    error("Load the Makie or Plots package first.")
end

function errorplot(::Any; kwargs...)
    error("Load the Makie package first.")
end

end  # module
