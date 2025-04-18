module RationalFunctionApproximation

using LinearAlgebra, Statistics, GenericLinearAlgebra, ComplexRegions, GenericSchur
using PyFormattedStrings
using Infiltrator

export unit_interval, unit_circle, unit_disk
include("utils.jl")

export nodes, weights, degree, degrees, poles, residues, roots
include("abstract.jl")

export PartialFractions, PartialFractionExpansion, pfe
include("parfrac.jl")

export Barycentric, Thiele
include("barycentric.jl")
include("thiele.jl")

export approximate, check, rewind, get_history, test_points
include("approximate.jl")

# legacy implementation of AAA
export aaa
include("aaa.jl")

export minimax
include("lawson.jl")

include("operations.jl")

# These are overloaded by plotting extensions.
export convergenceplot, errorplot, animate, poleplot

function animate(::Any)
    error("Load the Makie or Plots package first.")
end

function convergenceplot(::Any)
    error("Load the Makie or Plots package first.")
end

function errorplot(::Any; kwargs...)
    error("Load the Makie package first.")
end

function poleplot(::Any; kwargs...)
    error("Load the Makie package first.")
end

end  # module
