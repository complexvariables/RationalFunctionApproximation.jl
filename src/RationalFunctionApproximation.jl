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

# include("plotrecipes.jl")

end
