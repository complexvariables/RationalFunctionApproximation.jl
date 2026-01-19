module RationalFunctionApproximation

const RFA = RationalFunctionApproximation
export RFA

using LinearAlgebra, Statistics, GenericLinearAlgebra, ComplexRegions, GenericSchur
using PyFormattedStrings

export unit_interval, unit_circle, unit_disk, DiscretizedPath
include("utils.jl")

export DiscretizedPath
include("discretized_path.jl")

export nodes, weights, degree, degrees, poles, Res, residues, roots
include("abstract-rational.jl")

export approximate, get_function, domain, check, rewind, get_history, test_points
include("approximation.jl")

export Barycentric, AAA, Thiele, TCF, derivative, evaluate
include("barycentric.jl")
include("thiele.jl")

export ArnoldiBasis, ArnoldiPolynomial, PartialFractions
include("parfrac.jl")

# legacy implementation of AAA
export aaa
include("aaa.jl")

export minimax
include("lawson.jl")

# include("operations.jl")

# These are overloaded by plotting extensions.
export convergenceplot, convergenceplot!, errorplot, errorplot!, animate, poleplot, poleplot!

# COV_EXCL_START
function animate(::Any)
    error("Load the Makie or Plots package first.")
end

"""
    convergenceplot(r)

Plot the convergence history of a rational approximation.

Markers show the maximum error on (the boundary of) the domain as a function of the numerator/denominator degree. A red marker indicates that the approximation has disallowed poles in its domain. A gold halo highlights the best approximation.
"""
function convergenceplot(::Any)
    error("Load the Makie or Plots package first.")
end

function convergenceplot!(::Any)
    error("Load the Makie or Plots package first.")
end

"""
    errorplot(r; use_abs=false)

Plot the pointwise error of an `Approximation` on (the boundary of) its domain. If the error is not real, then the real and imaginary parts are plotted separately, unless `use_abs=true`.
"""
function errorplot(::Any; kwargs...)
    error("Load the Makie or Plots package first.")
end

function errorplot!(::Any; kwargs...)
    error("Load the Makie or Plots package first.")
end

"""
    poleplot(r)

Plot the poles of the rational approximant.
"""
function poleplot(::Any; kwargs...)
    error("Load the Makie package first.")
end

function poleplot!(::Any; kwargs...)
    error("Load the Makie package first.")
end
# COV_EXCL_END

end  # module
