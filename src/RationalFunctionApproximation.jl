module RationalFunctionApproximation

const RFA = RationalFunctionApproximation
export RFA

using LinearAlgebra, Statistics, GenericLinearAlgebra, ComplexRegions, GenericSchur, ArnoldiVandermonde
using PyFormattedStrings
using PrecompileTools
import ArnoldiVandermonde: degree

export unit_interval, unit_circle, unit_disk, DiscretizedPath
include("utils.jl")

export DiscretizedPath
include("discretized_path.jl")

export nodes, weights, degree, degrees, poles, Res, residues, roots
include("abstract-rational.jl")

export approximate, get_function, domain, check, rewind, get_history, test_points
include("approximation.jl")

export BarycentricInterpolant, Barycentric, AAA, derivative, evaluate
include("barycentric.jl")

export ContinuedFractionInterpolant, Thiele, TCF
include("thiele.jl")

export ArnoldiBasis, ArnoldiPolynomial, PartialFractionExpansion, PartialFractions
include("parfrac.jl")

# legacy implementation of AAA
export aaa
include("aaa.jl")

export minimax
include("lawson.jl")

@setup_workload begin
    x_interval = range(-1, 1, 200)
    @compile_workload begin
        for tag in (BarycentricInterpolant(), ContinuedFractionInterpolant())
            approximate(sin, unit_circle, tag)
            for domain in (unit_interval, x_interval)
                approximate(sin, domain, tag)
                approximate(cis, domain, tag)
                approximate(x -> 1/(x^2 + 4), domain, [2im, -2im], PartialFractionExpansion())
            end
        end
    end
end



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
