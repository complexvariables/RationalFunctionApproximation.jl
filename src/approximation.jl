#####
##### TYPES
#####

mutable struct IterationRecord{R,S,T}
    interpolant::R
    error::S
    poles::Union{Missing, Vector{T}}

    function IterationRecord(r::R, error, poles) where
        {S<:AbstractFloat, R<:AbstractRationalInterpolant{S}}
        return new{R,S,Complex{S}}(copy(r), error, poles)
    end
end

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", h::IterationRecord)
    print(IOContext(io, :compact=>true), "$(degrees(h.interpolant)) rational interpolant with error $(round(h.error, sigdigits=3))")
end
# COV_EXCL_STOP

"""
    Approximation (type)

Approximation of a function on a domain.

# Fields
- `original`: the original function
- `domain`: the domain of the approximation
- `fun`: the barycentric representation of the approximation
- `allowed`: function to determine if a pole is allowed
- `path`: a `DiscretizedPath` for the domain boundary
- `history`: all approximations in the iteration
"""
struct Approximation{T,S} <: Function
    original::Function
    domain::Domain{T}
    fun::Union{AbstractRationalInterpolant{T,S},AbstractRationalFunction{S}}
    allowed::Function
    path::DiscretizedPath
    history::Union{Vector{<:IterationRecord},Nothing}
end

function Approximation(
    f::Function,
    domain::Domain,
    fun::AbstractRationalInterpolant,
    allowed::Function,
    path::DiscretizedPath
    )
    return Approximation(f, domain, fun, allowed, path, nothing)
end

(f::Approximation)(z) = f.fun(z)

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", f::Approximation)
    print(io, f.fun)
    print(IOContext(io, :compact=>true), " on the domain: ", f.domain)
end
# COV_EXCL_STOP

nodes(r::Approximation, args...) = nodes(r.fun, args...)
Base.values(r::Approximation, args...) = values(r.fun, args...)
weights(r::Approximation, args...) = weights(r.fun, args...)
degree(r::Approximation) = degree(r.fun)
degrees(r::Approximation) = degrees(r.fun)
poles(F::Approximation) = poles(F.fun)
residues(f::Approximation, args...) = residues(f.fun, args...)
roots(f::Approximation) = roots(f.fun)

function test_points(r::Approximation; with_parameters=false)
    s, z = collect(r.path, :all)
    return with_parameters ? (s, z) : z
end

#####
##### Documentation strings
#####

@doc """
    approximate(f, domain)

Adaptively compute a rational interpolant on a continuous or discrete domain.

# Arguments
## Continuous domain
- `f::Function`: function to approximate
- `domain`: curve, path, or region from ComplexRegions

## Discrete domain
- `f::Function`: function to approximate
- `z::AbstractVector`: point set on which to approximate

# Keywords
- `max_iter::Integer=150`: maximum number of iterations on node addition
- `float_type::Type`: floating point type to use for the computation¹
- `tol::Real=1000*eps(float_type)`: relative tolerance for stopping
- `allowed::Function`: function to determine if a pole is allowed²
- `refinement::Integer=3`: number of test points between adjacent nodes (continuum only)
- `stagnation::Integer=5`: number of iterations to determine stagnation

¹Default of `float_type` is the promotion of `float(1)` and the float type of the domain.
²Default is to disallow poles on the curve or in the interior of a continuous domain, or to accept all poles on a discrete domain. Use `allowed=true` to allow all poles.

# Returns
- `r::Approximation`: the rational interpolant

See also [`Approximation`](@ref), [`check`](@ref), [`rewind`](@ref).

# Examples
```julia-repl
julia> f = x -> tanh( 40*(x - 0.15) );

julia> r = approximate(f, unit_interval)
Barycentric{Float64, Float64} rational function of type (22, 22) on the domain: Path{Float64} with 1 curve

julia> ( r(0.3), f(0.3) )
(0.9999877116508015, 0.9999877116507956)

julia> check(r);   # accuracy over the domain
[ Info: Max error is 1.58e-13
```
"""
approximate

@doc """
    approximate(f, domain, poles)

Computes a linear least-squares approximation with prescribed poles.

# Arguments
## Continuous domain
- `f::Function`: function to approximate
- `domain`: curve, path, or region from ComplexRegions

## Discrete domain
- `f::Function`: function to approximate
- `z::AbstractVector`: point set on which to approximate

# Keywords
- `degree`: degree of the polynomial part of the approximant (defaults to length(poles) ÷ 2)
- `init`: initial number of nodes on the path (continuum only)
- `refinement`: number of test points between adjacent nodes (continuum only)

# Returns
- `r::Approximation`: the rational approximant

See also [`Approximation`](@ref), [`check`](@ref), [`rewind`](@ref).

# Examples
```julia-repl
julia> f = tanh;

julia> ζ = 1im*[-π/2, π/2];

julia> r = approximate(f, unit_interval, ζ; degree=10)
PartialFractions{ComplexF64} rational function of type (12, 2) on the domain: Segment(-1.0,1.0)

julia> ( r(0.3), f(0.3) )
(0.2913126124509021 + 1.1102230246251565e-16im, 0.2913126124515909)

julia> check(r);   # accuracy over the domain
[ Info: Max error is 2.75e-12
```
"""
approximate

#####
##### Dispatch
#####

# Each rational type implements methods based on dispatch of Type for the first argument.
# Here, we can call those based on a method= keyword argument instead, giving a default.
function approximate(f::Function, domain::ComplexCurveOrPath; method::Type=Barycentric, kw...)
     approximate(method, f, domain; kw...)
end

function approximate(y::AbstractVector, z::AbstractVector; method::Type=Barycentric, kw...)
     approximate(method, y, z; kw...)
end

function approximate(
    f::Function, domain::ComplexCurveOrPath, ζ::AbstractVector;
    method::Type=PartialFractions,
    kw...
    )
    approximate(method, f, domain, ζ; kw...)
end

# Each rational type recognizes two signatures:
# - f::Function, domain::ComplexCurveOrPath, [poles::AbstractVector]
# - values::AbstractVector, test_points::AbstractVector, [poles::AbstractVector]

# We fill in other convenience cases here.

# ::Function, ::AbstractRegion
# Given a region as domain, we interpret poles as not being allowed in that region.
function approximate(f::Function, R::ComplexRegions.AbstractRegion; kw...)
    r = approximate(f, R.boundary; allowed=z->!in(z,R), kw...)
    return Approximation(f, R, r.fun, r.allowed, r.path, r.history)
end

# ::Function, ::AbstractVector
# Evaluate the function to call a fully discrete approximation.
function approximate(
    f::Function, z::AbstractVector;
    allowed = z -> true,
    kw...
    )
    r, history = approximate(f.(z), z; kw...)
    return Approximation(f, z, r, allowed, DiscretizedPath(), history)
end

# ::Function,::AbstractVector, ::AbstractVector
function approximate(
    f::Function, z::AbstractVector, ζ::AbstractVector;
    allowed = z -> true,
    method = PartialFractions,
    kw...
    )
    r = approximate(method, f.(z), z, ζ; kw...)
    return Approximation(f, z, r, allowed, DiscretizedPath(), nothing)
end

#####
##### Support functions
#####

"""
    rewind(r, index)

Rewind a rational approximation to a state encountered during an iteration.

# Arguments
- `r::Approximation`: the approximation to rewind
- `index::Integer`: the iteration number to rewind to

# Returns
- the rational function of the specified index (same type as input)

# Examples
```jldoctest
julia> r = approximate(x -> cos(20x), unit_interval)
Barycentric{Float64, Float64} rational interpolant of type (24, 24) on the domain: Path{Float64} with 1 curve

julia> rewind(r, 10)
Barycentric{Float64, Float64} rational interpolant of type (10, 10) on the domain: Path{Float64} with 1 curve
```
"""
function rewind(r::Approximation, idx::Integer)
    if isnothing(r.history)
        @error("No convergence history exists.")
    end
    return Approximation(r.original, r.domain, r.history[idx].interpolant, r.allowed, r.path, r.history)
end

"""
    check(r; quiet=false, prenodes=false)

Check the accuracy of a rational approximation `r` on its domain. Returns the test points and the error at those points.

# Arguments
- `r::Approximation`: rational approximation

# Keywords
- `quiet::Bool=false`: suppress @info output
- `prenodes::Bool=false`: return prenodes of the approximation as well

# Returns
- `τ::Vector`: test points
- `err::Vector`: error at test points

See also [`approximate`](@ref).
"""
function check(
    F::Approximation;
    quiet=false,
    prenodes=false,
    refinement=10
    )
    if F.domain isa AbstractVector    # discrete domain
        τ = F.domain
        t = collect(eachindex(τ))
    elseif refinement == :test
        t, τ = collect(F.path, :test)
    else
        q = DiscretizedPath(F.path, refinement)
        t, τ = collect(q, :all)
    end
    err = F.original.(τ) - F.fun.(τ)
    !quiet && @info f"Max error is {norm(err, Inf):.2e}"
    return prenodes ? (t, τ, err) : (τ, err)
end

"""
    get_history(r::Approximation)

Parse the convergence history of a rational approximation.
# Arguments
- `r::Approximation`: the approximation to get the history from
# Returns
- `::Vector`: degrees of the approximations
- `::Vector`: estimated maximum errors of the approximations
- `::Vector{Vector}`: poles of the approximations
- `::Vector{Vector}`: allowed poles of the approximations
- `::Integer`: index of the best approximation

See also [`convergenceplot`](@ref).
"""
function get_history(r::Approximation{T,S}) where {T,S}
    hist = r.history
    deg = Int[]
    zp = Vector{complex(S)}[]
    err = Float64[]
    allowed = BitVector[]
    best = 0
    for (idx, record) in enumerate(hist)
        fun = record.interpolant
        push!(deg, degree(fun))
        if ismissing(record.poles)
            record.poles = poles(fun)
        end
        push!(zp, record.poles)
        push!(allowed, r.allowed.(record.poles))
        push!(err, record.error)
        if length(nodes(r)) == length(nodes(fun))
            best = idx
        end
    end
    return deg, err, zp, allowed, best
end

# Return values for quitting_check:
#  -1: success
#   0: continue
#   n: iteration number to stop at
function quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
    # If allowed === true, do not check for allowed poles.
    n = length(history)
    err = [h.error for h in history]

    # Check for convergence
    status = 0
    if (err[end] <= tol*fmax)
        if (allowed === true)
            status = -1
        else
            zp = history[end].poles = poles(history[end].interpolant)
            status = all(allowed, zp) ? -1 : 0
        end
    end

    # Check for stagnation
    stagnant = false
    if n >= stagnation + 5
        old = n-stagnation-4:n-stagnation
        plateau = exp(median(log(err[i]) for i in old))
        stagnant = all(plateau < e < fmax/100 for e in last(err, stagnation))
    end

    # Decide on unsuccessful stopping
    if (n == max_iter) || stagnant
        if !(allowed === true)
            # Look for the best acceptable approximation:
            for k in sortperm(err)
                history[k].poles = @coalesce history[k].poles poles(history[k].interpolant)
                if all(allowed, history[k].poles)
                    n = k
                    break
                end
            end
        end
        status = n
    end

    return status
end
