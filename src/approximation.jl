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
    print(IOContext(io, :compact=>true), "$(degrees(h.interpolant)) rational approximation with error $(round(h.error, sigdigits=3))")
end
# COV_EXCL_STOP

abstract type AbstractApproximation{T,S} <: Function end

"""
    ContinuumApproximation (type)

Approximation of a function on a domain.

# Fields
- `original`: the original function
- `domain`: the domain of the approximation
- `fun`: the barycentric representation of the approximation
- `allowed`: function to determine if a pole is allowed
- `path`: a `DiscretizedPath` for the domain boundary
- `history`: all approximations in the iteration
"""
struct ContinuumApproximation{T,S,R} <: AbstractApproximation{T,S}
    original::Function
    domain::Domain{T}
    fun::R
    allowed::Union{Bool,Function}
    path::DiscretizedPath
    history::Union{Vector{<:IterationRecord},Nothing}
end

function ContinuumApproximation(
    f::Function,
    domain::Domain{T},
    fun::R,
    allowed::Union{Bool,Function},
    path::DiscretizedPath,
    history=nothing
    ) where {T,S,R<:AbstractRationalFunction{S}}
    return ContinuumApproximation{T,S,R}(f, domain, fun, allowed, path, history)
end

(f::ContinuumApproximation)(z) = f.fun(z)
domain(r::ContinuumApproximation) = r.domain
get_function(r::ContinuumApproximation) = r.fun
history(r::ContinuumApproximation) = r.history

"""
    DiscreteApproximation (type)

Approximation of a function on a domain.

# Fields
- `data`: sample values to be approximated
- `domain`: points for the sample values
- `fun`: the barycentric representation of the approximation
- `test_index`: indicator of which domain points remain as test points
- `allowed`: function to determine if a pole is allowed
- `history`: all approximations in the iteration
"""
struct DiscreteApproximation{T,S,R} <: AbstractApproximation{T,S}
    data::Vector{S}
    domain::Vector{T}
    fun::R
    test_index::BitVector
    allowed::Union{Bool,Function}
    history::Union{Vector{<:IterationRecord},Nothing}
    function DiscreteApproximation{T,S,R}(
        data::AbstractVector{S},
        domain::AbstractVector{T},
        args...
        ) where {T,S,R}
        @assert length(data) == length(domain)
        return new{T,S,R}(data, domain, args...)
    end
end

function DiscreteApproximation(
    data::AbstractVector{S},
    domain::AbstractVector{T},
    fun::R,
    test_index::BitVector,
    allowed::Union{Bool,Function}=true,
    history=nothing
    ) where {T,S,R<:AbstractRationalFunction}
    return DiscreteApproximation{T,float(S),typeof(fun)}(float(data), domain, fun, test_index, allowed, history)
end

(f::DiscreteApproximation)(z) = f.fun(z)
domain(r::DiscreteApproximation) = r.domain
get_function(r::DiscreteApproximation) = r.fun
history(r::DiscreteApproximation) = r.history

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", f::ContinuumApproximation)
    print(io, f.fun)
    print(IOContext(io, :compact=>true), " constructed on: ", f.domain)
end

function Base.show(io::IO, ::MIME"text/plain", f::AbstractApproximation)
    print(io, f.fun)
    print(IOContext(io, :compact=>true), " constructed from $(length(f.domain)) samples")
end
# COV_EXCL_STOP

nodes(r::AbstractApproximation, args...) = nodes(r.fun, args...)
Base.values(r::AbstractApproximation, args...) = values(r.fun, args...)
weights(r::AbstractApproximation, args...) = weights(r.fun, args...)
degree(r::AbstractApproximation) = degree(r.fun)
degrees(r::AbstractApproximation) = degrees(r.fun)
poles(F::AbstractApproximation) = poles(F.fun)
residues(f::AbstractApproximation, args...) = residues(f.fun, args...)
roots(f::AbstractApproximation) = roots(f.fun)

function test_points(r::ContinuumApproximation; with_parameters=false)
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
- `f::Function` or `y::AbstractVector`: function or discrete values to approximate
- `z::AbstractVector`: domain point set

# Keywords
- `method::Type`: type of rational interpolant to use (`AAA` default, `TCF`, `PartialFractions`)
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

See also [`ContinuumApproximation`](@ref), [`DiscreteApproximation`](@ref), [`check`](@ref), [`rewind`](@ref).

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
approximate(f::Function, domain)

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

See also [`ContinuumApproximation`](@ref), [`DiscreteApproximation`](@ref), [`check`](@ref), [`rewind`](@ref).

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
approximate(f::Function, domain, ζ::AbstractVector)

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
    return ContinuumApproximation(f, R, r.fun, r.allowed, r.path, r.history)
end

# ::Function, ::AbstractVector
# Evaluate the function to call a fully discrete approximation.
function approximate(
    f::Function, z::AbstractVector;
    allowed = true,
    kw...
    )
    y = f.(z)
    r = approximate(y, z; allowed, kw...)
    return DiscreteApproximation(y, z, r.fun, r.test_index, r.allowed, r.history)
end

# ::Function,::AbstractVector, ::AbstractVector
function approximate(
    f::Function, z::AbstractVector, ζ::AbstractVector;
    method = PartialFractions,
    kw...
    )
    y = f.(z)
    r = approximate(method, y, z, ζ; kw...)
    return DiscreteApproximation(y, z, r.fun, r.test_index, true, r.history)
end

#####
##### Support functions
#####

function Base.isapprox(r::AbstractApproximation, s::AbstractApproximation; kwargs...)
    return if isapprox(domain(r), domain(s))
        all( isapprox(r(z), s(z); kwargs...) for z in test_points(r) )
    else
        false
    end
end

function Base.isapprox(r::AbstractApproximation, f::Function; kwargs...)
    return all(isapprox(r(z), f(z); kwargs...) for z in test_points(r))
end

function Base.isapprox(r::AbstractApproximation, u::Number; kwargs...)
    return all(isapprox(r(z), u; kwargs...) for z in test_points(r))
end

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
function rewind(r::AbstractApproximation, idx::Integer)
    if isnothing(r.history)
        @error("No convergence history exists.")
    end
    return typeof(r)(r.original, r.domain, r.history[idx].interpolant, r.allowed, r.path, r.history)
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
function check(F::ContinuumApproximation;
    quiet=false,
    verbose=!quiet,
    prenodes=false,
    refinement=10
    )
    if refinement == :test
        t, τ = collect(F.path, :test)
    else
        q = DiscretizedPath(F.path, refinement)
        t, τ = collect(q, :all)
    end
    err = F.original.(τ) - F.(τ)
    verbose && @info f"Max error is {norm(err, Inf):.2e}"
    return prenodes ? (t, τ, err) : (τ, err)
end


function check(F::DiscreteApproximation;
    quiet=false,
    verbose=!quiet,
    )
    τ = F.domain[F.test_index]
    err = F.data[F.test_index] - F.(τ)
    verbose && @info f"Max error is {norm(err, Inf):.2e}"
    return τ, err
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
function get_history(r::AbstractApproximation{T,S}; get_poles=!(r.allowed == true)) where {T,S}
    hist = r.history
    deg = Int[]
    zp = Vector{complex(S)}[]
    err = Float64[]
    allowed = BitVector[]
    best = 0
    for (idx, record) in enumerate(hist)
        fun = record.interpolant
        push!(deg, degree(fun))
        res = []
        if get_poles && ismissing(record.poles)
            record.poles, res = residues(fun)
        end
        push!(zp, coalesce(record.poles, []))
        if !(r.allowed == true)
            allow = [r.allowed(z) || abs(R) < eps(S) for (z, R) in zip(zp[end], res)]
            push!(allowed, allow)
        else
            push!(allowed, fill(true, length(zp[end])))
        end
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
    n = length(history)
    err = [h.error for h in history]

    # Check for convergence
    # If allowed === true, do not check for allowed poles
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
    min_err, min_k = findmin(err)
    if (n >= stagnation + 5) && min_err < fmax/100
        plateau = 5min_err
        stagnant = all(plateau < e for e in last(err, stagnation)) || (min_k < n - 2stagnation)
    end

    # Decide on unsuccessful stopping
    if (n >= max_iter) || stagnant
        # Look for the best acceptable approximation:
        if (allowed === true)
            n = argmin(err)
        else
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

#####
##### Operations
#####

"""
    derivative(r::Approximation)

Create an approximation of the derivative of `r` on the same domain.
"""
function derivative(r::AbstractApproximation, order=1; kwargs...)
    # TODO: This ought to be handled by dispatch on a type parameter.
    return if isa(get_function(r), AbstractRationalInterpolant)
        approximate(derivative(get_function(r), order), domain(r); method=typeof(get_function(r)), kwargs...)
    else
        @error("Not supported. Take the derivative of the `.fun` field.")
    end
end

# Arithmetic Big 4

function Base.:+(r::AbstractApproximation, s::AbstractApproximation)
    if !(domain(r) ≈ domain(s))
        throw(DomainError("Approximation domains do not match"))
    end
    return r + s.fun
end

function Base.:+(r::AbstractApproximation, g::Function)
    f(z) = r(z) + g(z)
    return approximate(f, domain(r); method=typeof(get_function(r)))
end

function Base.:+(r::ContinuumApproximation, s::Number)
    rs = get_function(r) + s
    return ContinuumApproximation(rs, domain(r), rs, r.allowed, r.path, r.history)
end

function Base.:+(r::DiscreteApproximation, s::Number)
    rs = get_function(r) + s
    return DiscreteApproximation(r.data .+ s, domain(r), rs, r.test_index, r.allowed, r.history)
end

Base.:+(s::Union{Function,Number}, r::AbstractApproximation) = r + s

Base.:-(r::AbstractApproximation, s::AbstractApproximation) = -s + r
Base.:-(r::AbstractApproximation, s::Number) = -s + r
Base.:-(r::AbstractApproximation, s::Function) = r + ∘(-, s)
Base.:-(r::Union{Function,Number}, s::AbstractApproximation) = -s + r

# unary -
function Base.:-(r::ContinuumApproximation)
    rs = -get_function(r)
    return ContinuumApproximation(rs, domain(r), rs, r.allowed, r.path, r.history)
end

function Base.:-(r::DiscreteApproximation)
    rs = -get_function(r)
    return DiscreteApproximation(-r.data, domain(r), rs, r.test_index, r.allowed, r.history)
end

# * and / with 3 levels of generality
function Base.:*(r::AbstractApproximation, s::AbstractApproximation)
    if !(domain(r) ≈ domain(s))
        throw(DomainError("Approximation domains do not match"))
    end
    return r * s.fun
end

function Base.:*(r::AbstractApproximation, g::Function)
    f(z) = r(z) * g(z)
    return approximate(f, domain(r); method=typeof(get_function(r)))
end

function Base.:*(r::ContinuumApproximation, s::Number)
    rs = get_function(r) * s
    return ContinuumApproximation(rs, domain(r), rs, r.allowed, r.path, r.history)
end

function Base.:*(r::DiscreteApproximation, s::Number)
    rs = get_function(r) * s
    return DiscreteApproximation(r.data * s, domain(r), rs, r.test_index, r.allowed, r.history)
end

Base.:*(s::Union{Function,Number}, r::AbstractApproximation) = r * s

function Base.:/(r::AbstractApproximation, s::AbstractApproximation)
    if !(domain(r) ≈ domain(s))
        throw(DomainError("Approximation domains do not match"))
    end
    return r / s.fun
end

function Base.:/(r::AbstractApproximation, g::Function)
    f(z) = r(z) / g(z)
    return approximate(f, domain(r); method=typeof(get_function(r)))
end

function Base.:/(r::Function, s::AbstractApproximation)
    f(z) = r(z) / s.fun(z)
    return approximate(f, domain(s); method=typeof(s.fun))
end

Base.:/(r::AbstractApproximation, s::Number) = iszero(s) ? throw(DomainError("Division by zero")) : r * (1 / s)
Base.:/(r::Number, s::AbstractApproximation) = (z -> r) / s

# composition
function Base.:∘(f::Function, g::AbstractApproximation)
    # No domain checking is attempted.
    return approximate(f ∘ g.fun, g.domain; method=typeof(g.fun))
end
