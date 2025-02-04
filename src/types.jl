##
const ComplexPath = ComplexRegions.AbstractPath
const ComplexCurve = ComplexRegions.AbstractCurve
const ComplexClosedPath = ComplexRegions.AbstractClosedPath
const ComplexClosedCurve = ComplexRegions.AbstractClosedCurve
const ComplexSCRegion = ComplexRegions.AbstractSimplyConnectedRegion
isclosed(p::ComplexCurve) = isa(p, ComplexClosedCurve)
isclosed(p::ComplexPath) = isa(p, ComplexClosedPath)
const RealComplex{T} = Union{T, ComplexValues.AnyComplex{T}}
const VectorRealComplex{T} = Union{Vector{T}, Vector{Complex{T}}}
const VectorVectorRealComplex{T} = Union{Vector{Vector{T}},Vector{Vector{Complex{T}}}}

#####
##### abstract rational interpolant
#####

# types are T = float type, S = T or Complex{T}
abstract type AbstractRationalInterpolant{T,S} <: Function end

(r::AbstractRationalInterpolant)(z) = evaluate(r, z)
"nodes(r) returns a vector of the interpolation nodes of the rational interpolant."
nodes(::AbstractRationalInterpolant) = error("`nodes` not implemented for $(typeof(r))")
"values(r) returns a vector of the nodal values of the rational interpolant `r`."
Base.values(::AbstractRationalInterpolant) = error("`values` not implemented for $(typeof(r))")
weights(::AbstractRationalInterpolant) = error("`weights` not implemented for $(typeof(r))")
Base.eltype(r::AbstractRationalInterpolant) = eltype(values(r))
Base.length(r::AbstractRationalInterpolant) = length(nodes(r))
Base.isempty(r::AbstractRationalInterpolant) = isempty(nodes(r))

"degrees(r) returns the degrees of the numerator and denominator of the rational `r`."
degrees(r::AbstractRationalInterpolant) = error("`degrees` not implemented for $(typeof(r))")
"degree(r) returns the degree of the denominator of the rational `r`."
degree(r::AbstractRationalInterpolant) = error("`degree` not implemented for $(typeof(r))")

"poles(r) returns the poles of the rational interpolant `r`."
poles(::AbstractRationalInterpolant) = error("`poles` not implemented for $(typeof(r))")
"residues(r, z) returns the residues of the rational interpolant `r` at the poles in `z`."
residues(::AbstractRationalInterpolant, z) = error("`residues` not implemented for $(typeof(r))")
"roots(r) returns the roots of the rational interpolant `r`."
roots(::AbstractRationalInterpolant) = error("`roots` not implemented for $(typeof(r))")

"""
    decompose(r)

Return the roots, poles, and residues of the rational interpolant `r`.
"""
function decompose(r::AbstractRationalInterpolant)
    p = poles(r)
    return roots(r), p, residues(r, p)
end

function Base.show(io::IO, mimetype::MIME"text/plain", r::AbstractRationalInterpolant)
    ioc = IOContext(io,:compact=>get(io, :compact, true))
    len = length(r)
    if len==0
        println(ioc, "Empty $(typeof(r)) rational interpolant")
        return
    end
    println(ioc, "$(typeof(r)) rational interpolant with $len nodes and values:")
    # print out 3 nodes=>values
    nv, rest = Iterators.peel( zip(nodes(r), values(r)) )
    print(ioc, "    ", Pair(nv...))
    rest = collect(rest)
    next2 = Iterators.take(rest, 2)
    foreach(s->print(ioc, ",  ", Pair(s...)), next2)
    # if any left, just elide to the last one
    if length(rest) > 2
        print(ioc, ",  …  ")
        print(ioc, Pair(last(rest)...))
    end
end

function Base.show(io::IO, r::AbstractRationalInterpolant)
    m = length(r)
    if m==0
        println(io, "Empty $(typeof(r)) rational interpolant")
        return
    end
    print(
        IOContext(io,:compact=>true),
        "$(typeof(r)) rational interpolant on $m nodes"
        )
end

#####
##### convergence statistics
#####

"""
    ConvergenceStats{T}(bestidx, error, nbad, nodes, values, weights, poles)

Convergence statistics for a sequence of rational approximations.

# Fields
- `bestidx`: the index of the best approximation
- `error`: the error of each approximation
- `nodes`: the nodes of each approximation
- `values`: the values of each approximation
- `weights`: the weights of each approximation
- `poles`: the poles of each approximation
- `isbad`: vector indicating which poles are in the domain

See also: [`approximate`](@ref), [`Barycentric`](@ref)
"""
struct ConvergenceStats{T}
    bestidx::Int
    error::Vector{<:AbstractFloat}
    nodes::VectorVectorRealComplex{T}
    values::VectorVectorRealComplex{T}
    weights::VectorVectorRealComplex{T}
    poles::Vector{Vector{Complex{T}}}
    isbad::Vector{BitVector}
end

function Base.show(io::IO, ::MIME"text/plain", s::ConvergenceStats)
    print(io, f"Convergence stats on {length(s.error)} iterations, best error {s.error[s.bestidx]:.2e}")
end

function Base.show(io::IO, s::ConvergenceStats)
    print(
        IOContext(io,:compact=>true),
        f"Convergence stats on {length(s.error)} iterations, best error {s.error[s.bestidx]:.2e}"
        )
end


#####
##### Approximation on a domain
#####

const unit_interval = Segment(-1.,1.)
const unit_circle = Circle(0., 1.)
const unit_disk = disk(0., 1.)

Domain = Union{ComplexCurve, ComplexPath, ComplexSCRegion}
"""
    Approximation (type)

Approximation of a function on a domain.

# Fields
- `original`: the original function
- `domain`: the domain of the approximation
- `fun`: the barycentric representation of the approximation
- `prenodes`: the prenodes of the approximation
- `stats`: convergence statistics
"""
struct Approximation{T,S} <: Function
    original::Function
    domain::Domain
    fun::AbstractRationalInterpolant{T,S}
    prenodes::Vector{T}
    stats::Union{Missing,ConvergenceStats{T}}
end

(f::Approximation)(z) = f.fun(z)

function Base.show(io::IO, ::MIME"text/plain", f::Approximation)
    print(io, f.fun, " on the domain: ", f.domain)
end

nodes(r::Approximation, args...) = nodes(r.fun, args...)
Base.values(r::Approximation, args...) = values(r.fun, args...)
weights(r::Approximation, args...) = weights(r.fun, args...)
"stats(r) returns the convergence statistics of the rational interpolant `r`."
stats(r::Approximation) = r.stats
degree(r::Approximation) = degree(r.fun)
poles(F::Approximation) = poles(F.fun)
residues(f::Approximation, args...) = residues(f.fun, args...)
roots(f::Approximation) = roots(f.fun)


"""
    rewind(r, degree)

Rewind a rational approximation to a lower degree using stored convergence data.

# Arguments
- `r::Union{Barycentric,Approximation}`: the rational function to rewind
- `degree::Integer`: the degree to rewind to

# Returns
- the rational function of the specified degree (same type as input)

# Examples
```jldoctest
julia> r = aaa(x -> cos(20x))
Barycentric function with 25 nodes and values:
    -1.0=>0.408082,  -0.978022=>0.757786,  -0.912088=>0.820908,  …  1.0=>0.408082

julia> rewind(r, 10)
Barycentric function with 11 nodes and values:
    -1.0=>0.408082,  1.0=>0.408082,  -0.466667=>-0.995822,  …  0.898413=>0.636147
```
"""

function rewind(r::Approximation, degree::Integer)
    m = degree
    M = typeof(r.fun)
    new_r = M(r.stats.nodes[m], r.stats.values[m], r.stats.weights[m])
    return Approximation(r.original, r.domain, new_r, r.prenodes, r.stats)
end


"""
    check(r)

Check the accuracy of a rational approximation `r` on its domain.

# Arguments
- `r::Approximation`: rational approximation

# Returns
- `τ::Vector`: test points
- `err::Vector`: error at test points

See also [`approximate`](@ref), [`aaa`](@ref).
"""
function check(F::Approximation)
    p = F.domain
    if p isa ComplexSCRegion
        p = p.boundary
    end
    t, τ = refine(p, F.prenodes, 30)
    if isreal(F.fun.nodes)
        τ = real(τ)
    end
    err = F.original.(τ) - F.fun.(τ)
    @info f"Max error is {norm(err,Inf):.2e}"
    return τ, err
end
