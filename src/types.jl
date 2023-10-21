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
##### convergence statistics
#####
"""
    ConvergenceStats{T}(bestidx, error, nbad, nodes, values, weights, poles)

Convergence statistics for a sequence of rational approximations.

# Fields
- `bestidx`: the index of the best approximation
- `error`: the error of each approximation
- `nbad`: the number of bad nodes in each approximation
- `nodes`: the nodes of each approximation
- `values`: the values of each approximation
- `weights`: the weights of each approximation
- `poles`: the poles of each approximation

See also: [`approximate`](@ref), [`Barycentric`](@ref)
"""
struct ConvergenceStats{T}
    bestidx::Int
    error::Vector{<:AbstractFloat}
    nbad::Vector{Int}
    nodes::VectorVectorRealComplex{T}
    values::VectorVectorRealComplex{T}
    weights::VectorVectorRealComplex{T}
    poles::Vector{Vector{Complex{T}}}
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
##### barycentric representation
#####
"""
    Barycentric (type)

Barycentric representation of a rational function.

# Fields
- `node`: the nodes of the rational function
- `value`: the values of the rational function
- `weight`: the weights of the rational function
- `wf`: the weighted values of the rational function
- `stats`: convergence statistics
"""
struct Barycentric{T,S} <: Function
    nodes::Vector{S}
    values::Vector{S}
    weights::Vector{S}
    w_times_f::Vector{S}
    stats::Union{Missing,ConvergenceStats{T}}
    function Barycentric{T}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S},
        wf::AbstractVector{S} = value.*weight;
        stats::Union{Missing,ConvergenceStats{T}} = missing
        )  where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight) == length(wf)
        new{T,S}(node, value, weight, wf, stats)
    end
end

"""
    Barycentric(node, value, weight, wf=value.*weight; stats=missing)

Construct a `Barycentric` rational function.

# Arguments
- `node::AbstractVector`: interpolation nodes
- `value::AbstractVector`: values at the interpolation nodes
- `weight::AbstractVector`: barycentric weights
- `wf::AbstractVector`: weights times values (optional)
- `stats::ConvergenceStatistics``: convergence statistics (optional)

# Examples
```jldoctest
julia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])
Barycentric function with 3 nodes and values:
    1.0=>1.0,  2.0=>2.0,  3.0=>3.0

julia> r(1.5)
1.5
```
"""
function Barycentric(node, value, weight, wf=value.*weight; stats=missing)
    Barycentric( promote(float(node), float(value), float(weight))..., float(wf); stats )
end

function Barycentric(
    node::Vector{S}, value::Vector{S}, weight::Vector{S}, wf=value.*weight;
    stats=missing
    ) where {T<:AbstractFloat, S<:RealComplex{T}}
    return Barycentric{T}(node, value, weight, wf; stats)
end

# convenience accessors and overloads
"nodes(r) returns the nodes of the rational interpolant `r` as a vector."
nodes(r::Barycentric) = r.nodes
nodes(r::Barycentric, m::Integer) = r.stats.nodes[m]
"values(r) returns the nodal values of the rational interpolant `r` as a vector."
Base.values(r::Barycentric) = r.values
Base.values(r::Barycentric, m::Integer) = r.stats.values[m]
"weights(r) returns the weights of the rational interpolant `r` as a vector."
weights(r::Barycentric) = r.weights
weights(r::Barycentric, m::Integer) = r.stats.weights[m]
"stats(r) returns the convergence statistics of the rational interpolant `r`."
stats(r::Barycentric) = r.stats
Base.eltype(r::Barycentric{T}) where T = eltype(nodes(r))
Base.length(r::Barycentric) = length(nodes(r))
"degree(r) returns the degree of the numerator and denominator of the rational `r`."
degree(r::Barycentric) = length(r.nodes) - 1

"""
    rewind(r, degree)

Rewind a `Barycentric` rational function to a lower degree using stored convergence data.

# Arguments
- `r::Union{Barycentric,Approximation}`: the rational function to rewind
- `degree::Integer`: the degree to rewind to

# Returns
- the rational function of the specified degree (same type as input)

# Examples
```jldoctest
julia> r = aaa(x -> cos(20x), stats=true)
Barycentric function with 25 nodes and values:
    -1.0=>0.408082,  -0.978022=>0.757786,  -0.912088=>0.820908,  …  1.0=>0.408082

julia> rewind(r, 10)
Barycentric function with 11 nodes and values:
    -1.0=>0.408082,  1.0=>0.408082,  -0.466667=>-0.995822,  …  0.898413=>0.636147
```
"""
function rewind(r::Barycentric, degree::Integer)
    m = degree
    return Barycentric(r.stats.nodes[m], r.stats.values[m], r.stats.weights[m])
end

function Base.show(io::IO, mimetype::MIME"text/plain", r::Barycentric{T,S}) where {T,S}
    ioc = IOContext(io,:compact=>get(io, :compact, true))
    println(ioc, "Barycentric function with $(length(r)) nodes and values:")
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

function Base.show(io::IO, r::Barycentric{T,S}) where {T,S}
    m = length(nodes(r))
    print(
        IOContext(io,:compact=>true),
        "Barycentric rational function of type ($(m-1),$(m-1))"
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
"""
struct Approximation{T,S} <: Function
    original::Function
    domain::Domain
    fun::Barycentric{T,S}
    prenodes::Vector{T}
end

function Base.show(io::IO, ::MIME"text/plain", f::Approximation)
    print(io, f.fun, " on the domain: ", f.domain)
end

nodes(r::Approximation, args...) = nodes(r.fun, args...)
Base.values(r::Approximation, args...) = values(r.fun, args...)
weights(r::Approximation, args...) = weights(r.fun, args...)
stats(r::Approximation) = stats(r.fun)
degree(r::Approximation) = degree(r.fun)
rewind(f::Approximation, degree::Integer) = Approximation(f.original, f.domain, rewind(f.fun, degree), f.prenodes)
