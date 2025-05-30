
##### AbstractRationalInterpolant interface

# subtypes are T = float type, S = value type (T or Complex{T})
abstract type AbstractRationalInterpolant{T,S} <: Function end

# COV_EXCL_START
# Interface stubs
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

"""
    residues(r)

Returns two vectors of the poles and residues of the rational function `r`.
"""
residues(::AbstractRationalInterpolant) = error("`residues` not implemented for $(typeof(r))")

"""
    Res(r, z)

Returns the residue of the rational function `r` at the point `z`.
"""
function Res(f::Function, z::Number; radius=100eps(abs(z)), n=200)
    # Attempt to use contour integration to find a residue
    T = complex(typeof(float(z)))
    trap = 0
    for t in 0:n-1
        c = cispi(T(2t / n))
        trap += c * f(z + radius * c)
    end
    return radius * trap / n
end

"roots(r) returns the roots of the rational interpolant `r`."
roots(::AbstractRationalInterpolant) = error("`roots` not implemented for $(typeof(r))")
# COV_EXCL_STOP

# COV_EXCL_START
function Base.show(io::IO, mimetype::MIME"text/plain", r::AbstractRationalInterpolant)
    ioc = IOContext(io,:compact=>get(io, :compact, true))
    len = length(r)
    if len==0
        println(ioc, "Empty $(typeof(r)) rational interpolant")
        return
    end
    println(ioc, "$(typeof(r)) rational interpolant of type $(degrees(r)):")
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
        "$(typeof(r)) rational interpolant of type $(degrees(r))"
    )
end
# COV_EXCL_STOP

##### AbstractRFIVector for storing convergence history

abstract type AbstractRFIVector{T,S} <: AbstractVector{AbstractRationalInterpolant{T,S}} end
"""
Sequence of rational interpolants produced by an iteration.

# Fields
- `nodes`: vector of interpolation nodes
- `values`: vector of interpolation values
- `weights`: matrix of all weights (upper triangle)
- `len`: the number of nodes for each approximation
- `best`: the index of the best approximation

See also: [`approximate`](@ref)
"""
struct RFIVector{R} <: AbstractVector{R}
    nodes::AbstractVector
    values::AbstractVector
    weights::AbstractMatrix
    len::Vector{Int}
    best::Int
end

RFIVector{R}() where R<:AbstractRationalInterpolant = RFIVector{R}([], [], zeros(0, 0), Int[], 0)

# must construct RFIVector with type explicitly, since otherwise it doesn't know what kind
# of rational interpolant to use/make
Base.size(v::RFIVector) = (length(v.len),)
Base.IndexStyle(::Type{<:RFIVector}) = IndexLinear()
function Base.getindex(v::RFIVector{R}, idx::Int) where {R <: AbstractRationalInterpolant}
    L = v.len[idx]
    return R(promote(v.nodes[1:L], v.values[1:L], v.weights[1:L, L])...)
end
