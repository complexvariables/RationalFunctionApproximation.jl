
##### AbstractRationalInterpolant interface

# parameter S is the value type
abstract type AbstractRationalFunction{S} <: Function end
(r::AbstractRationalFunction)(z) = evaluate(r, z)
Base.eltype(::AbstractRationalFunction{S}) where {S} = S

# COV_EXCL_START
"degrees(r) returns the degrees of the numerator and denominator of the rational `r`."
degrees(r::AbstractRationalFunction) = error("`degrees` not implemented for $(typeof(r))")

"degree(r) returns the degree of the denominator of the rational `r`."
degree(r::AbstractRationalFunction) = error("`degree` not implemented for $(typeof(r))")

"poles(r) returns the poles of the rational interpolant `r`."
poles(::AbstractRationalFunction) = error("`poles` not implemented for $(typeof(r))")

"residues(r) returns two vectors of the poles and residues of the rational function `r`."
residues(::AbstractRationalFunction) = error("`residues` not implemented for $(typeof(r))")

Base.isempty(r::AbstractRationalFunction) = degrees(r) == (0, 0)

"roots(r) returns the roots of the rational function `r`."
roots(::AbstractRationalFunction) = error("`roots` not implemented for $(typeof(r))")

function Base.show(io::IO, mimetype::MIME"text/plain", r::AbstractRationalFunction)
    ioc = IOContext(io, :compact=>get(io, :compact, true))
    println(ioc, "$(typeof(r)) rational function of type $(degrees(r)):")
end

function Base.show(io::IO, r::AbstractRationalFunction)
    print(
        IOContext(io, :compact=>true),
        "$(typeof(r)) rational function of type $(degrees(r))"
    )
end
#COV_EXCL_STOP

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

# parameters are T = float type, S = value type (T or Complex{T})
abstract type AbstractRationalInterpolant{T,S} <: AbstractRationalFunction{S} end

# COV_EXCL_START
# Interface stubs
"nodes(r) returns a vector of the interpolation nodes of the rational interpolant."
nodes(::AbstractRationalInterpolant) = error("`nodes` not implemented for $(typeof(r))")
"values(r) returns a vector of the nodal values of the rational interpolant `r`."
Base.values(::AbstractRationalInterpolant) = error("`values` not implemented for $(typeof(r))")
Base.length(r::AbstractRationalInterpolant) = length(nodes(r))
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
