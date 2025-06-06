# Stable basis for polynomials via the Arnoldi process on a given set of nodes.
struct ArnoldiBasis{T}
    nodes::Vector{T}
    Q::Matrix{T}
    H::Matrix{T}
    function ArnoldiBasis{T}(nodes::AbstractVector{T}, Q::AbstractMatrix{T}, H::AbstractMatrix{T}) where {T}
        @assert size(Q, 1) == length(nodes)
        @assert size(Q, 2) == size(H, 1) == size(H, 2) + 1
        return new{T}(nodes, Q, H)
    end
    function ArnoldiBasis{T}(nodes::AbstractVector, Q::AbstractMatrix, H::AbstractMatrix) where {T}
        nodes = convert.(T, nodes)
        Q = convert.(T, Q)
        H = convert.(T, H)
        return ArnoldiBasis{T}(nodes, Q, H)
    end
    function ArnoldiBasis{T}(B::ArnoldiBasis{S}) where {T,S}
        return ArnoldiBasis{T}(B.nodes, B.Q, B.H)
    end
end
Base.eltype(B::ArnoldiBasis) = eltype(B.nodes)
nodes(b::ArnoldiBasis) = b.nodes

function ArnoldiBasis(z::AbstractVector=ComplexF64[], m::Integer=0)
    n = length(z)
    T = eltype(float(z))
    v = Vector{T}(undef, n)
    Q = similar(v, n, m+1)
    H = similar(v, m+1, m)
    Q[:, 1] .= 1
    for m in 1:m
        v = z .* Q[:, m]
        for k in 1:m
            Qk = view(Q, :, k)
            H[k, m] = dot(Qk, v) / n
            v .-= H[k, m] * Qk
        end
        H[m+1, m] = norm(v) / sqrt(n)
        @. Q[:, m+1] = v / H[m+1, m]
    end
    return ArnoldiBasis{T}(z, Q, H)
end

# Polynomial represented in an Arnoldi basis.
struct ArnoldiPolynomial{T} <: Function
    coeff::Vector{T}
    basis::ArnoldiBasis{T}
    function ArnoldiPolynomial{T}(coeff::AbstractVector{T}, basis::ArnoldiBasis{T}) where {T}
        @assert length(coeff) == size(basis.Q, 2)
        return new{T}(coeff, basis)
    end
    function ArnoldiPolynomial{T}(p::ArnoldiPolynomial{S}) where {T,S}
        return ArnoldiPolynomial{T}(convert.(T, p.coeff), ArnoldiBasis{T}(p.basis))
    end
end

function ArnoldiPolynomial(
    coeff::AbstractVector=[0im],
    basis::ArnoldiBasis{S}=ArnoldiBasis()
    ) where {S}
    T = promote_type(eltype(coeff), S)
    coeff = convert.(T, coeff)
    return ArnoldiPolynomial{T}(coeff, ArnoldiBasis{T}(basis))
end

Base.eltype(p::ArnoldiPolynomial) = eltype(p.coeff)
Base.length(p::ArnoldiPolynomial) = length(p.coeff)
degree(p::ArnoldiPolynomial) = length(p.coeff) - 1
nodes(p::ArnoldiPolynomial) = nodes(p.basis)

# COV_EXCL_START
function Base.show(io::IO, mimetype::MIME"text/plain", p::ArnoldiPolynomial)
    ioc = IOContext(io,:compact=>get(io, :compact, true))
    len = length(p)
    if len==0
        print(ioc, "Empty $(typeof(p)) Arnoldi polynomial")
    else
        print(ioc, "Arnoldi polynomial of degree $(len-1)")
    end
end
# COV_EXCL_STOP

(p::ArnoldiPolynomial)(z) = evaluate(p, z)
function evaluate(p::ArnoldiPolynomial, z::Number)
    g = p.coeff[1]
    H = p.basis.H
    n = size(H, 2)
    Q = fill(one(g), n+1)
    for k in 1:n-1
        v = z .* Q[k]
        for j in 1:k
            v -= H[j, k] * Q[j]
        end
        Q[k+1] = v / H[k+1, k]
        g += p.coeff[k+1] * Q[k+1]
    end
    return g
end

# Partial fraction expansion of a rational function.
struct PartialFractions{S} <: AbstractRationalFunction{S}
    polynomial::ArnoldiPolynomial{S}
    poles::Vector{S}
    residues::Vector{S}
    function PartialFractions{S}(p::ArnoldiPolynomial{S}, poles::AbstractVector{S}, residues::AbstractVector{S}) where {S}
        @assert length(poles) == length(residues)
        return new{S}(p, poles, residues)
    end
end

function PartialFractions(
    p::ArnoldiPolynomial = ArnoldiPolynomial(),
    poles::AbstractVector = ComplexF64[],
    residues::AbstractVector = ComplexF64[]
    )
    @assert length(poles) == length(residues)
    if isempty(poles)
        poles = residues = Vector{eltype(p)}()
    end
    T = promote_type(eltype(p), eltype(poles), eltype(residues))
    poles = convert.(T, poles)
    residues = convert.(T, residues)
    return PartialFractions{T}(ArnoldiPolynomial{T}(p), poles, residues)
end

function PartialFractions(z::AbstractVector, y::AbstractVector, ζ::AbstractVector, degree::Int)
    B = ArnoldiBasis(z, degree)
    C = [1 / (z - zp) for z in z, zp in ζ]
    c = isempty(C) ? B.Q \ y : [B.Q C] \ y
    p = ArnoldiPolynomial(c[1:degree+1], B)
    return PartialFractions(p, ζ, c[degree+2:end])
end

function degrees(r::PartialFractions)
    d = degree(r.polynomial)
    n = length(r.poles)
    return d + n, n
end
degree(r::PartialFractions) = length(r.poles)

(f::PartialFractions)(z) = evaluate(f, z)
function evaluate(r::PartialFractions, z::Number)
    u = r.polynomial(z)
    for (r, p) in zip(r.residues, r.poles)
        u += r / (z - p)
    end
    return u
end

poles(r::PartialFractions) = r.poles
residues(r::PartialFractions) = (r.poles, r.residues)

# COV_EXCL_START
function Base.show(io::IO, mimetype::MIME"text/plain", r::PartialFractions{T}) where {T}
    ioc = IOContext(io,:compact=>get(io, :compact, true))
    print(ioc, "$T partial fraction expansion of type $(degrees(r))")
    if length(r.poles) > 0
        println(ioc, " with poles=>residues:")
        # print out 3 nodes=>values
        nv, rest = Iterators.peel( zip(poles(r), residues(r)) )
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
end
# COV_EXCL_STOP

function refine_by_singularity(d::ComplexCurveOrPath, ζ::AbstractVector;
    init=100,
    refinement::Int=2,
    maxpoints::Int=20_000
    )
    # Iteratively refine a discretization of a curve/path such that the distance between adjacent points is no more than 1/2 the distance to any singularity.
    path = DiscretizedPath(d, range(0, length(d), init+1); refinement, maxpoints)
    isempty(ζ) && return path
    z = [1.]
    while length(z) < 19_000
        Δ = spacing(path)
        m, n = size(Δ)
        z = path.points[1:m, 2:n]
        s = [minimum(abs(z - w) for w in ζ) for z in z]
        ρ, idx = findmax(Δ[:, 2:n] ./ s)
        if ρ <= 0.5
            return path
        end
        add_node!(path, (idx[1], idx[2]+1))
    end
    @warn "Refinement was not successful"
    return path
end
