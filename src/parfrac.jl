# Partial fraction expansion of a rational function.
"""
    PartialFractions (type)

Well-conditioned representation for polynomials on a discrete point set.

# Fields
- `polynomial::ArnoldiPolynomial`: analytic part
- `poles::Vector`: poles
- `residues::Vector`: residues (i.e., coefficients of the partial fractions)
"""
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
    c = isempty(C) ? vectors(B) \ y : [vectors(B) C] \ y
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

# Can't use gradient with mutation in ArnoldiPolynomial.
# TODO: not working
# derivative(r::PartialFractions) = z -> conj(Zygote.forwarddiff(real ∘ r, complex(z))[1])

poles(r::PartialFractions) = r.poles
residues(r::PartialFractions) = (r.poles, r.residues)

# COV_EXCL_START
function Base.show(io::IO, mimetype::MIME"text/plain", r::PartialFractions{T}) where {T}
    ioc = IOContext(io,:compact=>get(io, :compact, true))
    print(ioc, "$T partial fraction expansion of type $(degrees(r))")
    if length(r.poles) > 0
        println(ioc, " with poles=>residues:")
        # print out 3 nodes=>values
        nv, rest = Iterators.peel( zip(residues(r)...) )
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

# Iteratively refine a discretization of a curve/path such that the distance between adjacent
# points is no more than 1/2 the distance to any singularity.
function refine_by_singularity(d::ComplexCurveOrPath, ζ::AbstractVector;
    init=100,
    refinement::Int=2,
    maxpoints::Int=20_000
    )
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

function approximate(::Type{PartialFractions},
    f::Function, d::ComplexCurveOrPath, ζ::AbstractVector;
    degree = max(1, div(length(ζ), 2)),
    init =  max(400, length(d) * 100),
    refinement = 3,
    )

    path = refine_by_singularity(d, ζ; refinement, init)
    _, σ = collect(path, :nodes)
    fσ = f.(σ)
    r = PartialFractions(σ, fσ, ζ, degree)
    return ContinuumApproximation(f, d, r, true, path, nothing)
end

function approximate(::Type{PartialFractions},
    y::AbstractVector, z::AbstractVector, ζ::AbstractVector;
    degree = max(1, div(length(ζ), 2)),
    )
    r = PartialFractions(z, y, ζ, degree)
    return DiscreteApproximation(y, z, r, trues(length(y)), true, nothing)
end

function Base.convert(::Type{PartialFractions}, r::Union{Barycentric, Thiele})
    zp, res = residues(r)
    m, n = degrees(r)
    d = m - n    # degree of the polynomial part
    if d >= 0
        z = first(roots(r), d + 1)
        pf(x) = -sum(res[i] / (x - zp[i]) for i in eachindex(zp))
        p = ArnoldiBasis(z, d) \ pf
    else
        p = ArnoldiPolynomial([0], ArnoldiBasis([1], 0))
    end
    return PartialFractions(p, zp, res)
end
