


import Base: values

export QuadraticThiele, QTCF, SymmetricQTCF, Symmetric, symmetric

struct Symmetric{D,R}
    domain::D
    reflection::R
end

function Symmetric(d; reflection=nothing)
    return Symmetric{typeof(d),typeof(reflection)}(d, reflection)
end

symmetric(d; kw...) = Symmetric(d; kw...)




struct QuadraticThiele{T,S} <: RFA.AbstractRationalInterpolant{T,S}
    P::Vector{S}
    Q::Vector{S}
    pairs::Vector{Tuple{S,S}}
    bpolys::Vector{Vector{S}}
    apolys::Vector{Vector{S}}
    nodes::Vector{S}
    vals::Vector{S}
end

const QTCF = QuadraticThiele
const SymmetricQTCF = :symmetric

function _qtcf_strategy_symbol(strategy)
    strategy isa Symbol && return strategy
    strategy isa AbstractString && return Symbol(strategy)
    return strategy
end



_rtyp(::Type{T}) where {T<:Real} = T
_rtyp(::Type{Complex{T}}) where {T<:Real} = T





function QuadraticThiele(P, Q, pairs, bpolys, apolys, z, y)
    S = promote_type(eltype(P), eltype(Q), eltype(z), eltype(y))
    if any(_qtcf_infinite_block, bpolys)
        P, Q = _qtcf_continuants([S.(b) for b in bpolys], [S.(a) for a in apolys])
    end
    return QuadraticThiele{_rtyp(S),S}(
        S.(P), S.(Q), Tuple{S,S}.(pairs),
        Vector{S}.(bpolys), Vector{S}.(apolys), S.(z), S.(y))
end


QuadraticThiele() = QuadraticThiele(
    [0.0], [1.0], Tuple{Float64,Float64}[], Vector{Float64}[],
    Vector{Float64}[], Float64[], Float64[])



function Base.copy(g::QuadraticThiele)
    return QuadraticThiele(copy(g.P), copy(g.Q), copy(g.pairs),
        copy.(g.bpolys), copy.(g.apolys), copy(g.nodes), copy(g.vals))
end



nodes(g::QuadraticThiele) = g.nodes
values(g::QuadraticThiele) = g.vals
weights(g::QuadraticThiele) = g.bpolys
Base.isreal(g::QuadraticThiele) = all(isreal, g.P) && all(isreal, g.Q)

_is_selected_node(g::QuadraticThiele, z) =
    any(w -> _same(convert(eltype(g), z), w), nodes(g))


function _ptrim(p; atol=0.0)
    q = copy(p)
    while length(q) > 1 && abs(last(q)) <= atol
        pop!(q)
    end
    return q
end



_pdegree(p) = length(_ptrim(p)) - 1



degrees(g::QuadraticThiele) = (_pdegree(g.P), _pdegree(g.Q))


degree(g::QuadraticThiele) = _pdegree(g.Q)

function _peval(p, z)
    y = zero(promote_type(eltype(p), typeof(z)))
    @inbounds for i in length(p):-1:1
        y = muladd(y, z, p[i])
    end
    return y
end

_qtcf_infinite_block(b) = length(b) == 1 && isinf(only(b))

function _qtcf_ratio(numer, denom)
    n, d = promote(numer, denom)
    (isnan(n) || isnan(d)) && return oftype(d, NaN)
    iszero(d) && return oftype(d, iszero(n) ? NaN : Inf)
    isinf(d) && return oftype(d, isinf(n) ? NaN : 0)
    return n / d
end

function _qtcf_update_continuants(b, a, P, Q, Pold, Qold)
    if _qtcf_infinite_block(b)
        return copy(P), copy(Q), zero.(Pold), zero.(Qold)
    end
    Pnew = _padd(_pmul(b, P), _pmul(a, Pold))
    Qnew = _padd(_pmul(b, Q), _pmul(a, Qold))
    return Pnew, Qnew, P, Q
end

function _qtcf_continuants(bpolys, apolys)
    S = eltype(first(bpolys))
    Pold, Qold = S[1], S[0]
    if _qtcf_infinite_block(first(bpolys))
        P, Q = S[1], S[0]
        Pold, Qold = S[0], S[0]
    else
        P, Q = copy(first(bpolys)), S[1]
    end
    for j in 2:length(bpolys)
        P, Q, Pold, Qold = _qtcf_update_continuants(
            bpolys[j], apolys[j-1], P, Q, Pold, Qold)
    end
    return P, Q
end

function _qtcf_eval_b(b, pair, z)
    length(b) == 1 && return only(b)
    if pair !== nothing && length(b) == 2
        midpoint = _qtcf_ldexp(pair[1], -1) + _qtcf_ldexp(pair[2], -1)
        midpoint_value = b[1] + b[2] * midpoint
        return midpoint_value + b[2] * (z - midpoint)
    end
    return _peval(b, z)
end

function _qtcf_eval_a(a, pair, z)
    if pair !== nothing && length(a) == 2
        return a[2] * (z - pair[1])
    elseif pair !== nothing && length(a) == 3
        return a[3] * (z - pair[1]) * (z - pair[2])
    end
    return _peval(a, z)
end

function _qtcf_eval_a_deriv(a, pair, z)
    if pair !== nothing && length(a) == 2
        return a[2]
    elseif pair !== nothing && length(a) == 3
        return a[3] * ((z - pair[1]) + (z - pair[2]))
    end
    return _peval(_pderiv(a), z)
end




function _padd(p, q)
    n = max(length(p), length(q))
    r = zeros(promote_type(eltype(p), eltype(q)), n)
    r[1:length(p)] .+= p
    r[1:length(q)] .+= q
    return _ptrim(r)
end




function _pmul(p, q)
    r = zeros(promote_type(eltype(p), eltype(q)), length(p) + length(q) - 1)
    @inbounds for i in eachindex(p), j in eachindex(q)
        r[i+j-1] += p[i] * q[j]
    end
    return _ptrim(r)
end




function _pderiv(p)
    length(p) <= 1 && return [zero(eltype(p))]
    return [k * p[k+1] for k in 1:length(p)-1]
end

_qtcf_component_abs(x::Real) = abs(x)
_qtcf_component_abs(x::Complex) = max(abs(real(x)), abs(imag(x)))
_qtcf_ldexp(x::Real, exponent::Integer) = ldexp(x, exponent)
_qtcf_ldexp(x::Complex, exponent::Integer) =
    complex(ldexp(real(x), exponent), ldexp(imag(x), exponent))

function _qtcf_normalize_pair(numer, denom)
    scale = max(_qtcf_component_abs(numer), _qtcf_component_abs(denom))
    if isfinite(scale) && scale > 0
        shift = -exponent(scale)
        return _qtcf_ldexp(numer, shift), _qtcf_ldexp(denom, shift)
    end
    return numer, denom
end

function _evaluate_numden(g::QuadraticThiele, z::Number)
    @assert isfinite(z)
    m = length(g.bpolys)
    if m == 0
        return _peval(g.P, z), _peval(g.Q, z)
    end
    S = promote_type(eltype(g.P), eltype(g.Q), typeof(z))
    pairm = m <= length(g.pairs) ? g.pairs[m] : nothing
    infinite = _qtcf_infinite_block(g.bpolys[m])
    numer = infinite ? one(S) : S(_qtcf_eval_b(g.bpolys[m], pairm, z))
    denom = infinite ? zero(S) : one(S)
    numer, denom = _qtcf_normalize_pair(numer, denom)
    @inbounds for k in (m-1):-1:1
        pair = k <= length(g.pairs) ? g.pairs[k] : nothing
        b = _qtcf_eval_b(g.bpolys[k], pair, z)
        a = _qtcf_eval_a(g.apolys[k], pair, z)
        scale = max(_qtcf_component_abs(a), _qtcf_component_abs(b))
        if isfinite(scale) && scale > floatmax(typeof(scale)) * 0.125
            numer = _qtcf_ldexp(numer, -3)
            denom = _qtcf_ldexp(denom, -3)
        end
        numer, denom = _qtcf_infinite_block(g.bpolys[k]) ?
            (numer, zero(S)) : (b * numer + a * denom, numer)
        numer, denom = _qtcf_normalize_pair(numer, denom)
    end
    return numer, denom
end



function _evaluate_classic(g::QuadraticThiele, z::Number)
    isempty(g.bpolys) && return _qtcf_ratio(_peval(g.P, z), _peval(g.Q, z))
    m = length(g.bpolys)
    pair = m <= length(g.pairs) ? g.pairs[m] : nothing
    u = _qtcf_eval_b(g.bpolys[m], pair, z)
    for k in m-1:-1:1
        pair = k <= length(g.pairs) ? g.pairs[k] : nothing
        b = _qtcf_eval_b(g.bpolys[k], pair, z)
        a = _qtcf_eval_a(g.apolys[k], pair, z)
        if _qtcf_infinite_block(g.bpolys[k])
            u = oftype(u, iszero(u) || isnan(u) ? NaN : Inf)
        else
            u = b + _qtcf_ratio(a, u)
        end
    end
    return u
end

function _evaluate_onediv(g::QuadraticThiele, z::Number)
    numer, denom = _evaluate_numden(g, z)
    return if iszero(denom)
        @debug "Evaluation produced a division by zero at " z
        _evaluate_classic(g, z)
    else
        numer / denom
    end
end

_evaluate(g::QuadraticThiele, z, ::OneDiv) = _evaluate_onediv(g, z)
_evaluate(g::QuadraticThiele, z, ::Classic) = _evaluate_classic(g, z)

function evaluate(g::QuadraticThiele, z::Number,
    method::ThieleMethod=default_eval_method())
    if isinf(z)
        m, n = degrees(g)
        return m > n ? oftype(last(g.P), Inf) :
            m < n ? zero(last(g.P)) : _qtcf_ratio(last(g.P), last(g.Q))
    elseif isnan(z)
        return oftype(z, NaN)
    end
    return _evaluate(g, z, method)
end

function evaluate(g::QuadraticThiele, z::AbstractArray{<:Number},
    method::ThieleMethod=default_eval_method())
    S = promote_type(eltype(g.P), eltype(g.Q), eltype(z))
    return evaluate!(similar(z, S), g, z, method)
end

function evaluate!(t::AbstractArray, g::QuadraticThiele,
    z::AbstractArray{<:Number}, method::ThieleMethod=default_eval_method())
    return evaluate!(t, g, z, similar(t), similar(t), method)
end

function evaluate!(t::AbstractArray, g::QuadraticThiele,
    z::AbstractArray{<:Number}, a::AbstractArray, b::AbstractArray,
    method::ThieleMethod=default_eval_method())
    axes(t) == axes(z) == axes(a) == axes(b) ||
        throw(DimensionMismatch("evaluation input, output, and scratch arrays must have matching axes"))
    if Base.mightalias(a, b) || Base.mightalias(a, t) || Base.mightalias(b, t) ||
        Base.mightalias(a, z) || Base.mightalias(b, z)
        throw(ArgumentError("evaluation scratch arrays must not overlap each other, the input, or the output"))
    end
    source = Base.mightalias(t, z) ? copy(z) : z
    if !(method isa OneDiv)
        for i in eachindex(t, source)
            t[i] = evaluate(g, source[i], method)
        end
    else
        for i in eachindex(a, b, source)
            if isfinite(source[i])
                a[i], b[i] = _evaluate_numden(g, source[i])
            else
                a[i] = evaluate(g, source[i], method)
                b[i] = one(eltype(b))
            end
        end
        for i in eachindex(t, a, b, source)
            t[i] = !isfinite(source[i]) ? a[i] :
                iszero(b[i]) ? _evaluate_classic(g, source[i]) : a[i] / b[i]
        end
    end
    return t
end

function _derivative!(z, g::QuadraticThiele, order::Integer,
    A::AbstractVector, B::AbstractVector, vals::AbstractVector)
    order >= 0 || throw(ArgumentError("derivative order must be nonnegative"))
    min(length(A), length(B), length(vals)) >= order + 1 ||
        throw(DimensionMismatch("derivative workspaces must contain at least order + 1 entries"))
    firstindex(A) == firstindex(B) == firstindex(vals) == 1 ||
        throw(ArgumentError("derivative workspaces must use one-based indexing"))
    if Base.mightalias(A, B) || Base.mightalias(A, vals) || Base.mightalias(B, vals)
        throw(ArgumentError("derivative workspaces must not overlap"))
    end
    numer, denom = _qtcf_numden_derivatives(g, z, order)
    copyto!(A, 1, numer, 1, order + 1)
    copyto!(B, 1, denom, 1, order + 1)
    recovered = _qtcf_quotient_derivatives(numer, denom)
    copyto!(vals, 1, recovered, 1, order + 1)
    return vals
end

function _qtcf_b_derivatives(b, pair, z, order)
    S = promote_type(ComplexF64, eltype(b), typeof(z))
    jet = zeros(S, order + 1)
    jet[1] = S(_qtcf_eval_b(b, pair, z))
    order > 0 && length(b) > 1 && (jet[2] = S(b[2]))
    return jet
end

function _qtcf_a_derivatives(a, pair, z, order)
    S = promote_type(ComplexF64, eltype(a), typeof(z))
    jet = zeros(S, order + 1)
    if pair !== nothing && length(a) == 2
        c, z1 = S(a[2]), S(pair[1])
        jet[1] = c * (z - z1)
        order > 0 && (jet[2] = c)
    elseif pair !== nothing && length(a) == 3
        c, z1, z2 = S(a[3]), S(pair[1]), S(pair[2])
        jet[1] = c * (z - z1) * (z - z2)
        order > 0 && (jet[2] = c * ((z - z1) + (z - z2)))
        order > 1 && (jet[3] = 2c)
    else
        p = S.(a)
        for k in 0:order
            jet[k + 1] = _peval(p, z)
            p = _pderiv(p)
        end
    end
    return jet
end

function _qtcf_shift_jets!(numer, denom, shift)
    for jet in (numer, denom), k in eachindex(jet)
        jet[k] = _qtcf_ldexp(jet[k], shift)
    end
    return numer, denom
end

function _qtcf_normalize_jets!(numer, denom)
    scale = max(maximum(_qtcf_component_abs, numer; init=0),
        maximum(_qtcf_component_abs, denom; init=0))
    if isfinite(scale) && scale > 0
        _qtcf_shift_jets!(numer, denom, -exponent(scale))
    end
    return numer, denom
end

function _qtcf_numden_derivatives(g::QuadraticThiele, z::Number, order::Integer)
    S = promote_type(ComplexF64, eltype(g.P), eltype(g.Q), typeof(z))
    m = length(g.bpolys)
    if m == 0 || any(_qtcf_infinite_block, g.bpolys)
        numer = zeros(S, order + 1)
        denom = zeros(S, order + 1)
        p, q = S.(g.P), S.(g.Q)
        for k in 0:order
            numer[k + 1] = _peval(p, z)
            denom[k + 1] = _peval(q, z)
            p, q = _pderiv(p), _pderiv(q)
        end
        return numer, denom
    end

    pairm = m <= length(g.pairs) ? g.pairs[m] : nothing
    numer = _qtcf_b_derivatives(g.bpolys[m], pairm, z, order)
    denom = zeros(S, order + 1)
    denom[1] = one(S)
    _qtcf_normalize_jets!(numer, denom)

    for k in (m - 1):-1:1
        pair = k <= length(g.pairs) ? g.pairs[k] : nothing
        bjet = _qtcf_b_derivatives(g.bpolys[k], pair, z, order)
        ajet = _qtcf_a_derivatives(g.apolys[k], pair, z, order)
        scale = max(maximum(_qtcf_component_abs, bjet; init=0),
            maximum(_qtcf_component_abs, ajet; init=0))
        if isfinite(scale) && scale > 0
            margin = 5 + 2ndigits(max(order, 1); base=2)
            shift = min(0, exponent(floatmax(typeof(scale))) - margin - exponent(scale))
            shift < 0 && _qtcf_shift_jets!(numer, denom, shift)
        end
        next_numer, next_denom = numer, denom
        numer = zeros(S, order + 1)
        for q in 0:order, j in 0:min(q, 2)
            c = binomial(q, j)
            numer[q + 1] += c * (
                bjet[j + 1] * next_numer[q - j + 1] +
                ajet[j + 1] * next_denom[q - j + 1])
        end
        denom = next_numer

        _qtcf_normalize_jets!(numer, denom)
    end
    return numer, denom
end

function _qtcf_quotient_derivatives(numer, denom)
    order = length(numer) - 1
    vals = similar(numer)
    inv_denom = inv(denom[1])
    vals[1] = numer[1] * inv_denom
    for m in 1:order
        correction = zero(eltype(vals))
        for j in 0:m-1
            correction += binomial(m, j) * vals[j + 1] * denom[m - j + 1]
        end
        vals[m + 1] = (numer[m + 1] - correction) * inv_denom
    end
    return vals
end

function derivative(g::QuadraticThiele, order::Integer=1)
    f = derivative(g, [order])
    return z -> only(f(z))
end

function derivative(g::QuadraticThiele, orders::AbstractVector{<:Integer})
    isempty(orders) && throw(ArgumentError("at least one derivative order is required"))
    any(<(0), orders) && throw(ArgumentError("derivative orders must be nonnegative"))
    order = maximum(orders)
    index = 1 .+ orders
    return function(z)
        numer, denom = _qtcf_numden_derivatives(g, z, order)
        vals = _qtcf_quotient_derivatives(numer, denom)
        isreal(z) && isreal(g) && (vals = real(vals))
        return vals[index]
    end
end



function _evaluate_numden_derivs(g::QuadraticThiele, z::Number)
    numer, denom = _qtcf_numden_derivatives(g, z, 1)
    S = promote_type(eltype(g.P), eltype(g.Q), typeof(z))
    return S(numer[1]), S(denom[1]), S(numer[2]), S(denom[2])
end


function _add_affine_entry!(C, D, i, j, p)
    T = eltype(C)
    pp = T.(p)
    C[i, j] += pp[1]
    length(pp) >= 2 && (D[i, j] += pp[2])
    return nothing
end

function _block_offdiag_factors(g::QuadraticThiele, k::Int, ::Type{S}) where {S}
    a = S.(g.apolys[k])
    if length(a) <= 2
        return a, S[-one(S)]
    end
    z1, z2 = S(g.pairs[k][1]), S(g.pairs[k][2])
    return S[-z1, one(S)], S[z2, -one(S)]
end

function _qtcf_denominator_linear_pencil(g::QuadraticThiele)
    S = promote_type(eltype(g.Q), eltype(g.nodes))
    n = max(length(g.bpolys) - 1, 0)
    C = zeros(S, n, n)
    D = zeros(S, n, n)
    @inbounds for j in 1:n
        b = S.(g.bpolys[j+1])
        _add_affine_entry!(C, D, j, j, b)
        if j >= 2
            ell, u = _block_offdiag_factors(g, j, S)
            _add_affine_entry!(C, D, j-1, j, ell)
            _add_affine_entry!(C, D, j, j-1, u)
        end
    end
    return C, D
end

function poles(g::QuadraticThiele)
    (isempty(g.bpolys) || any(_qtcf_infinite_block, g.bpolys)) &&
        return _companion_roots(g.Q)
    C, D = _qtcf_denominator_linear_pencil(g)
    n = size(C, 1)
    n <= 0 && return eltype(g.nodes)[]
    z = try
        filter(isfinite, eigvals(-C, D))
    catch
        _, _, _, _, α, β = schur(complex(-C), complex(D))
        filter(isfinite, α ./ β)
    end
    return z
end

function residues(g::QuadraticThiele)
    ζ = poles(g)
    res = similar(complex(ζ))
    for i in eachindex(ζ)
        p, q, pʹ, qʹ = _evaluate_numden_derivs(g, ζ[i])
        candidate = iszero(qʹ) ? oftype(p, NaN) : p / qʹ
        if isfinite(candidate)
            res[i] = candidate
        else
            @debug "Fallback for residue at pole " ζ[i]
            radius = max(one(abs(ζ[i])), abs(ζ[i])) * 0.25
            res[i] = Res(g, ζ[i]; avoid=ζ, radius)
        end
    end
    return ζ, res
end

function _pdivrem(P, Q)
    P, Q = _ptrim(P), _ptrim(Q)
    all(iszero, Q) && throw(DivideError())
    S = promote_type(eltype(P), eltype(Q))
    p, q = S.(P), S.(Q)
    dp, dq = _pdegree(p), _pdegree(q)
    dp < dq && return S[0], p

    quotient = zeros(S, dp - dq + 1)
    remainder = copy(p)
    @inbounds for shift in (dp-dq):-1:0
        coefficient = remainder[dq + shift + 1] / q[dq + 1]
        quotient[shift + 1] = coefficient
        for k in 0:dq
            remainder[k + shift + 1] -= coefficient * q[k + 1]
        end
    end
    return _ptrim(quotient), _ptrim(remainder[1:min(dq, length(remainder))])
end

function _companion_roots(p)
    p = _ptrim(p)
    n = _pdegree(p)
    S = promote_type(ComplexF64, eltype(p))
    n == 0 && return S[]
    C = zeros(S, n, n)
    @inbounds for j in 1:n-1
        C[j+1, j] = one(S)
    end
    C[:, n] .= -S.(p[1:n]) ./ S(p[n+1])
    return eigvals(C)
end

function _qtcf_polynomial_part(quotient, g::QuadraticThiele)
    S = promote_type(ComplexF64, eltype(quotient), eltype(g.nodes))
    if all(iszero, quotient)
        return ArnoldiPolynomial(S[0], ArnoldiBasis(S[0], 0))
    end
    d = _pdegree(quotient)
    center = isempty(g.nodes) ? zero(S) : sum(S.(g.nodes)) / length(g.nodes)
    radius = max(1.0, maximum(z -> Float64(abs(S(z) - center)), g.nodes; init=0.0))
    sample = d == 0 ? S[center] :
        S[center + radius * cis(2pi * k / (d + 1)) for k in 0:d]
    basis = ArnoldiBasis(sample, d)
    values = [_peval(quotient, z) for z in sample]
    coefficients = vectors(basis) \ values
    return ArnoldiPolynomial(coefficients, basis)
end

function Base.convert(::Type{PartialFractions}, g::QuadraticThiele)
    quotient, _ = _pdivrem(g.P, g.Q)
    polynomial = _qtcf_polynomial_part(quotient, g)
    ζ, residue = residues(g)
    any(ζ[i] == ζ[j] for i in eachindex(ζ) for j in firstindex(ζ):i-1) &&
        throw(ArgumentError("PartialFractions conversion requires simple poles"))
    return PartialFractions(polynomial, ζ, residue)
end

function Base.convert(::Type{F}, g::QuadraticThiele{T,S}) where {F<:AbstractFloat,T,S}
    V = S <: Complex ? Complex{F} : F
    return QuadraticThiele{F,V}(
        V.(g.P), V.(g.Q), Tuple{V,V}.(g.pairs),
        [V.(b) for b in g.bpolys], [V.(a) for a in g.apolys],
        V.(g.nodes), V.(g.vals))
end

function Base.:+(g::QuadraticThiele, s::Number)
    S = promote_type(eltype(g), typeof(s))
    bpolys = [S.(b) for b in g.bpolys]
    apolys = [S.(a) for a in g.apolys]
    if isempty(bpolys)
        P, Q = _padd(S.(g.P), s .* S.(g.Q)), S.(g.Q)
    else
        bpolys[1][1] += s
        P, Q = _qtcf_continuants(bpolys, apolys)
    end
    return QuadraticThiele(P, Q, copy(g.pairs), bpolys, apolys,
        copy(g.nodes), g.vals .+ s)
end

function Base.:-(g::QuadraticThiele)
    return g * -one(eltype(g))
end

function Base.:*(g::QuadraticThiele, s::Number)
    S = promote_type(eltype(g), typeof(s))
    if iszero(s)
        z = S.(g.nodes[1:min(1, length(g.nodes))])
        pairs = isempty(z) ? Tuple{S,S}[] : [(only(z), only(z))]
        bpolys = isempty(z) ? Vector{S}[] : [S[0]]
        apolys = isempty(z) ? Vector{S}[] : [S[-only(z), 1]]
        return QuadraticThiele(S[0], S[1], pairs, bpolys, apolys, z, zero.(z))
    end
    bpolys = [isodd(j) ? s .* S.(b) : S.(b) ./ s for (j, b) in enumerate(g.bpolys)]
    apolys = [S.(a) for a in g.apolys]
    P, Q = isempty(bpolys) ? (s .* S.(g.P), S.(g.Q)) :
        _qtcf_continuants(bpolys, apolys)
    return QuadraticThiele(P, Q, copy(g.pairs), bpolys, apolys,
        copy(g.nodes), s .* g.vals)
end

function _qtcf_numerator_linear_pencil(g::QuadraticThiele)
    S = promote_type(eltype(g.P), eltype(g.nodes))
    n = length(g.bpolys)
    C = zeros(S, n, n)
    D = zeros(S, n, n)
    for j in 1:n
        _add_affine_entry!(C, D, j, j, g.bpolys[j])
        if j >= 2
            ell, u = _block_offdiag_factors(g, j - 1, S)
            _add_affine_entry!(C, D, j - 1, j, ell)
            _add_affine_entry!(C, D, j, j - 1, u)
        end
    end
    return C, D
end

function roots(g::QuadraticThiele)
    if isempty(g.bpolys) || any(_qtcf_infinite_block, g.bpolys)
        return filter(isfinite, _companion_roots(g.P))
    end
    C, D = _qtcf_numerator_linear_pencil(g)
    z = try
        filter(isfinite, eigvals(-C, D))
    catch
        _, _, _, _, α, β = schur(complex(-C), complex(D))
        filter(isfinite, α ./ β)
    end
    return z
end

function Base.getproperty(g::QuadraticThiele, name::Symbol)
    name === :values && return getfield(g, :vals)
    name === :weights && return getfield(g, :bpolys)
    return getfield(g, name)
end

Base.propertynames(g::QuadraticThiele, private::Bool=false) =
    (fieldnames(typeof(g))..., :values, :weights)

function QuadraticThiele{T}(args...; kwargs...) where {T<:AbstractFloat}
    g = convert(T, QuadraticThiele(args...; kwargs...))
    all(isfinite, nodes(g)) || throw(ArgumentError("interpolation nodes must be finite"))
    all(isfinite, values(g)) || throw(ArgumentError("QTCF block target values must be finite"))
    all(b -> all(isfinite, b), weights(g)) ||
        throw(ArgumentError("QTCF block weights must be finite"))
    any(_same(g.nodes[i], g.nodes[j]) for i in eachindex(g.nodes)
        for j in firstindex(g.nodes):i-1) &&
        throw(ArgumentError("interpolation nodes must remain distinct after conversion"))
    return g
end

function _qtcf_block_sizes(block_sizes, n)
    sizes = block_sizes === nothing ? ones(Int, n) : collect(block_sizes)
    all(d -> d isa Integer && d in (1, 2), sizes) ||
        throw(ArgumentError("block_sizes must contain only 1 or 2"))
    sum(sizes; init=0) == n ||
        throw(DimensionMismatch("block_sizes must sum to the number of nodes"))
    return Int.(sizes)
end

function _qtcf_empty(::Type{S}) where {S}
    return QuadraticThiele(S[0], S[1], Tuple{S,S}[],
        Vector{S}[], Vector{S}[], S[], S[])
end

function _qtcf_append_block!(g::QuadraticThiele, z, y, b)
    S = eltype(g)
    zz, yy, bb = S.(z), S.(y), S.(b)
    all(isfinite, yy) || throw(ArgumentError("QTCF block target values must be finite"))
    any(isnan, bb) && throw(RFA.NaNException("Adding QTCF block caused NaN weight"))
    all(isfinite, bb) || throw(ArgumentError("QTCF block weights must be finite"))
    pair = (first(zz), last(zz))
    a = S.(_factor_poly(pair...))
    P, Q = _qtcf_continuants(vcat(g.bpolys, [bb]), vcat(g.apolys, [a]))
    resize!(g.P, length(P))
    copyto!(g.P, P)
    resize!(g.Q, length(Q))
    copyto!(g.Q, Q)
    push!(g.pairs, pair)
    push!(g.bpolys, bb)
    push!(g.apolys, a)
    append!(g.nodes, zz)
    append!(g.vals, yy)
    return g
end

function _new_weight_classic(g::QuadraticThiele, z_new, y_new)
    S = promote_type(eltype(g), typeof(z_new), typeof(y_new))
    u = S(y_new)
    for k in eachindex(g.bpolys)
        pair = g.pairs[k]
        a = _qtcf_eval_a(g.apolys[k], pair, z_new)
        b = _qtcf_eval_b(g.bpolys[k], pair, z_new)
        u = _qtcf_ratio(a, u - b)
    end
    return u
end

function _new_weight_onediv(g::QuadraticThiele, z_new, y_new)
    S = promote_type(eltype(g), typeof(z_new), typeof(y_new))
    return S(only(_reduced_values([z_new], [y_new], g.bpolys, g.apolys, g.pairs)))
end

_new_weight(::Classic, g::QuadraticThiele, z_new, y_new) =
    _new_weight_classic(g, z_new, y_new)
_new_weight(::OneDiv, g::QuadraticThiele, z_new, y_new) =
    _new_weight_onediv(g, z_new, y_new)

function add_node!(g::QuadraticThiele, z_new::Number, y_new::Number,
    method::ThieleMethod=default_weight_method())
    return add_node!(g, [z_new], [y_new], method)
end

function add_node!(g::QuadraticThiele, z_new::AbstractVector,
    y_new::AbstractVector, method::ThieleMethod=default_weight_method())
    length(z_new) == length(y_new) ||
        throw(DimensionMismatch("node and value counts must agree"))
    length(z_new) in (1, 2) ||
        throw(ArgumentError("a QTCF update requires one or two nodes"))
    isempty(g.bpolys) && (!isempty(g.nodes) || g.P != [0] || g.Q != [1]) &&
        throw(ArgumentError("node updates require a block representation"))
    z, y = eltype(g).(z_new), eltype(g).(y_new)
    all(isfinite, z) || throw(ArgumentError("interpolation nodes must be finite"))
    all(isfinite, y) || throw(ArgumentError("QTCF block target values must be finite"))
    any(t -> any(w -> _same(t, w), g.nodes), z) &&
        throw(ArgumentError("interpolation nodes must be distinct"))
    length(z) == 2 && _same(first(z), last(z)) &&
        throw(ArgumentError("a two-point block requires distinct nodes"))
    h = [_new_weight(method, g, z[k], y[k]) for k in eachindex(z)]
    any(isnan, h) && throw(RFA.NaNException("Adding QTCF block caused NaN reduced value"))
    all(isfinite, h) || throw(ArgumentError("QTCF block reduced values must be finite"))
    b = length(z) == 1 ? [only(h)] :
        _greedy_block_line(first(z), first(h), last(z), last(h))
    return _qtcf_append_block!(g, z, y, b)
end

function QuadraticThiele(x::AbstractVector, y::AbstractVector;
    block_sizes=nothing, method::ThieleMethod=default_weight_method())
    length(x) == length(y) ||
        throw(DimensionMismatch("node and value counts must agree"))
    sizes = _qtcf_block_sizes(block_sizes, length(x))
    if isempty(x)
        T = eltype(x) <: Number ? float(eltype(x)) : Float64
        U = eltype(y) <: Number ? float(eltype(y)) : Float64
        return _qtcf_empty(promote_type(T, U))
    end
    z, f = promote(float.(collect(x)), float.(collect(y)))
    g = _qtcf_empty(eltype(z))
    k = 1
    for d in sizes
        add_node!(g, z[k:k+d-1], f[k:k+d-1], method)
        k += d
    end
    return g
end

function QuadraticThiele(x::AbstractVector, y::AbstractVector,
    w::AbstractVector; block_sizes=nothing)
    length(x) == length(y) ||
        throw(DimensionMismatch("node and value counts must agree"))
    sizes = _qtcf_block_sizes(block_sizes, length(x))
    length(w) == length(sizes) ||
        throw(DimensionMismatch("provide one weight per block and its block_sizes"))
    if isempty(x)
        g = QuadraticThiele(x, y)
        W = eltype(w)
        W <: AbstractVector && (W = eltype(W))
        S = W <: Number ? promote_type(eltype(g), float(W)) : eltype(g)
        return _qtcf_empty(S)
    end
    bpolys = [b isa Number ? [float(b)] : float.(collect(b)) for b in w]
    all(1 <= length(b) <= d for (b, d) in zip(bpolys, sizes)) ||
        throw(ArgumentError("a block weight must be constant or affine for its block size"))
    S = promote_type(eltype(float.(x)), eltype(float.(y)), map(eltype, bpolys)...)
    z, f = S.(x), S.(y)
    all(isfinite, z) || throw(ArgumentError("interpolation nodes must be finite"))
    any(_same(z[i], z[j]) for i in eachindex(z) for j in firstindex(z):i-1) &&
        throw(ArgumentError("interpolation nodes must be distinct"))
    g = _qtcf_empty(S)
    k = 1
    for (d, b) in zip(sizes, bpolys)
        _qtcf_append_block!(g, z[k:k+d-1], f[k:k+d-1], b)
        k += d
    end
    return g
end

function _same_reduced_value(y1, y2; rtol=100eps(Float64))
    scale = max(1.0, Float64(abs(y1)), Float64(abs(y2)))
    return Float64(abs(y2 - y1)) <= rtol * scale
end

function _qtcf_checked_block(b, z1, h1, z2, h2)
    all(isfinite, b) && return b
    S = eltype(b)
    if !all(isfinite, (z1, h1, z2, h2))
        return S[NaN]
    end
    wide = setprecision(BigFloat, max(256, precision(BigFloat))) do
        if length(b) == 1
            [(big(h1) + big(h2)) / 2]
        else
            slope = (big(h2) - big(h1)) / (big(z2) - big(z1))
            [(big(h1) + big(h2) - slope * (big(z1) + big(z2))) / 2, slope]
        end
    end
    result = S.(wide)
    return all(isfinite, result) ? result : S[NaN]
end

function _line_poly(z1, y1, z2, y2)
    if !isfinite(y1) || !isfinite(y2)
        S = promote_type(typeof(y1), typeof(y2))
        return S[NaN]
    end
    if _same(z1, z2)
        return [y1]
    end
    if _same_reduced_value(y1, y2)
        S = promote_type(typeof(y1), typeof(y2))
        return _qtcf_checked_block(S[(S(y1) + S(y2)) / S(2)], z1, y1, z2, y2)
    end
    if _same(z2, conj(z1)) && _same(y2, conj(y1))
        R = promote_type(typeof(real(z1)), typeof(real(y1)))
        γ = R(imag(y1) / imag(z1))
        β = R(real(y1) - γ * real(z1))
        return _qtcf_checked_block(R[β, γ], z1, y1, z2, y2)
    end
    s = (y2 - y1) / (z2 - z1)
    zmid = (z1 + z2) / 2
    ymid = (y1 + y2) / 2
    return _qtcf_checked_block([ymid - s*zmid, s], z1, y1, z2, y2)
end




function _factor_poly(z1, z2)
    S = promote_type(ComplexF64, typeof(z1), typeof(z2))
    if _same(z1, z2)
        z = S(z1)
        return S[-z, one(S)]
    end
    z1, z2 = S(z1), S(z2)
    return S[z1*z2, -(z1 + z2), one(S)]
end

function _block_line(z1, h1, z2, h2)
    (isfinite(h1) && isfinite(h2)) || return [oftype(h1, NaN)]
    if _same(z1, z2)
        return [h1]
    end
    return _line_poly(z1, h1, z2, h2)
end

function _greedy_block_line(z1, h1, z2, h2)
    if !isfinite(h1) || !isfinite(h2)
        S = promote_type(typeof(h1), typeof(h2))
        return S[NaN]
    end
    _same(z1, z2) && return [h1]
    slope = (h2 - h1) / (z2 - z1)
    midpoint = (z1 + z2) / 2
    midpoint_value = (h1 + h2) / 2
    return _qtcf_checked_block([midpoint_value - slope * midpoint, slope], z1, h1, z2, h2)
end

function _block_nodes(z1, y1, z2, y2)
    if _same(z1, z2)
        return [z1], [y1]
    end
    return [z1, z2], [y1, y2]
end



function _canonical_param(d, s)
    L = convert(typeof(s), length(d))
    t = mod(s, L)
    tol = 100eps(typeof(t)) * max(one(t), abs(L))
    return (t <= tol || abs(t - L) <= tol) ? zero(t) : t
end

function _param_point(d, s)
    T = typeof(s)
    L = convert(T, length(d))
    if isclosed(d)
        return point(d, _canonical_param(d, s))
    end
    return point(d, min(max(s, zero(T)), L))
end

function _clean_param(d, s)
    T = typeof(s)
    L = convert(T, length(d))
    if isclosed(d)
        return _canonical_param(d, s)
    end
    return min(max(s, zero(T)), L)
end

function _reflection_param(d, s, reflection)
    T = typeof(s)
    L = convert(T, length(d))
    closed = isclosed(d)
    target = reflection(_param_point(d, s))
    _same(_param_point(d, zero(T)), target) && return zero(T)
    !closed && _same(_param_point(d, L), target) && return L
    if applicable(ComplexRegions.arg, d, target)
        return _clean_param(d, T(ComplexRegions.arg(d, target)))
    end
    if d isa ComplexRegions.AbstractPath
        distances = [ComplexRegions.dist(target, side) for side in d]
        k = argmin(distances)
        side = d[k]
        if applicable(ComplexRegions.arg, side, target)
            return _clean_param(d, T(k - 1) + T(ComplexRegions.arg(side, target)))
        end
    end
    ngrid = 768
    step = L / T(ngrid)
    best = zero(T)
    besterr = Inf
    for k in 0:ngrid
        t = T(k) * step
        tt = closed ? _canonical_param(d, t) : min(max(t, zero(T)), L)
        e = abs(point(d, tt) - target)
        if e < besterr
            best, besterr = t, e
        end
    end
    a, b = -step, step
    ϕ = (sqrt(T(5)) - one(T)) / T(2)
    c = b - ϕ * (b - a)
    dloc = a + ϕ * (b - a)
    dist(u) = begin
        tt = closed ? _canonical_param(d, best + u) : min(max(best + u, zero(T)), L)
        abs(point(d, tt) - target)
    end
    fc, fd = dist(c), dist(dloc)
    for _ in 1:60
        if fc > fd
            a = c; c = dloc; fc = fd; dloc = a + ϕ * (b - a); fd = dist(dloc)
        else
            b = dloc; dloc = c; fd = fc; c = b - ϕ * (b - a); fc = dist(c)
        end
    end
    return _clean_param(d, best + (a + b) / T(2))
end

function _reflection_partner_param(d, s, reflection)
    target = reflection(_param_point(d, s))
    sr = _reflection_param(d, s, reflection)
    return _same(_param_point(d, sr), target) ? sr : nothing
end

function _is_reflection_closed_domain(d, reflection; samples=24)
    T = typeof(float(length(d)))
    L = T(length(d))
    for k in 0:samples-1
        s = T(k) * L / T(samples)
        _reflection_partner_param(d, s, reflection) === nothing &&
            return false
    end
    return true
end

function _reflection_fixes_domain(d, reflection; samples=128)
    T = typeof(float(length(d)))
    L = T(length(d))
    for k in 0:samples
        z = _param_point(d, T(k) * L / T(samples))
        _same(z, reflection(z)) || return false
    end
    return true
end

function _reflection_fixes_grid(z, reflection)
    return all(node -> _same(node, reflection(node)), z)
end

function _with_reflection_params(d, params, s, reflection)
    T = typeof(s)
    L = convert(T, length(d))
    tol = T(100eps(Float64)) * max(one(T), abs(L))
    vals = T[zero(T), L]
    function add_break!(v)
        vv = _clean_param(d, T(v))
        (abs(vv) <= tol || abs(vv - L) <= tol) && return nothing
        push!(vals, vv)
        return nothing
    end
    foreach(add_break!, params)
    sr = _reflection_partner_param(d, s, reflection)
    add_break!(s)
    sr === nothing || add_break!(sr)
    sort!(vals)
    out = T[]
    for v in vals
        w = abs(v) <= tol ? zero(T) : abs(v - L) <= tol ? L : v
        if isempty(out) || abs(w - last(out)) > tol * max(one(T), abs(w), abs(last(out)))
            push!(out, w)
        end
    end
    first(out) != zero(T) && pushfirst!(out, zero(T))
    last(out) != L && push!(out, L)
    return out
end

function _with_pair_params(d, params, s1, s2)
    T = promote_type(eltype(params), typeof(s1), typeof(s2))
    L = T(length(d))
    tol = T(100eps(Float64)) * max(one(T), abs(L))
    vals = T[zero(T), L]
    function add_break!(s)
        value = _clean_param(d, T(s))
        (abs(value) <= tol || abs(value - L) <= tol) && return nothing
        push!(vals, value)
        return nothing
    end
    foreach(add_break!, params)
    add_break!(s1)
    add_break!(s2)
    sort!(vals)
    out = T[]
    for value in vals
        clean = abs(value) <= tol ? zero(T) : abs(value - L) <= tol ? L : value
        if isempty(out) ||
                abs(clean - last(out)) > tol * max(one(T), abs(clean), abs(last(out)))
            push!(out, clean)
        end
    end
    first(out) != zero(T) && pushfirst!(out, zero(T))
    last(out) != L && push!(out, L)
    return out
end



function _same(z, w)
    return abs(z - w) <= 100eps(Float64) * max(1, abs(z), abs(w))
end

function _node_hash_step(zgrid)
    scale = max(1.0, maximum(z -> Float64(abs(z)), zgrid; init=0.0))
    return 64eps(Float64) * scale
end

function _node_hash_key(z, h)
    return (round(Int, Float64(real(z)) / h), round(Int, Float64(imag(z)) / h))
end

function _node_hash_table(zgrid)
    h = _node_hash_step(zgrid)
    table = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in eachindex(zgrid)
        push!(get!(table, _node_hash_key(zgrid[i], h), Int[]), i)
    end
    return table, h
end

function _lookup_node_index(table, zgrid, z, h)
    k0 = _node_hash_key(z, h)
    for di in -1:1, dj in -1:1
        bucket = get(table, (k0[1] + di, k0[2] + dj), nothing)
        bucket === nothing && continue
        for idx in bucket
            _same(zgrid[idx], z) && return idx
        end
    end
    return nothing
end

function _insert_node_hash!(table, z, idx, h)
    push!(get!(table, _node_hash_key(z, h), Int[]), idx)
    return nothing
end

function _mark_reflected!(used, zgrid, z0, reflection)
    zr = reflection(z0)
    for i in eachindex(zgrid)
        if _same(zgrid[i], z0) || _same(zgrid[i], zr)
            used[i] = true
        end
    end
end

function _projective_reduce!(hnum, hden, a, b, zgrid, pair=nothing)
    @inbounds for k in eachindex(zgrid)
        N = hnum[k]
        D = hden[k]
        q = _qtcf_eval_a(a, pair, zgrid[k])
        bj = _qtcf_eval_b(b, pair, zgrid[k])
        hnum[k] = _qtcf_infinite_block(b) ? zero(N) : q * D
        hden[k] = _qtcf_infinite_block(b) ? -D : N - bj * D
        s = max(abs(hnum[k]), abs(hden[k]))
        if isfinite(s) && s > 0
            hnum[k] /= s
            hden[k] /= s
        end
    end
    return hnum, hden
end

function _admissible_reduced_score(score, hnum, hden, residual=score)
    out = copy(score)
    @inbounds for k in eachindex(out)
        if !isfinite(out[k]) || !isfinite(residual[k]) ||
                !isfinite(_qtcf_ratio(hnum[k], hden[k]))
            out[k] = -Inf
        end
    end
    return out
end

function _qtcf_linearized_residual(fvals, zvals, P, Q)
    scale = max(maximum(abs, P; init=0.0), maximum(abs, Q; init=0.0))
    (!isfinite(scale) || iszero(scale)) && return fill(-Inf, length(zvals))
    Ps = P ./ scale
    Qs = Q ./ scale
    score = Vector{Float64}(undef, length(zvals))
    @inbounds for k in eachindex(zvals)
        value = fvals[k] * _peval(Qs, zvals[k]) - _peval(Ps, zvals[k])
        score[k] = isfinite(value) ? Float64(abs(value)) : -Inf
    end
    return score
end

function _qtcf_projective_value(z, bpolys, apolys, pairs)
    m = length(bpolys)
    S = promote_type(ComplexF64, typeof(z), eltype(last(bpolys)))
    infinite = _qtcf_infinite_block(bpolys[m])
    numer = infinite ? one(S) : S(_qtcf_eval_b(bpolys[m], pairs[m], z))
    denom = infinite ? zero(S) : one(S)
    logscale = 0.0
    @inbounds for k in (m-1):-1:1
        b = S(_qtcf_eval_b(bpolys[k], pairs[k], z))
        a = S(_qtcf_eval_a(apolys[k], pairs[k], z))
        numer, denom = _qtcf_infinite_block(bpolys[k]) ?
            (numer, zero(S)) : (b * numer + a * denom, numer)
        scale = max(abs(numer), abs(denom))
        (!isfinite(scale) || iszero(scale)) && return S(NaN), S(NaN), NaN
        numer /= scale
        denom /= scale
        logscale += Float64(log(scale))
    end
    return numer, denom, logscale
end

function _qtcf_linearized_residual(fvals, zvals, bpolys, apolys, pairs)
    isempty(bpolys) && return ones(Float64, length(zvals))
    logscore = fill(-Inf, length(zvals))
    valid = falses(length(zvals))
    @inbounds for k in eachindex(zvals)
        numer, denom, logscale =
            _qtcf_projective_value(zvals[k], bpolys, apolys, pairs)
        cross = abs(fvals[k] * denom - numer)
        if isfinite(cross) && isfinite(logscale)
            valid[k] = true
            iszero(cross) || (logscore[k] = Float64(log(cross)) + logscale)
        end
    end
    any(valid) || return fill(-Inf, length(zvals))
    finite_logscore = filter(isfinite, logscore)
    if isempty(finite_logscore)
        score = fill(-Inf, length(zvals))
        score[valid] .= 0.0
        return score
    end
    offset = maximum(finite_logscore)
    score = fill(-Inf, length(zvals))
    @inbounds for k in eachindex(score)
        valid[k] && (score[k] = isfinite(logscore[k]) ?
            exp(logscore[k] - offset) : 0.0)
    end
    return score
end

function _qtcf_residual_rule(residual, name)
    residual in (:true_residual, :linearized_residual) || throw(ArgumentError(
        "$name must be :true_residual or :linearized_residual"))
    return residual
end

function _qtcf_residual_max(residual)
    value = 0.0
    found = false
    for r in residual
        r == -Inf && continue
        found = true
        isfinite(r) || return Inf
        value = max(value, Float64(r))
    end
    return found ? value : 0.0
end

function _reflected_continuum_samples(f, d, zgrid, fgrid, sgrid, reflection)
    S = promote_type(ComplexF64, eltype(zgrid), eltype(fgrid))
    T = promote_type(typeof(float(sgrid[firstindex(sgrid)])), typeof(float(length(d))))
    zout = S[]
    fout = S[]
    sout = T[]
    out_h = _node_hash_step(S.(zgrid))
    out_table = Dict{Tuple{Int,Int},Vector{Int}}()
    function add_sample!(z, y, s)
        zz = S(z)
        _lookup_node_index(out_table, zout, zz, out_h) !== nothing && return nothing
        push!(zout, zz)
        push!(fout, S(y))
        push!(sout, T(_clean_param(d, s)))
        _insert_node_hash!(out_table, zz, length(zout), out_h)
        return nothing
    end
    for k in eachindex(zgrid)
        z = S(zgrid[k])
        y = S(fgrid[k])
        s = T(sgrid[k])
        add_sample!(z, y, s)
        sr = _reflection_partner_param(d, s, reflection)
        if sr !== nothing
            zr = S(point(d, sr))
            add_sample!(zr, S(f(zr)), sr)
        end
    end
    return zout, fout, sout
end

function _qtcf_fixed_midpoint_params(d, full_interval, reflection)
    a, b = full_interval
    T = typeof(a)
    smid = _clean_param(d, (a + b) / T(2))
    z = _param_point(d, smid)
    return _same(z, reflection(z)) ? T[smid] : T[]
end

function _qtcf_fixed_params(d, full_interval, reflection)
    fixed = _qtcf_fixed_midpoint_params(d, full_interval, reflection)
    a, b = full_interval
    T = typeof(a)
    tol = T(100eps(Float64)) * max(one(T), abs(b - a))
    function add_fixed!(s)
        ss = _clean_param(d, T(s))
        z = _param_point(d, ss)
        _same(z, reflection(z)) || return nothing
        any(t -> abs(t - ss) <= tol, fixed) || push!(fixed, ss)
        return nothing
    end
    gap(s) = begin
        z = _param_point(d, _clean_param(d, T(s)))
        abs(reflection(z) - z)
    end
    nprobe = max(128, 32ceil(Int, Float64(b - a)))
    step = (b - a) / T(nprobe)
    probes = collect(range(a, b; length=nprobe + 1))
    gaps = gap.(probes)
    scale = maximum(z -> abs(_param_point(d, z)), probes; init=one(T))
    fixed_tol = T(1000eps(Float64)) * max(one(T), scale)
    if maximum(gaps; init=zero(T)) <= fixed_tol
        add_fixed!((a + b) / T(2))
        return fixed
    end
    for k in eachindex(probes)
        left_gap = k == firstindex(probes) ? Inf : gaps[k - 1]
        right_gap = k == lastindex(probes) ? Inf : gaps[k + 1]
        gaps[k] <= left_gap || continue
        gaps[k] <= right_gap || continue
        lo = probes[k] - step
        hi = probes[k] + step
        ϕ = (sqrt(T(5)) - one(T)) / T(2)
        c = hi - ϕ * (hi - lo)
        e = lo + ϕ * (hi - lo)
        fc, fe = gap(c), gap(e)
        for _ in 1:70
            if fc > fe
                lo = c
                c, fc = e, fe
                e = lo + ϕ * (hi - lo)
                fe = gap(e)
            else
                hi = e
                e, fe = c, fc
                c = hi - ϕ * (hi - lo)
                fc = gap(c)
            end
        end
        candidate = _clean_param(d, (lo + hi) / T(2))
        gap(candidate) <= fixed_tol && add_fixed!(candidate)
    end
    sort!(fixed)
    return fixed
end

function _qtcf_initial_mesh_params(d, full_interval, reflection)
    T = typeof(full_interval[1])
    a, b = full_interval
    vals = T[a, b]
    nseg = floor(Int, Float64(b))
    if nseg > 1 && abs(T(nseg) - b) <= sqrt(eps(T)) * max(one(T), abs(b))
        for k in 1:nseg-1
            vals = _with_reflection_params(d, vals, T(k), reflection)
        end
    end
    for s in _qtcf_fixed_params(d, full_interval, reflection)
        vals = _with_reflection_params(d, vals, s, reflection)
    end
    return vals
end

function _qtcf_direct_value_mode(y, yr)
    isempty(y) && return (:generic, zero(ComplexF64), Inf)
    S = promote_type(ComplexF64, eltype(y), eltype(yr))
    yy, yyR = S.(y), S.(yr)
    scale = max(1.0, maximum(abs, yy; init=0.0),
        maximum(abs, yyR; init=0.0))
    equality_defect = maximum(abs, yyR .- yy; init=0.0) / scale
    center = sum(yy .+ yyR) / S(2 * length(yy))
    odd_defect = maximum(abs, yyR .- (S(2) .* center .- yy); init=0.0) / scale
    tolerance = 1000eps(Float64)
    equality_defect <= tolerance && return (:equality, center, equality_defect)
    odd_defect <= tolerance && return (:odd, center, odd_defect)
    return (:generic, center, min(equality_defect, odd_defect))
end

function _reduced_values(ztest, ftest, bpolys, apolys, pairs=nothing)
    S = promote_type(ComplexF64, eltype(ztest), eltype(ftest))
    hnum = S.(ftest)
    hden = ones(S, length(ftest))
    for k in eachindex(hnum)
        if isinf(hnum[k]) && !isnan(hnum[k])
            hnum[k], hden[k] = one(S), zero(S)
        end
    end
    for j in eachindex(bpolys)
        pair = pairs === nothing || j > length(pairs) ? nothing : pairs[j]
        _projective_reduce!(hnum, hden, apolys[j], bpolys[j], ztest, pair)
    end
    return _qtcf_ratio.(hnum, hden)
end

function _qtcf_converged(record, threshold, allowed)
    isfinite(record.error) && record.error <= threshold || return false
    allowed === true && return true
    if ismissing(record.poles)
        record.poles = poles(record.interpolant)
    end
    return all(allowed, record.poles)
end

function _qtcf_result(history, reason, allowed, iterations)
    assessed = findall(record -> !isnan(record.error), history)
    best = reason === :converged ? lastindex(history) :
        assessed[RFA.best_acceptable(history[assessed], allowed)]
    return copy(history[best].interpolant),
        RFA.ConvergenceStatus(reason, best, iterations, history[best].error)
end

function _quadratic_from_grid(zgrid, fgrid, scalar_seed; tol=1000eps(Float64), max_iter=80,
    stagnation=5, allowed=true, reflection=conj,
    value_reflection=nothing, noise_floor_factor=1000,
    float_type::Type=promote_type(RFA.real_type(eltype(zgrid)), RFA.real_type(eltype(fgrid)), Float64))

    S = complex(float_type)
    zgrid = S.(zgrid)
    fgrid = S.(fgrid)
    used = falses(length(zgrid))
    reflected_index = Vector{Int}(undef, length(zgrid))
    reflected_table, reflected_h = _node_hash_table(zgrid)
    for i in eachindex(zgrid)
        j = _lookup_node_index(reflected_table, zgrid, reflection(zgrid[i]), reflected_h)
        reflected_index[i] = j === nothing ? 0 : j
    end
    if value_reflection === nothing
        missing_index = findfirst(i -> reflected_index[i] == 0 &&
            !_same(zgrid[i], reflection(zgrid[i])), eachindex(zgrid))
        missing_index === nothing || throw(ArgumentError(
            "A distinct reflected mate lies outside the discrete grid. " *
            "Provide value_reflection=(z, y, z_reflected) -> f(z_reflected)."))
    end
    fmax = maximum(abs, fgrid)
    stop_tol = max(tol * fmax, noise_floor_factor * eps(float_type) * max(1, fmax))

    Pold = S[1]
    Qold = S[0]
    P = S[]
    Q = S[]
    aprev = S[]
    pairs = Tuple{S,S}[]
    bpolys = Vector{S}[]
    apolys = Vector{S}[]
    zact = S[]
    yact = S[]
    hnum = copy(fgrid)
    hden = ones(S, length(fgrid))
    reflected_index_at(idx) = reflected_index[idx] == 0 ? nothing : reflected_index[idx]

    function initial_grid_index()
        z0 = S(scalar_seed)
        idx = _lookup_node_index(reflected_table, zgrid, z0, reflected_h)
        idx !== nothing && return idx
        _, nearest = findmin(abs.(zgrid .- z0))
        return _same(zgrid[nearest], z0) ? nearest : nothing
    end

    function choose_reflected_cached(score)
        best = 0
        best_score = -Inf
        for i in eachindex(zgrid)
            used[i] && continue
            isfinite(score[i]) || continue
            ip = reflected_index[i]
            block_score = Float64(score[i])
            if ip != 0 && !used[ip] && isfinite(score[ip])
                block_score = max(block_score, Float64(score[ip]))
            end
            if block_score > best_score
                best, best_score = i, block_score
            end
        end
        best == 0 && return nothing
        return best
    end

    function accept_grid_block!(idx)
        (idx === nothing || idx == 0) && return false
        znew = zgrid[idx]
        ynew = fgrid[idx]
        hnew = isempty(bpolys) ? ynew : _qtcf_ratio(hnum[idx], hden[idx])
        (isfinite(ynew) && isfinite(hnew)) || return false
        idxr = reflected_index_at(idx)
        zc = idxr === nothing ? S(reflection(znew)) : zgrid[idxr]
        yc = idxr === nothing ? S(value_reflection(znew, ynew, zc)) : fgrid[idxr]
        hc = if isempty(bpolys)
            yc
        elseif idxr === nothing
            _reduced_values(S[zc], S[yc], bpolys, apolys, pairs)[1]
        else
            _qtcf_ratio(hnum[idxr], hden[idxr])
        end
        (isfinite(yc) && isfinite(hc)) || return false
        b = S.(_block_line(znew, hnew, zc, hc))
        all(isfinite, b) || return false
        anew = S.(_factor_poly(znew, zc))
        if isempty(bpolys)
            P = _qtcf_infinite_block(b) ? S[1] : copy(b)
            Q = _qtcf_infinite_block(b) ? S[0] : S[1]
            if _qtcf_infinite_block(b)
                Pold, Qold = S[0], S[0]
            end
        else
            P, Q, Pold, Qold = _qtcf_update_continuants(b, aprev, P, Q, Pold, Qold)
        end
        aprev = copy(anew)
        push!(pairs, (znew, zc))
        push!(bpolys, copy(b))
        push!(apolys, copy(anew))
        zb, yb = _block_nodes(znew, ynew, zc, yc)
        append!(zact, S.(zb))
        append!(yact, S.(yb))
        used[idx] = true
        idxr === nothing ? _mark_reflected!(used, zgrid, znew, reflection) : (used[idxr] = true)
        _projective_reduce!(hnum, hden, anew, b, zgrid, (znew, zc))
        return true
    end

    i0 = initial_grid_index()
    i0 === nothing && error("Initial QTCF node matching the scalar TCF start was not found in the supplied grid.")
    accept_grid_block!(i0) ||
        error("Initial QTCF node matching the scalar TCF start was not admissible.")





    g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
    err = abs.(fgrid .- g.(zgrid))
    err[used] .= -Inf
    errmax = _qtcf_residual_max(err)
    history = [RFA.IterationRecord(g, errmax, missing)]
    used_history = [copy(used)]
    stop_reason = :max_degree
    best_err = errmax

    blocks = length(bpolys)
    last_improve = blocks
    while blocks < max_iter && !_qtcf_converged(history[end], stop_tol, allowed)
        score = _admissible_reduced_score(err, hnum, hden)
        idx = choose_reflected_cached(score)
        if idx === nothing
            stop_reason = all(used) ? :exhausted : :node_failure
            break
        end
        while idx !== nothing && !accept_grid_block!(idx)
            score[idx] = -Inf
            idxr = reflected_index_at(idx)
            idxr !== nothing && (score[idxr] = -Inf)
            idx = choose_reflected_cached(score)
        end
        if idx === nothing
            stop_reason = all(used) ? :exhausted : :node_failure
            break
        end

        g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
        err = abs.(fgrid .- g.(zgrid))
        err[used] .= -Inf
        errmax = _qtcf_residual_max(err)
        push!(history, RFA.IterationRecord(g, errmax, missing))
        push!(used_history, copy(used))
        new_blocks = length(bpolys)
        if isfinite(errmax) && errmax < best_err
            errmax < 0.999 * best_err && (last_improve = new_blocks)
            best_err = errmax
        end
        blocks = new_blocks
        if blocks - last_improve >= stagnation
            stop_reason = :stagnated
            break
        end
        if length(history) >= stagnation + 5
            recent = [h.error for h in history[end-stagnation+1:end]]
            if minimum(recent) > 5best_err
                stop_reason = :stagnated
                break
            end
        end
    end
    _qtcf_converged(history[end], stop_tol, allowed) && (stop_reason = :converged)
    best_g, stop = _qtcf_result(history, stop_reason, allowed, length(bpolys))
    return best_g, used_history[stop.best], history, stop
end




function _quadratic_greedy_from_grid(zgrid, fgrid, scalar_seed;
    tol=1000eps(Float64), max_iter=80, stagnation=5, allowed=true, noise_floor_factor=1000,
    first_residual=:true_residual, second_residual=:true_residual,
    float_type::Type=promote_type(RFA.real_type(eltype(zgrid)), RFA.real_type(eltype(fgrid)), Float64))

    first_residual = _qtcf_residual_rule(first_residual, "first_residual")
    second_residual = _qtcf_residual_rule(second_residual, "second_residual")
    S = complex(float_type)
    zgrid = S.(zgrid)
    fgrid = S.(fgrid)
    used = falses(length(zgrid))
    table, hash_step = _node_hash_table(zgrid)
    fmax = maximum(abs, fgrid)
    stop_tol = max(tol * fmax,
        noise_floor_factor * eps(float_type) * max(1, fmax))

    Pold = S[1]
    Qold = S[0]
    P = S[]
    Q = S[]
    aprev = S[]
    pairs = Tuple{S,S}[]
    bpolys = Vector{S}[]
    apolys = Vector{S}[]
    zact = S[]
    yact = S[]
    hnum = copy(fgrid)
    hden = ones(S, length(fgrid))

    function current_reduced(idx)
        isempty(bpolys) && return fgrid[idx]
        return _qtcf_ratio(hnum[idx], hden[idx])
    end

    function accept_pair!(idx1, idx2)
        z1, y1 = zgrid[idx1], fgrid[idx1]
        z2, y2 = zgrid[idx2], fgrid[idx2]
        h1, h2 = current_reduced(idx1), current_reduced(idx2)
        all(isfinite, (y1, y2, h1, h2)) || return false
        b = S.(_greedy_block_line(z1, h1, z2, h2))
        all(isfinite, b) || return false
        a = S.(_factor_poly(z1, z2))
        if isempty(bpolys)
            P = _qtcf_infinite_block(b) ? S[1] : copy(b)
            Q = _qtcf_infinite_block(b) ? S[0] : S[1]
            if _qtcf_infinite_block(b)
                Pold, Qold = S[0], S[0]
            end
        else
            P, Q, Pold, Qold = _qtcf_update_continuants(b, aprev, P, Q, Pold, Qold)
        end
        aprev = copy(a)
        push!(pairs, (z1, z2))
        push!(bpolys, copy(b))
        push!(apolys, copy(a))
        zb, yb = _block_nodes(z1, y1, z2, y2)
        append!(zact, S.(zb))
        append!(yact, S.(yb))
        used[idx1] = true
        used[idx2] = true
        _projective_reduce!(hnum, hden, a, b, zgrid, (z1, z2))
        return true
    end

    function choose(score, skip=0)
        best = 0
        best_score = -Inf
        for i in eachindex(score)
            (used[i] || i == skip) && continue
            isfinite(score[i]) || continue
            if score[i] > best_score
                best, best_score = i, Float64(score[i])
            end
        end
        return best == 0 ? nothing : best
    end

    seed = _lookup_node_index(table, zgrid, S(scalar_seed), hash_step)
    seed === nothing && error("Initial QTCF node matching the scalar TCF start was not found.")
    accept_pair!(seed, seed) ||
        error("Initial QTCF node matching the scalar TCF start was not admissible.")

    g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
    err = abs.(fgrid .- g.(zgrid))
    err[used] .= -Inf
    errmax = _qtcf_residual_max(err)
    history = [RFA.IterationRecord(g, errmax, missing)]
    used_history = [copy(used)]
    stop_reason = :max_degree
    best_err = errmax
    blocks = length(bpolys)
    last_improve = blocks

    while blocks < max_iter && !_qtcf_converged(history[end], stop_tol, allowed)
        true_score = _admissible_reduced_score(err, hnum, hden)
        linearized_score = if first_residual === :linearized_residual ||
                second_residual === :linearized_residual
            _admissible_reduced_score(
                _qtcf_linearized_residual(
                    fgrid, zgrid, bpolys, apolys, pairs), hnum, hden, err)
        else
            nothing
        end
        first_score = copy(first_residual === :linearized_residual ?
            linearized_score : true_score)
        idx1 = choose(first_score)
        accepted = false
        while idx1 !== nothing
            second_score = copy(second_residual === :linearized_residual ?
                linearized_score : true_score)
            for i in eachindex(second_score)
                _same(zgrid[i], zgrid[idx1]) && (second_score[i] = -Inf)
            end
            idx2 = choose(second_score)
            while idx2 !== nothing
                if accept_pair!(idx1, idx2)
                    accepted = true
                    break
                end
                second_score[idx2] = -Inf
                idx2 = choose(second_score, idx1)
            end
            accepted && break
            first_score[idx1] = -Inf
            idx1 = choose(first_score)
        end
        if !accepted
            stop_reason = all(used) ? :exhausted : :node_failure
            break
        end

        g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
        err = abs.(fgrid .- g.(zgrid))
        err[used] .= -Inf
        errmax = _qtcf_residual_max(err)
        push!(history, RFA.IterationRecord(g, errmax, missing))
        push!(used_history, copy(used))
        blocks = length(bpolys)
        if isfinite(errmax) && errmax < best_err
            errmax < 0.999 * best_err && (last_improve = blocks)
            best_err = errmax
        end
        if blocks - last_improve >= stagnation
            stop_reason = :stagnated
            break
        end
        if length(history) >= stagnation + 5
            recent = [h.error for h in history[end-stagnation+1:end]]
            if minimum(recent) > 5best_err
                stop_reason = :stagnated
                break
            end
        end
    end
    _qtcf_converged(history[end], stop_tol, allowed) && (stop_reason = :converged)
    best_g, stop = _qtcf_result(history, stop_reason, allowed, length(bpolys))
    return best_g, used_history[stop.best], history, stop
end



function approximate(
    y::AbstractVector{T}, z::AbstractVector{S}, ::QuadraticThiele;
    float_type::Type = promote_type(RFA.real_type(eltype(z)), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = length(y),
    stagnation::Int = 5,
    reflection::Union{Nothing,Function} = nothing,
    value_reflection::Union{Nothing,Function} = nothing,
    noise_floor_factor::Real = 1000,
    first_residual::Symbol = :true_residual,
    second_residual::Symbol = :true_residual
    ) where {T<:Number,S<:Number}

    if reflection === nothing
        _, idx_min = findmin(abs, y)
        g, used, history, stop = _quadratic_greedy_from_grid(z, y, z[idx_min];
            tol=convert(float_type, tol), max_iter, stagnation, allowed, noise_floor_factor,
            first_residual, second_residual, float_type)
        return RFA.DiscreteApproximation(y, z, g, .!used, allowed, history, stop)
    end

    (first_residual === :true_residual && second_residual === :true_residual) ||
        throw(ArgumentError(
            "residual selection options apply only to QTCF without a reflection map"))

    _, idx_min = findmin(abs, y)
    scalar_seed = z[idx_min]
    g, used, history, stop = _quadratic_from_grid(z, y, scalar_seed;
        tol=convert(float_type, tol), max_iter, stagnation, allowed,
        reflection, value_reflection, noise_floor_factor, float_type)
    return RFA.DiscreteApproximation(y, z, g, .!used, allowed, history, stop)
end




function _quadratic_greedy_continuum(f, d;
    float_type, tol, allowed, max_iter, refinement, initial_refinement, stagnation,
    first_residual, second_residual)

    first_residual = _qtcf_residual_rule(first_residual, "first_residual")
    second_residual = _qtcf_residual_rule(second_residual, "second_residual")
    max_iter = something(max_iter, 80)
    Tpar = typeof(float(length(d)))
    L = Tpar(length(d))
    num_ref = max(refinement, initial_refinement)
    s_nodes = Tpar[zero(Tpar), L]
    maxpoints = max(50, 2max_iter + 20)
    path = RFA.DiscretizedPath(
        d, s_nodes; refinement=num_ref, maxpoints=maxpoints)
    S = complex(float_type)

    Pold = S[1]
    Qold = S[0]
    P = S[]
    Q = S[]
    aprev = S[]
    pairs = Tuple{S,S}[]
    bpolys = Vector{S}[]
    apolys = Vector{S}[]
    zact = S[]
    yact = S[]
    active_h = _node_hash_step(S[point(d, zero(Tpar)), point(d, L)])
    active_table = Dict{Tuple{Int,Int},Vector{Int}}()
    fmax = zero(float_type)

    active_node_seen(z) =
        _lookup_node_index(active_table, zact, S(z), active_h) !== nothing

    function active_node_push!(z, y)
        push!(zact, S(z))
        push!(yact, S(y))
        _insert_node_hash!(active_table, S(z), length(zact), active_h)
        return nothing
    end

    function audit_samples(audit_path, breaks, nref)
        indices = CartesianIndices((1:length(breaks)-1, 1:nref+1))
        z = S.(vec(audit_path.points[indices]))
        s = Tpar.(vec(audit_path.params[indices]))
        push!(z, S(point(d, L)))
        push!(s, L)
        return z, S.(f.(z)), s
    end

    function block_data(s1raw, s2raw)
        s1, s2 = _clean_param(d, Tpar(s1raw)), _clean_param(d, Tpar(s2raw))
        z1, z2 = S(_param_point(d, s1)), S(_param_point(d, s2))
        y1, y2 = S(f(z1)), S(f(z2))
        h1 = isempty(bpolys) ? y1 :
            _reduced_values(S[z1], S[y1], bpolys, apolys, pairs)[1]
        h2 = isempty(bpolys) ? y2 :
            _reduced_values(S[z2], S[y2], bpolys, apolys, pairs)[1]
        all(isfinite, (y1, y2, h1, h2)) || return nothing
        b = S.(_greedy_block_line(z1, h1, z2, h2))
        all(isfinite, b) || return nothing
        a = S.(_factor_poly(z1, z2))
        return (; s1, s2, z1, z2, y1, y2, b, a)
    end

    function accept_pair!(s1, s2)
        data = block_data(s1, s2)
        data === nothing && return false
        if isempty(bpolys)
            P = _qtcf_infinite_block(data.b) ? S[1] : copy(data.b)
            Q = _qtcf_infinite_block(data.b) ? S[0] : S[1]
            if _qtcf_infinite_block(data.b)
                Pold, Qold = S[0], S[0]
            end
        else
            P, Q, Pold, Qold = _qtcf_update_continuants(data.b, aprev, P, Q, Pold, Qold)
        end
        aprev = copy(data.a)
        push!(pairs, (data.z1, data.z2))
        push!(bpolys, copy(data.b))
        push!(apolys, copy(data.a))
        zb, yb = _block_nodes(data.z1, data.y1, data.z2, data.y2)
        for i in eachindex(zb)
            active_node_push!(zb[i], yb[i])
        end
        fmax = max(fmax, abs(data.y1), abs(data.y2))
        s_nodes = _with_pair_params(d, s_nodes, data.s1, data.s2)
        return true
    end

    function choose(score)
        _, index = findmax(score)
        return isfinite(score[index]) ? index : nothing
    end

    accept_pair!(zero(Tpar), zero(Tpar)) ||
        error("Initial QTCF continuum node was not admissible.")
    g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
    history = [RFA.IterationRecord(g, NaN, missing)]
    best_err = Inf
    stop_reason = :max_degree
    blocks = length(bpolys)
    last_improve = blocks

    RFA.reset!(path, s_nodes; refinement=num_ref)
    ztest, ftest, stest = audit_samples(path, s_nodes, num_ref)
    fmax = max(fmax, maximum(abs, ftest))

    while true
        used = active_node_seen.(ztest)
        err = abs.(ftest .- g.(ztest))
        err[used] .= -Inf
        errmax = _qtcf_residual_max(err)
        reduced = _reduced_values(ztest, ftest, bpolys, apolys, pairs)
        true_score = copy(err)
        ineligible = .!isfinite.(err) .| .!isfinite.(reduced)
        true_score[ineligible] .= -Inf
        linearized_score = if first_residual === :linearized_residual ||
                second_residual === :linearized_residual
            score = _qtcf_linearized_residual(
                ftest, ztest, bpolys, apolys, pairs)
            score[ineligible] .= -Inf
            score
        else
            nothing
        end
        first_score = copy(first_residual === :linearized_residual ?
            linearized_score : true_score)

        history[end].error = errmax
        if errmax < best_err
            errmax < 0.999 * best_err && (last_improve = blocks)
            best_err = errmax
        end
        if _qtcf_converged(history[end], tol * fmax, allowed)
            stop_reason = :converged
            break
        end
        blocks >= max_iter && break

        accepted = false
        next_num_ref = num_ref
        idx1 = choose(first_score)
        while idx1 !== nothing
            s1 = stest[idx1]
            second_score = copy(second_residual === :linearized_residual ?
                linearized_score : true_score)
            for i in eachindex(second_score)
                _same(ztest[i], ztest[idx1]) && (second_score[i] = -Inf)
            end
            idx2 = choose(second_score)
            while idx2 !== nothing
                if accept_pair!(s1, stest[idx2])
                    accepted = true
                    next_num_ref = num_ref > refinement ?
                        max(refinement, num_ref - 1) : num_ref
                    break
                end
                second_score[idx2] = -Inf
                idx2 = choose(second_score)
            end
            accepted && break
            first_score[idx1] = -Inf
            idx1 = choose(first_score)
        end
        if !accepted
            stop_reason = :node_failure
            break
        end

        g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
        push!(history, RFA.IterationRecord(g, NaN, missing))
        blocks = length(bpolys)
        if blocks - last_improve >= stagnation
            stop_reason = :stagnated
            break
        end
        num_ref = next_num_ref
        try
            RFA.reset!(path, s_nodes; refinement=num_ref)
        catch
            stop_reason = :refinement
            break
        end
        ztest, ftest, stest = audit_samples(path, s_nodes, num_ref)
        fmax = max(fmax, maximum(abs, ftest))
        if length(history) >= stagnation + 5
            recent = [h.error for h in history[end-stagnation:end-1] if isfinite(h.error)]
            if !isempty(recent) && minimum(recent) > 5best_err
                stop_reason = :stagnated
                break
            end
        end
    end
    best_g, stop = _qtcf_result(history, stop_reason, allowed, length(bpolys))
    return RFA.ContinuumApproximation(f, d, best_g, allowed, path, history, stop)
end


function approximate(
    f::Function, d::ComplexCurveOrPath, ::QuadraticThiele;
    float_type::Type = promote_type(RFA.real_type(d), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed = true,
    max_degree::Int = 100,
    max_iter::Int = max_degree,
    refinement::Int = 3,
    initial_refinement::Int = 15,
    stagnation::Int = 5,
    reflection::Union{Nothing,Function} = nothing,
    value_reflection::Union{Nothing,Function} = nothing,
    first_residual::Symbol = :true_residual,
    second_residual::Symbol = :true_residual,
    qtcf_strategy = SymmetricQTCF
    )

    if allowed == :strict
        allowed = z -> ComplexRegions.dist(z, d) > tol
    end
    reflection === nothing && value_reflection !== nothing && throw(ArgumentError(
        "value_reflection requires a reflection map"))

    strategy = _qtcf_strategy_symbol(qtcf_strategy)
    if !(strategy === SymmetricQTCF || strategy === :symmetric_qtcf)
        error("Unknown QTCF strategy $(repr(qtcf_strategy)). Use SymmetricQTCF.")
    end
    if reflection === nothing
        return _quadratic_greedy_continuum(f, d;
            float_type, tol, allowed, max_iter, refinement, initial_refinement, stagnation,
            first_residual, second_residual)
    end
    (first_residual === :true_residual && second_residual === :true_residual) ||
        throw(ArgumentError(
            "residual selection options apply only to QTCF without a reflection map"))
    max_iter = something(max_iter, 80)

    Tpar = typeof(float(length(d)))
    full_interval = (zero(Tpar), Tpar(length(d)))
    num_ref = max(refinement, initial_refinement)
    s_nodes = _qtcf_initial_mesh_params(d, full_interval, reflection)
    path = RFA.DiscretizedPath(d, s_nodes; refinement=num_ref,
        maxpoints=max(50, 2max_iter * num_ref + length(s_nodes) + 20))
    τ = path.points
    idx_test = CartesianIndices((1:length(s_nodes)-1, 1:num_ref+1))
    S = complex(float_type)
    fτ = Matrix{S}(undef, size(τ))
    fτ[idx_test] .= f.(τ[idx_test])
    fmax = maximum(abs, view(fτ, idx_test))

    Pold = S[1]
    Qold = S[0]
    P = S[]
    Q = S[]
    aprev = S[]
    pairs = Tuple{S,S}[]
    bpolys = Vector{S}[]
    apolys = Vector{S}[]
    zact = S[]
    yact = S[]
    parity_fixed_index = 0
    parity_center = zero(S)
    active_h = _node_hash_step(vec(τ[idx_test]))
    active_table = Dict{Tuple{Int,Int},Vector{Int}}()

    probe_values = S[]
    reflected_probe_values = S[]
    for I in idx_test
        sprobe = path.params[I]
        zprobe = S(_param_point(d, sprobe))
        zrprobe = S(reflection(zprobe))
        _same(zprobe, zrprobe) && continue
        yprobe = S(f(zprobe))
        yrprobe = S(f(zrprobe))
        if value_reflection !== nothing
            ymapped = S(value_reflection(zprobe, yprobe, zrprobe))
            isapprox(ymapped, yrprobe; rtol=1e-10, atol=1e-12) || throw(ArgumentError(
                "value_reflection does not match f at reflected continuum points"))
            yrprobe = ymapped
        end
        push!(probe_values, yprobe)
        push!(reflected_probe_values, yrprobe)
        length(probe_values) == 16 && break
    end
    value_mode, value_center, _ =
        _qtcf_direct_value_mode(probe_values, reflected_probe_values)
    value_equality_mode = value_mode === :equality

    active_node_seen(z) = _lookup_node_index(active_table, zact, z, active_h) !== nothing
    function active_node_push!(z, y)
        push!(zact, S(z))
        push!(yact, S(y))
        _insert_node_hash!(active_table, S(z), length(zact), active_h)
        return nothing
    end

    function accept_continuum_block!(sraw)
        s1 = _clean_param(d, sraw)
        z1 = S(_param_point(d, s1))
        reflected_s2 = _reflection_partner_param(d, s1, reflection)
        z2 = S(reflection(z1))
        y1 = S(f(z1))
        y2 = if _same(z1, z2)
            y1
        elseif value_reflection === nothing
            S(f(z2))
        else
            S(value_reflection(z1, y1, z2))
        end
        h1 = isempty(bpolys) ? y1 : _reduced_values(S[z1], S[y1], bpolys, apolys, pairs)[1]
        h2 = _same(z1, z2) ? h1 :
            (isempty(bpolys) ? y2 : _reduced_values(S[z2], S[y2], bpolys, apolys, pairs)[1])
        all(isfinite, (y1, y2, h1, h2)) || return false
        fmax = reflected_s2 === nothing ? max(fmax, abs(y1)) :
            max(fmax, abs(y1), abs(y2))
        fixed_block = _same(z1, z2)
        b = S.(_block_line(z1, h1, z2, h2))
        if parity_fixed_index > 0 && !fixed_block && isfinite(h1) && isfinite(h2) &&
                _same((z1 + z2) / S(2), parity_center) &&
                _same_reduced_value(h1, -h2)
            slope = (h2 - h1) / (z2 - z1)
            midpoint = (z1 + z2) / S(2)
            b = _qtcf_checked_block(S[-slope * midpoint, slope], z1, h1, z2, h2)
        end
        all(isfinite, b) || return false
        anew = S.(_factor_poly(z1, z2))
        if isempty(bpolys)
            P = _qtcf_infinite_block(b) ? S[1] : copy(b)
            Q = _qtcf_infinite_block(b) ? S[0] : S[1]
            if _qtcf_infinite_block(b)
                Pold, Qold = S[0], S[0]
            end
        else
            P, Q, Pold, Qold = _qtcf_update_continuants(b, aprev, P, Q, Pold, Qold)
        end
        aprev = copy(anew)
        if fixed_block && value_equality_mode && parity_fixed_index == 0
            parity_fixed_index = length(bpolys) + 1
            parity_center = z1
        end
        push!(pairs, (S(z1), S(z2)))
        push!(bpolys, copy(b))
        push!(apolys, copy(anew))
        zb, yb = _block_nodes(z1, y1, z2, y2)
        for i in eachindex(zb)
            active_node_push!(zb[i], yb[i])
        end
        s_nodes = _with_reflection_params(d, s_nodes, s1, reflection)
        return true
    end

    parity_complete() = !value_equality_mode || parity_fixed_index == 0 ||
        length(bpolys) <= parity_fixed_index || iseven(length(bpolys) - parity_fixed_index)

    function choose_continuum_candidate(choice)
        _, jmax = findmax(choice)
        return isfinite(choice[jmax]) ? jmax : nothing
    end

    initial_param = zero(Tpar)
    if value_mode === :odd
        fixed_params = _qtcf_fixed_params(d, full_interval, reflection)
        if !isempty(fixed_params)
            _, fixed_index = findmin(s -> abs(S(f(_param_point(d, s))) - value_center),
                fixed_params)
            initial_param = fixed_params[fixed_index]
        end
    end
    accept_continuum_block!(initial_param) ||
        error("Initial QTCF continuum node was not admissible.")
    g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
    history = [RFA.IterationRecord(g, NaN, missing)]
    best_err = Inf
    stop_reason = :max_degree
    last_improve = 0

    RFA.reset!(path, s_nodes; refinement=num_ref)
    τ = path.points
    idx_test = CartesianIndices((1:length(s_nodes)-1, 1:num_ref+1))
    fτ = Matrix{S}(undef, size(τ))
    fτ[idx_test] .= f.(τ[idx_test])
    fmax = max(fmax, maximum(abs, view(fτ, idx_test)))

    blocks = length(bpolys)
    while true
        zraw = vec(τ[idx_test])
        fraw = vec(fτ[idx_test])
        sraw = vec(path.params[idx_test])
        ztest, ftest, stest = _reflected_continuum_samples(f, d, zraw, fraw, sraw, reflection)
        used = falses(length(ztest))
        for k in eachindex(ztest)
            active_node_seen(ztest[k]) && (used[k] = true)
        end
        err = abs.(ftest .- g.(ztest))
        err[used] .= -Inf
        choice = copy(err)
        choice[.!isfinite.(choice)] .= -Inf
        htest = _reduced_values(ztest, ftest, bpolys, apolys, pairs)
        choice[.!isfinite.(htest)] .= -Inf
        err_max = _qtcf_residual_max(err)
        complete = parity_complete()
        if complete
            history[end].error = err_max
            if err_max < best_err
                err_max < 0.999 * best_err && (last_improve = blocks)
                best_err = err_max
            end
        end
        if complete && _qtcf_converged(history[end], tol * fmax, allowed)
            stop_reason = :converged
            break
        end
        blocks >= max_iter && break

        accepted = false
        while true
            jmax = choose_continuum_candidate(choice)
            jmax === nothing && break
            snew1 = stest[jmax]
            znew1 = ztest[jmax]
            znew2 = S(reflection(znew1))
            if accept_continuum_block!(snew1)
                accepted = true
                break
            end
            for k in eachindex(choice)
                if _same(ztest[k], znew1) || _same(ztest[k], znew2)
                    choice[k] = -Inf
                end
            end
        end
        if !accepted
            stop_reason = :node_failure
            break
        end
        g = QuadraticThiele(P, Q, pairs, bpolys, apolys, zact, yact)
        parity_complete() && push!(history, RFA.IterationRecord(g, NaN, missing))

        if num_ref > refinement
            num_ref = max(refinement, num_ref - 1)
        end
        try
            RFA.reset!(path, s_nodes; refinement=num_ref)
        catch
            stop_reason = :refinement
            break
        end
        blocks = length(bpolys)
        if parity_complete() && blocks - last_improve >= stagnation
            stop_reason = :stagnated
            break
        end
        τ = path.points
        idx_test = CartesianIndices((1:length(s_nodes)-1, 1:num_ref+1))
        fτ = Matrix{S}(undef, size(τ))
        fτ[idx_test] .= f.(τ[idx_test])
        fmax = max(fmax, maximum(abs, view(fτ, idx_test)))
        if length(history) >= stagnation + 5
            recent = [h.error for h in history[end-stagnation:end-1] if isfinite(h.error)]
            if !isempty(recent) && minimum(recent) > 5best_err
                stop_reason = :stagnated
                break
            end
        end
    end

    best_g, stop = _qtcf_result(history, stop_reason, allowed, length(bpolys))
    return RFA.ContinuumApproximation(f, d, best_g, allowed, path, history, stop)
end

function approximate(f::Function, d::Symmetric; kw...)
    if haskey(kw, :method)
        method, rest = RFA._pop_deprecated_method_kw(kw)
        return approximate(f, d, method; rest...)
    end
    return approximate(f, d, QuadraticThiele(); kw...)
end

function approximate(f::Function, d::Symmetric, ::QuadraticThiele; kw...)
    return approximate(f, d.domain, QuadraticThiele();
        reflection=d.reflection, kw...)
end

function approximate(f::Function, d::Symmetric, method::RFA.AbstractRationalFunction; kw...)
    throw(ArgumentError("Symmetric dispatch requires QuadraticThiele/QTCF"))
end

function approximate(y::AbstractVector, d::Symmetric; kw...)
    if haskey(kw, :method)
        method, rest = RFA._pop_deprecated_method_kw(kw)
        return approximate(y, d, method; rest...)
    end
    return approximate(y, d, QuadraticThiele(); kw...)
end

function approximate(y::AbstractVector, d::Symmetric, ::QuadraticThiele; kw...)
    d.domain isa AbstractVector ||
        throw(ArgumentError("sampled values require a vector of domain points"))
    return approximate(y, d.domain, QuadraticThiele(); reflection=d.reflection, kw...)
end

function approximate(y::AbstractVector, d::Symmetric, method::RFA.AbstractRationalFunction; kw...)
    throw(ArgumentError("Symmetric dispatch requires QuadraticThiele/QTCF"))
end
