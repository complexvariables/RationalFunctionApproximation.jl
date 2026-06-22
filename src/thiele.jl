struct Thiele{T,S} <: AbstractRationalInterpolant{T,S}
    nodes::Vector{S}
    values::Vector{S}
    weights::Vector{S}
    function Thiele{T}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S}
        ) where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight)
        new{T,S}(node, value, weight)
    end
    function Thiele{T,S}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S}
        ) where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight)
        new{T,S}(node, value, weight)
    end
end

# Convert the numeric type:
function Base.convert(::Type{F}, r::Thiele{T,S}) where {F<:AbstractFloat,T,S<:Real}
    return Thiele{F,F}( F.(nodes(r)), F.(values(r)), F.(weights(r)))
end

function Base.convert(::Type{F}, r::Thiele{T,S}) where {F<:AbstractFloat,T,S<:Complex}
    CF = complex(F)
    return Thiele{F,CF}( CF.(nodes(r)), CF.(values(r)), CF.(weights(r)))
end

# convenience accessors amd overloads
nodes(r::Thiele) = r.nodes
Base.values(r::Thiele) = r.values
weights(r::Thiele) = r.weights
function degrees(r::Thiele)
    n = length(r.nodes) - 1
    q, s = divrem(n, 2)
    return q + s, q
end
degree(r::Thiele) = div(length(r.nodes) - 1, 2)

Base.isreal(r::Thiele) = isreal(r.nodes) && isreal(r.weights)

function Base.copy(r::Thiele)
    return Thiele(copy(r.nodes), copy(r.values), copy(r.weights))
end

const TCF = Thiele    # alias

"""
    ThieleMethod

Abstract supertype for the recurrence strategies a [`Thiele`](@ref) can use, namely
[`OneDiv`](@ref) and [`Classic`](@ref). Pass an instance to `evaluate(r, z, method)` or
`add_node!(r, z, y, method)` for a one-off choice, or set the global default with
[`set_eval_method`](@ref) / [`set_weight_method`](@ref).

These let external or test code A/B the two algorithms. The default selectors return
singletons, so they constant-fold and the production paths stay static, allocation-free calls.
"""
abstract type ThieleMethod end

"""
    OneDiv() <: ThieleMethod

Evaluate a [`Thiele`](@ref) via the numerator/denominator pair recurrence, deferring the
single division to the end. Faster and the default for point evaluation; falls back to
[`Classic`](@ref) if the denominator underflows to zero.
"""
struct OneDiv  <: ThieleMethod end

"""
    Classic() <: ThieleMethod

Evaluate a [`Thiele`](@ref) by the textbook backward continued-fraction recurrence (one
division per node). Slower than [`OneDiv`](@ref) but the default for new-weight computation in
`add_node!`.
"""
struct Classic <: ThieleMethod end

default_eval_method()   = OneDiv()    # point evaluation
default_weight_method() = Classic()   # new-weight computation in `add_node!`

"""
    set_eval_method(m::ThieleMethod)
    set_weight_method(m::ThieleMethod)

Globally select the recurrence used for Thiele point evaluation / new-weight computation
(`OneDiv()` or `Classic()`). Intended for testing; each call recompiles the affected path.
"""
set_eval_method(m::ThieleMethod)   = (@eval default_eval_method()   = $m; m)
set_weight_method(m::ThieleMethod) = (@eval default_weight_method() = $m; m)

# Evaluation at a point
function evaluate(r::Thiele, z::Number, method::ThieleMethod=default_eval_method())
    return if isinf(z)
        if isodd(length(r.nodes))
            sum(r.weights[1:2:end])
        else
            Inf
        end
    elseif isnan(z)
        NaN
    else
        _evaluate(r, z, method)
    end
end

# faster for array-valued evaluation
function evaluate(r::Thiele, z::AbstractArray{<:Number})
    t = similar(z, eltype(r.values))
    return evaluate!(t, r, z)
end

function evaluate!(t::AbstractArray, r::Thiele, z::AbstractArray{<:Number})
    return evaluate!(t, r, z, similar(t), similar(t))
end

# In-place evaluation with caller-supplied scratch `a`, `b` (same shape as `t`),
# so a hot loop can reuse them instead of allocating on every call.
function evaluate!(t::AbstractArray, r::Thiele, z::AbstractArray{<:Number},
    a::AbstractArray, b::AbstractArray)
    # use 3-term pair recurrence to avoid division until the end
    n = length(r.weights)
    if n == 1
        t .= r.weights[1]
    else
        a .= r.weights[n]
        b .= z .- r.nodes[n-1]
        @inbounds for k in n-1:-1:2
            t .= b
            axpy!(r.weights[k], a, t)
            b .= a .* (z .- r.nodes[k-1])
            a .= t
        end
        t .= r.weights[1] .+ b ./ a
    end
    return t
end

function _evaluate_classic(r::Thiele, z::Number)
    n = length(r.nodes)
    u = last(r.weights)
    @inbounds for k in n-1:-1:1
        u = r.weights[k] + (z - r.nodes[k]) / u
    end
    return u
end

function _evaluate_onediv(r::Thiele, z::Number)
    numer, denom = _evaluate_numden(r, z)
    return if iszero(denom)
        @debug "Evaluation produced a division by zero at " z
        _evaluate_classic(r, z)  # fallback to slow evaluation
    else
        numer / denom
    end
end

_evaluate(r::Thiele, z, ::OneDiv)  = _evaluate_onediv(r, z)
_evaluate(r::Thiele, z, ::Classic) = _evaluate_classic(r, z)

function _evaluate_numden(r::Thiele, z::Number)
    # use 3-term pair recurrence to avoid division until the end
    @assert isfinite(z)
    n = length(r.weights)
    return if n == 1
        r.weights[1], 1
    else
        a = r.weights[n]
        b = z - r.nodes[n-1]
        @inbounds for k in n-1:-1:2
            t = r.weights[k] * a + b
            b = a * (z - r.nodes[k-1])
            a = t
            if abs2(a) < 1e-200   # prevent underflow
                @debug "scaling to avoid underflow at " z
                a *= 1e100
                b *= 1e100
            end
        end
        r.weights[1] * a + b, a
    end
end

function _derivative!(ζ, w, z, order, A, B, vals)
    n = length(w)
    if n == 1
        vals[1] = w[1]
        vals[2:order+1] .= 0
        return vals
    end
    A[1] = w[n]
    A[2:order+1] .= 0
    B[1] = ζ - z[n-1]
    B[2] = 1
    B[3:order+1] .= 0
    @inbounds for k in n-1:-1:2
        d = ζ - z[k-1]
        for m in order:-1:0
            t = w[k] * A[m+1] + B[m+1]
            B[m+1] = d * A[m+1] + (m > 0 ? m * A[m] : 0)
            A[m+1] = t
        end
    end
    # numer derivs = w[1] * A + B
    # denom derivs = A
    vals[1] = w[1]  + B[1] / A[1]
    for m in 1:order
        s = sum( binomial(m, j) * vals[j+1] * A[m-j+1] for j in 0:m-1 )
        vals[m+1] = (w[1] * A[m+1] + B[m+1] - s) / A[1]
    end
    return vals
end

function derivative(r::Thiele{T,S}, order::Integer=1) where {T,S}
    f = derivative(r, [order])
    return ζ -> only(f(ζ))
end

function derivative(r::Thiele{T,S}, orders::AbstractVector{<:Integer}) where {T,S}
    w, z = weights(r), nodes(r)
    order = maximum(orders)
     A = similar(complex(z), order+1)
    B = similar(A)
    vals = similar(A)
    index = 1 .+ orders
    return function(ζ)
        d = _derivative!(ζ, w, z, order, A, B, vals)
        if isreal(ζ) && isreal(r)
            d = real(d)
        end
        return d[index]
    end
end


function _evaluate_numden_derivs(r::Thiele, z::Number)
    n = length(r.weights)
    return if n == 1
        r.weights[1], 1, 0, 1
    else
        a = r.weights[n]
        b = z - r.nodes[n-1]
        aʹ = 0
        bʹ = 1
        @inbounds for k in n-1:-1:2
            d = z - r.nodes[k-1]
            t = r.weights[k] * aʹ + bʹ
            bʹ = a + d * aʹ
            aʹ = t
            t = r.weights[k] * a + b
            b = a * d
            a = t
        end
        r.weights[1] * a + b, a, r.weights[1] * aʹ + bʹ, aʹ
    end
end

function poles(r::Thiele{T,S}) where {T,S}
    n = length(r.nodes)
    (n < 3) && return S[]
    C = diagm(1 => -ones(T, n-2))
    xi = view(r.nodes, 2:n-1)
    D = diagm(0 => view(r.weights, 2:n), -1 => -ones(T, n-2), 1 => -xi)
    z = try
        filter(isfinite, eigvals(D, C))
    catch
        _, _, _, _, ⍺, β = schur(complex(D), complex(C))
        filter(isfinite, ⍺ ./ β)
    end
    return z
end

function roots(r::Thiele{S,T}) where {S,T}
    n = length(r.nodes)
    C = diagm(1 => -ones(T, n-1))
    xi = view(r.nodes, 1:n-1)
    D = diagm(0 => r.weights, -1 => -ones(T, n-1), 1 => -xi)
    return filter(isfinite, eigvals(D, C))
end

function residues(r::Thiele)
    ζ = poles(r)
    res = similar(complex(ζ))
    T = real_type(eltype(ζ))
    for i in eachindex(ζ)
        # is it a simple pole?
        p, q, pʹ, qʹ = _evaluate_numden_derivs(r, ζ[i])
        if !iszero(qʹ)
            # simple pole
            res[i] = p / qʹ
        else
            @debug "Fallback for residue at pole " ζ[i]
            # not a simple pole; use fallback contour integral
            res[i] = Res(r, ζ[i]; avoid=ζ)
        end
    end
    return ζ, res
end

function Thiele(nodes::AbstractVector, values::AbstractVector, weights::AbstractVector)
    if isempty(nodes) && eltype(nodes) == Any
        nodes = values = weights = Float64[]
    end
    nodes, values, weights = promote(nodes, values, weights)
    T = eltype(values)
    return Thiele{real_type(T)}(nodes, values, weights)
end

function Thiele(x::AbstractVector, y::AbstractVector)
    idx = axes(x, 1)
    d = copy(y)
    for (len, i) in enumerate(idx)
        for k in first(idx, len-1)
            d[i] = (x[i] - x[k]) / (d[i] - d[k])
        end
    end
    return Thiele(collect(x), y, d)
end

# update in-place the nodes and weights
function add_node!(r::Thiele, z_new, y_new, method::ThieleMethod=default_weight_method())
    w = _new_weight(method, r.nodes, r.weights, z_new, y_new)
    if isnan(w)
        throw(NaNException("Adding node at $z_new caused NaN weight"))
        @debug("Adding node at $z_new caused NaN weight")
    end
    push!(r.weights, w)
    push!(r.nodes, z_new)
    push!(r.values, y_new)
    return r
end

function _new_weight_onediv(z, w, z_new, y_new)
    a = 1
    b = y_new
    @inbounds for k in eachindex(z)
        t = -w[k] * a + b
        b = a * (z_new - z[k])
        a = t
    end
    return b / a
end

function _new_weight_classic(z, w, z_new, y_new)
    u = y_new
    @inbounds for k in eachindex(z)
        D = u - w[k]
        N = z_new - z[k]
        u =
        if iszero(D)
            if iszero(N)
                convert(typeof(D), NaN)
            else
                convert(typeof(D), Inf)
            end
        else
            N / D
        end
    end
    return u
end

_new_weight(::OneDiv,  z, w, z_new, y_new) = _new_weight_onediv(z, w, z_new, y_new)
_new_weight(::Classic, z, w, z_new, y_new) = _new_weight_classic(z, w, z_new, y_new)

# Per-iteration kernels for the continuum `approximate` loop below, gathered under one name.
# They are kept as separate functions so each specializes on the concrete element types of its
# array arguments: even though `DiscretizedPath.points` is parameterized, ComplexRegions' curve
# evaluation (`point(d, ·)` / `path(t)`) infers as `Any`, so `σ`, `fσ`, and everything derived
# from them (`fτ`, `rbuf`, `zbuf`) have non-inferable element types inside `approximate`. Inline
# per-element loops over them would dispatch dynamically on every element; routing through these
# barriers replaces that with one dispatch per call plus a type-stable loop.

# Evaluate r at the active test points (reusing the scratch buffers) and return the maximum
# error and the index (into `active`) of the worst point; mirrors `findmax` (NaN wins).
function _sweep!(rbuf, zbuf, abuf, bbuf, r::Thiele, τ, fτ, active)
    m = length(active)
    resize!(rbuf, m); resize!(zbuf, m); resize!(abuf, m); resize!(bbuf, m)
    @inbounds for i in 1:m
        zbuf[i] = τ[active[i]]
    end
    evaluate!(rbuf, r, zbuf, abuf, bbuf)
    @inbounds for i in 1:m       # array evaluation skips the underflow check
        isnan(rbuf[i]) && (rbuf[i] = r(zbuf[i]))
    end
    err_max = abs(fτ[active[1]] - rbuf[1])
    kmax = 1
    @inbounds for k in 2:m
        e = abs(fτ[active[k]] - rbuf[k])
        if isnan(e) || e > err_max
            err_max, kmax = e, k
            isnan(e) && break
        end
    end
    return err_max, kmax
end

# Sample f at the points indexed by `idx`, storing into `fτ`; return max |f| over them.
function _sweep!(fτ, f::Function, τ, idx)
    fmax = abs(zero(eltype(fτ)))
    @inbounds for i in idx
        v = f(τ[i])
        fτ[i] = v
        a = abs(v)
        (a > fmax) && (fmax = a)
    end
    return fmax
end

# TODO: This should probably enforce parameters S and T
approximate(::Type{Thiele{S,T}}, args...; kw...) where {S,T} = approximate(Thiele, args...; kw...)

function approximate(::Type{Thiele},
    f::Function, d::Union{ComplexPath,ComplexCurve};
    float_type::Type = promote_type(real_type(d), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed::Union{Function,Bool} = z -> dist(z, d) > tol,
    max_iter::Int = 240,
    refinement::Int = 3,
    stagnation::Int = 5
    )

    num_ref = 15    # initial number of test points between nodes; decreases to `refinement`
    path = DiscretizedPath(d, [0, 1]; refinement=num_ref, maxpoints=max_iter * refinement)
    σ = [point(d, 0)]
    fσ = f.(σ)         # f at nodes

    # Test points have one matrix row per node (except the last). `active` holds *linear*
    # indices of the currently live test points; it need not be rectangular. We grow it
    # incrementally (each accepted node appends one row's worth of test points) and only
    # rebuild it wholesale when the path is re-discretized in the initial phase. Linear
    # integer indexing is much cheaper here than CartesianIndex indexing, which dominated
    # the per-point loops below via `to_indices`.
    # `path.points` is concretely typed at runtime, but inference here can't prove it (its
    # element type flows from `point(d, ·)`, which infers as `Any`), so the per-point work
    # goes through the kernel barriers above rather than relying on the type in this scope.
    τ = path.points
    R = size(τ, 1)        # rows; linear index of (i, j) is (j-1)*R + i
    active = Int[]
    sizehint!(active, max_iter * refinement)
    for j in 2:num_ref+1
        push!(active, (j - 1) * R + 1)       # row 1, columns 2:num_ref+1
    end
    fτ = Matrix{eltype(fσ)}(undef, size(τ))        # f at test points
    fmax = _sweep!(fτ, f, τ, active)        # sample f; fmax is the scale of f

    # Reusable scratch, aligned with `active`: the gathered test points `zbuf`, r there
    # (`rbuf`), and the two work arrays the recurrence in `evaluate!` would otherwise
    # allocate on every call. Gathering into a dense `zbuf` lets `evaluate!` run on
    # contiguous vectors. Use the natural recurrence type (values promoted with the
    # nodes/points) so real problems stay real-valued instead of paying for complex
    # arithmetic, while complex domains (e.g. an imaginary interval) still get a wide
    # enough buffer.
    S = promote_type(eltype(fσ), eltype(σ))
    zbuf = Vector{eltype(σ)}()
    rbuf = S[]
    abuf = S[]
    bbuf = S[]

    # Initialize rational approximation
    r = Thiele(σ, fσ)
    history = [IterationRecord(r, NaN, missing)]

    # Main iteration
    n = 1       # iteration counter
    while true
        err_max, idx_max = _sweep!(rbuf, zbuf, abuf, bbuf, r, τ, fτ, active)
        history[n].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping with estimated error $(round(history[status].error, sigdigits=4)) after $n iterations")
            r = history[status].interpolant
        end
        (status != 0) && break

        # Add node to approximant
        idx_new = active[idx_max]      # location of worst test point
        try
            add_node!(r, τ[idx_new], fτ[idx_new])
            push!(history, IterationRecord(r, NaN, missing))
        catch(e)
            # look for the best acceptable case
            status = quitting_check(history, stagnation, tol, fmax, 1, allowed)
            r = history[status].interpolant
            @warn("NaN weight encountered; stopping with estimated error $(round(history[status].error, sigdigits=4))")
            @debug("Error $e")
            break
        end

        # Add node to path (needs the (row, column) form of the worst point)
        local idx_new_test
        try
            idx_new_test = add_node!(path, CartesianIndices(τ)[idx_new])
        catch
            # look for the best acceptable case
            status = quitting_check(history, stagnation, tol, fmax, 1, allowed)
            r = history[status].interpolant
            @warn("Maximum path refinement exceeded; stopping with estimated error $(round(history[status].error, sigdigits=4))")
            break
        end

        n += 1
        # In the initial phase, we throw out the old test points.
        if num_ref > refinement    # initial phase
            num_ref = max(refinement, num_ref - 1)    # gradually decrease refinement level
            s = first(collect(path))
            reset!(path, s; refinement=num_ref)
            empty!(active)
            for j in 2:num_ref+1, i in 1:n
                push!(active, (j - 1) * R + i)
            end
            fmax = _sweep!(fτ, f, τ, active)
        else
            # Evaluate f at the new test points only, and extend `active` by the new row.
            _sweep!(fτ, f, τ, idx_new_test)
            new_row = idx_new_test[end][1]
            for j in 2:num_ref+1
                push!(active, (j - 1) * R + new_row)
            end
        end
    end
    return ContinuumApproximation(f, d, r, allowed, path, history)
end

function approximate(::Type{Thiele},
    y::AbstractVector{T}, z::AbstractVector{S};
    float_type::Type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::AbstractFloat = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = length(y),
    stagnation::Int = 5,
    ) where {T<:Number,S<:Number}

    m = length(z)
    idx_test = trues(m)
    fmax = maximum(abs, y)     # scale of f
    _, idx_min = findmin(abs, y)
    r = Thiele([z[idx_min]], [y[idx_min]])
    y_test = copy(collect(y))
    z_test = copy(collect(z))
    deleteat!(y_test, idx_min)
    deleteat!(z_test, idx_min)
    r_test = similar(y_test)

    history = [IterationRecord(r, NaN, missing)]
    n = 1    # iteration counter
    while length(z) > 0
        evaluate!(r_test, r, z_test)
        r_test .-= y_test
        err_max, idx_max = findmax(abs(e) for e in r_test)
        history[n].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            if isinf(err_max)
                @warn("Used all sample values without convergence")
                status = max_iter
            else
                @warn("Stopping with estimated error $(round(history[status].error, sigdigits=4)) after $n iterations")
            end
            r = history[status].interpolant
        end
        (status != 0) && break

        # Add new node:
        try
            add_node!(r, z_test[idx_max], y_test[idx_max])
            push!(history, IterationRecord(r, NaN, missing))
            deleteat!(y_test, idx_max)
            deleteat!(z_test, idx_max)
            deleteat!(r_test, idx_max)
        catch(e)
            # look for the best acceptable case
            status = quitting_check(history, stagnation, tol, fmax, 1, allowed)
            r = history[status].interpolant
            @warn("Adding node failed; stopping with estimated error $(round(history[status].error, sigdigits=4))")
            @debug("Error $e")
            break
        end
        n += 1
    end
    return DiscreteApproximation(y, z, r, idx_test, allowed, history)
end

# Operations with scalars that can be done quickly.

function Base.:+(r::Thiele, s::Number)
    w = copy(r.weights)
    w = [w[1] + s; w[2:end]]
    return Thiele(r.nodes, r.values .+ s, w)
end

function Base.:-(r::Thiele)
    return Thiele(r.nodes, -r.values, -r.weights)
end

function Base.:*(r::Thiele, s::Number)
    return if iszero(s)
        zer = zero(eltype(r.values))
        Thiele(r.nodes[1:1], [zer], [zer])
    else
        y = s * r.values
        w = Vector{eltype(y)}(undef, length(r.weights))
        w[1:2:end] .= r.weights[1:2:end] * s
        w[2:2:end] .= r.weights[2:2:end] / s
        Thiele(r.nodes, y, w)
    end
end
