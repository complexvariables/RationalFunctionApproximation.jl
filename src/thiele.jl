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
function Base.convert(::Type{Thiele{T}}, r::Thiele{S}) where {S,T}
    return Thiele{T}( T.(nodes(r)), T.(values(r)), T.(weights(r)) )
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

# Evaluation at a point
function evaluate(r::Thiele, z::Number)
    return if isinf(z)
        if isodd(length(r.nodes))
            sum(r.weights[1:2:end])
        else
            Inf
        end
    elseif isnan(z)
        NaN
    else
        _evaluator(r, z)
    end
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

_evaluator = _evaluate_onediv    # default choice

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
        end
        r.weights[1] * a + b, a
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

derivative(r::Thiele, ζ::Number) = derivative(r)(ζ)
function derivative(r::Thiele)
    return function(ζ)
        p, q, pʹ, qʹ = _evaluate_numden_derivs(r, ζ)
        return (pʹ * q - p * qʹ) / (q^2)
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
        if abs(q) < 100eps(T) * abs(p) && !iszero(qʹ)
            # simple pole
            res[i] = p / qʹ
        else
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
function add_node!(r::Thiele, z_new, y_new)
    w = _new_weight(r.nodes, r.weights, z_new, y_new)
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
        u = (z_new - z[k]) / (u - w[k])
    end
    return u
end

_new_weight = _new_weight_onediv   # default choice

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

    # Arrays of test points have one row per node (except the last)
    τ = path.points
    idx_test = CartesianIndices((1:1, 2:num_ref+1))
    idx_new_test = idx_test
    fτ = Matrix{eltype(fσ)}(undef, size(τ))        # f at test points
    fτ[idx_test] .= f.(τ[idx_test])
    fmax = maximum(abs, view(fτ, idx_test))        # scale of f

    # Initialize rational approximation
    r = Thiele(σ, fσ)
    history = [IterationRecord(r, NaN, missing)]

    # Main iteration
    n = 1       # iteration counter
    while true
        # test_actual = view(fτ, idx_test)      # f at test points
        err = @. abs(fτ[idx_test] - r(τ[idx_test]))
        err_max, idx_max = findmax(err)
        history[n].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping with estimated error $(round(history[status].error, sigdigits=4)) after $n iterations")
            r = history[status].interpolant
        end
        (status != 0) && break

        # Add node to approximant
        idx_new = idx_test[idx_max]      # location of worst test point
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

        # Add node to path
        try
            idx_new_test = add_node!(path, idx_new)
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
            idx_test = CartesianIndices((1:n, 2:num_ref+1))
            @. fτ[idx_test] = f(τ[idx_test])
            fmax = maximum(abs, view(fτ, idx_test))
        else
            # At new test points only, evaluate f.
            @. fτ[idx_new_test] = f(τ[idx_new_test])
            idx_test = CartesianIndices((1:n, 2:num_ref+1))
        end
    end
    return ContinuumApproximation(f, d, r, allowed, path, history)
end

function approximate(::Type{Thiele},
    y::AbstractVector{T}, z::AbstractVector{S};
    float_type::Type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::AbstractFloat = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = 240,
    stagnation::Int = 5,
    ) where {T<:Number,S<:Number}

    m = length(z)
    idx_test = trues(m)
    fmax = maximum(abs, y)     # scale of f
    _, idx_min = findmin(abs, y)
    r = Thiele([z[idx_min]], [y[idx_min]])
    idx_test[idx_min] = false

    history = [IterationRecord(r, NaN, missing)]
    n = 1    # iteration counter
    while length(z) > 0
        idx_max, err_max = 0, -Inf
        for i in eachindex(y)
            if idx_test[i]
                err = abs(y[i] - r(z[i]))
                if err > err_max
                    err_max = err
                    idx_max = i
                end
            end
        end
        history[n].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping with estimated error $(round(history[status].error, sigdigits=4)) after $n iterations")
            r = history[status].interpolant
        end
        (status != 0) && break

        # Add new node:
        try
            add_node!(r, z[idx_max], y[idx_max])
            push!(history, IterationRecord(r, NaN, missing))
            idx_test[idx_max] = false
        catch(e)
            # look for the best acceptable case
            status = quitting_check(history, stagnation, tol, fmax, 1, allowed)
            r = history[status].interpolant
            @warn("NaN weight encountered; stopping with estimated error $(round(history[status].error, sigdigits=4))")
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
