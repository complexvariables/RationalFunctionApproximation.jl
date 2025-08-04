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

# convenience accessors
nodes(r::Thiele) = r.nodes
Base.values(r::Thiele) = r.values
weights(r::Thiele) = r.weights
function degrees(r::Thiele)
    n = length(r.nodes) - 1
    q, s = divrem(n, 2)
    return q + s, q
end
degree(r::Thiele) = div(length(r.nodes) - 1, 2)

function Base.copy(r::Thiele)
    return Thiele(copy(r.nodes), copy(r.values), copy(r.weights))
end

function evaluate(r::Thiele, z::Number)
    return if isinf(z)
         if isodd(n)
            sum(r.weights[1:2:end])
        else
            Inf
        end
    elseif isnan(z)
        NaN
    else
        numer, denom = _evaluate(r, z)
        numer / denom
    end
end

function _evaluate(r::Thiele, z::Number)
    @assert isfinite(z)
    n = length(r.weights)
    return if n == 1
        r.weights[1], 1
    else
        # use 3-term pair recurrence to avoid division until the end
        a = r.weights[n]
        b = z - r.nodes[n-1]
        @fastmath for k in n-1:-1:2
            t = r.weights[k] * a + b
            b = a * (z - r.nodes[k-1])
            a = t
        end
        r.weights[1] * a + b, a
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
    # This is a low-accuracy version for speed.
    ζ = poles(r)
    res = similar( complex(ζ) )
    T = eltype(res)
    circ = [cispi(T(2k // 9)) for k in 0:8]
    for (i, z) in enumerate(ζ)
        δ = minimum(abs, ζ[[1:i-1;i+1:end]] .- z) / 2
        res[i] = (δ / 9) * sum(u * r(z + δ*u) for u in circ)
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
function add_node!(r::Thiele, new_σ, new_f)
    σ = r.nodes
    push!(σ, new_σ)
    push!(r.values, new_f)

    # Compute the new weight
    φ = r.weights
    push!(φ, new_f)
    n = length(σ)
    for k in 1:n-1
        d = φ[n] - φ[k]
        @inbounds if iszero(d)
            φ[n] = Inf
        else
            φ[n] = (σ[n] - σ[k]) / d
        end
        isnan(φ[n]) && @error "Inverse difference produced NaN." σ φ
    end
    return r
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
        history[end].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping at estimated error $(round(err_max, sigdigits=4)) after $n iterations")
            r = history[status].interpolant
        end
        (status != 0) && break

        ### Refinement
        idx_new = idx_test[idx_max]      # location of worst test point
        add_node!(r, τ[idx_new], fτ[idx_new])
        push!(history, IterationRecord(r, NaN, missing))
        n += 1

        idx_new_test = add_node!(path, idx_new)
        # In the initial phase, we throw out the old test points.
        if num_ref > refinement    # initial phase
            num_ref -= 3    # gradually decrease refinement level
            s = first(collect(path))
            path = DiscretizedPath(d, s; refinement=num_ref, maxpoints=max_iter * refinement)
            τ = path.points
            idx_test = CartesianIndices((1:n, 2:num_ref+1))
            @. fτ[idx_test] = f(τ[idx_test])
            fmax = maximum(abs, view(fτ, idx_test))
        else
            # At new test points only, evaluate f.
            @. fτ[idx_new_test] = f(τ[idx_new_test])
            idx_test = CartesianIndices((1:n, 2:num_ref+1))
        end
    end
    return Approximation(f, d, r, allowed, path, history)
end

function approximate(::Type{Thiele},
    y::AbstractVector{T}, z::AbstractVector{S};
    float_type::Type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::AbstractFloat = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = 240,
    stagnation::Int = 5,
    ) where {T<:Number,S<:Number}

    y = copy(y)
    z = collect(copy(z))
    fmax = maximum(abs, y)     # scale of f
    _, idx_min = findmin(abs, y)
    r = Thiele([z[idx_min]], [y[idx_min]])
    deleteat!(y, idx_min)    # remove the minimum value
    deleteat!(z, idx_min)    # remove the minimum value
    history = [IterationRecord(r, NaN, missing)]
    n = 1    # iteration counter
    while length(z) > 0
        err = @. abs(y - r(z))
        err_max, idx_new = findmax(err)
        history[n].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping at estimated error $(round(err_max, sigdigits=4)) after $n iterations")
            r = history[n].interpolant
        end
        (status != 0) && break

        # Add new node:
        add_node!(r, z[idx_new], y[idx_new])
        deleteat!(y, idx_new)    # remove the minimum value
        deleteat!(z, idx_new)    # remove the minimum value
        push!(history, IterationRecord(r, NaN, missing))
        n += 1
    end
    return r, history
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
