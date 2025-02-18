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

nodes(r::Thiele) = r.nodes
Base.values(r::Thiele) = r.values
weights(r::Thiele) = r.weights
function degrees(r::Thiele)
    n = length(r.nodes) - 1
    q, s = divrem(n, 2)
    return q + s, q
end
degree(r::Thiele) = div(length(r.nodes) - 1, 2)

function evaluate(r::Thiele, z::Number)
    n = length(r.nodes)
    u = last(r.weights)
    for k in n-1:-1:1
        u = r.weights[k] + (z - r.nodes[k]) / u
    end
    return u
end

function poles(r::Thiele{T}) where {T}
    n = length(r.nodes)
    (n < 3) && return T[]
    C = diagm(1 => -ones(T, n-2))
    xi = view(r.nodes, 2:n-1)
    D = diagm(0 => view(r.weights, 2:n), -1 => -ones(T, n-2), 1 => -xi)
    z = try
        filter(isfinite, eigvals(D, C))
    catch
        @warn "Problem with generalized eigvals"
        T[]
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

function Thiele(nodes::AbstractVector, values::AbstractVector, weights::AbstractVector)
    if isempty(nodes) && eltype(nodes) == Any
        nodes = values = weights = Float64[]
    end
    T = eltype(values)
    return Thiele{real_type(T)}(nodes, values, weights)
end

function Thiele(x::AbstractVector, y::AbstractVector)
    idx = axes(x, 1)
    @assert idx == axes(y, 1)
    d = copy(y)
    for (len, i) in enumerate(idx)
        for k in first(idx, len-1)
            d[i] = (x[i] - x[k]) / (d[i] - d[k])
        end
    end
    return Thiele(x, y, d)
end

function update_test_values!(::Type{Thiele}, numeric_type::Type, num_refine::Integer, max_degree::Integer)
    Δ = Array{numeric_type}(undef, num_refine, max_degree+1, max_degree+3)
    return Δ
end

function update_test_values!(r::Thiele, Δ, τ, fτ, idx_new_test)
    σ, φ = r.nodes, r.weights
    n = length(σ)

    @inbounds @fastmath for i in idx_new_test, j in eachindex(σ)
        Δ[i, j] = τ[i] - σ[j]
    end

    # Evaluate at test points, using stored deltas
    D = view(Δ, CartesianIndices(τ), 1:n)
    u = fill(φ[n], size(τ))
    for k in n-1:-1:1
        u = φ[k] .+ D[:, :, k] ./ u
    end
    return u
end

function add_nodes!(r::Thiele, Δ, τ, fτ, new_σ, new_f)
    for (σ, fσ) in zip(new_σ, new_f)
        add_node!(r, σ, fσ)
    end
    σ = r.nodes
    n = length(σ)
    @inbounds @fastmath for i in CartesianIndices(τ)
        Δ[i, n] = τ[i] - σ[n]
    end
    return r
end

function add_node!(r::Thiele, new_σ, new_f)
    σ = r.nodes
    n = length(σ)
    φ = r.weights
    push!(φ, new_f)
    push!(σ, new_σ)
    push!(r.values, new_f)
    n += 1
    for k in 1:n-1
        φ[n] = (σ[n] - σ[k]) / (φ[n] - φ[k])
    end
    return r
end
