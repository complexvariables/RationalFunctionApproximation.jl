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

# First call to allocate space for the pairwise difference matrix
function update_test_values!(
    ::Type{Thiele},
    numeric_type::Type,
    num_refine::Integer,
    max_iter::Integer,
    max_test::Integer=max_iter+2
    )
    Δ = Array{numeric_type}(undef, max_test, num_refine+1, max_iter + 2)
    return Δ
end

# Update the difference matrix at new test points, then evaluate at all the test points.
function update_test_values!(values, r::Thiele, Δ, τ, fτ, idx_test, idx_new_test)
    σ, φ = r.nodes, r.weights
    n = length(σ)

    # Update the difference matrix at new test points for all nodes
    @inbounds @fastmath for i in idx_new_test, j in eachindex(σ)
        Δ[i, j] = τ[i] - σ[j]
    end

    # Evaluate at all test points
    V = view(values, idx_test)
    V .= φ[n]
    @inbounds @fastmath for k in n-1:-1:1, i in eachindex(idx_test)
        V[i] = φ[k] + Δ[idx_test[i], k] / V[i]
    end
    return V
end

# add multiple interpolation nodes, updating the matrix of node-test distances
function add_nodes!(r::Thiele, Δ, τ, fτ, idx_test, new_σ, new_f)
    for (σ, fσ) in zip(new_σ, new_f)
        add_node!(r, σ, fσ)
    end
    σ = r.nodes
    n = length(σ)
    @inbounds @fastmath for i in idx_test
        Δ[i, n] = τ[i] - σ[n]
    end
    return r
end

# update in-place the nodes and weights
function add_node!(r::Thiele, new_σ, new_f)
    σ = r.nodes
    n = length(σ)
    φ = r.weights
    push!(φ, new_f)
    push!(σ, new_σ)
    push!(r.values, new_f)
    n += 1
    for k in 1:n-1
        d = φ[n] - φ[k]
        if iszero(d)
            error("Infinite inverse difference. Try breaking symmetry of the function.")
        end
        @inbounds @fastmath φ[n] = (σ[n] - σ[k]) / d
    end
    return r
end
