
"""
    Barycentric (type)

Barycentric representation of a rational function.

# Fields
- `node`: the nodes of the rational function
- `value`: the values of the rational function
- `weight`: the weights of the rational function
- `wf`: the weighted values of the rational function
- `stats`: convergence statistics
"""
struct Barycentric{T,S} <: AbstractRationalInterpolant{T,S}
    nodes::Vector{S}
    values::Vector{S}
    weights::Vector{S}
    w_times_f::Vector{S}
    stats::Union{Missing,ConvergenceStats{T}}
    function Barycentric{T}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S},
        wf::AbstractVector{S} = value.*weight;
        stats::Union{Missing,ConvergenceStats{T}} = missing
        )  where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight) == length(wf)
        new{T,S}(node, value, weight, wf, stats)
    end
end

"""
    Barycentric(node, value, weight, wf=value .* weight; stats=missing)

Construct a `Barycentric` rational function.

# Arguments
- `node::Vector`: interpolation nodes
- `value::Vector`: values at the interpolation nodes
- `weight::Vector`: barycentric weights
- `wf::Vector`: weights times values (optional)
- `stats::ConvergenceStatistics``: convergence statistics (optional)

# Examples
```jldoctest
julia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])
Barycentric function with 3 nodes and values:
    1.0=>1.0,  2.0=>2.0,  3.0=>3.0

julia> r(1.5)
1.5
```
"""
function Barycentric(node, value, weight, wf=value.*weight; stats=missing)
    Barycentric( promote(float(node), float(value), float(weight))..., float(wf); stats )
end

function Barycentric(
    node::Vector{S}, value::Vector{S}, weight::Vector{S}, wf=value.*weight;
    stats=missing
    ) where {T<:AbstractFloat, S<:RealComplex{T}}
    return Barycentric{T}(node, value, weight, wf; stats)
end

# convenience accessors and overloads
nodes(r::Barycentric) = r.nodes
Base.values(r::Barycentric) = r.values
weights(r::Barycentric) = r.weights
stats(r::Barycentric) = r.stats
"degree(r) returns the degree of the numerator and denominator of the rational `r`."
degree(r::Barycentric) = length(r.nodes) - 1

"""
    r(z)
    evaluate(r, z)

Evaluate the rational function at `z`.
"""

(r::Barycentric)(z) = evaluate(r, z)

function evaluate(r::Barycentric, z::Number)
    if isinf(z)
        return sum(r.w_times_f) / sum(r.weights)
    end
    k = findfirst(z .== r.nodes)
    if isnothing(k)         # not at a node
        C = @. 1 / (z - r.nodes)
        return sum(C .* r.w_times_f) / sum(C .* r.weights)
    else                    # interpolation at node
        return r.values[k]
    end
end

"""
    poles(r)

Return the poles of the rational function `r`.
"""
function poles(r::Barycentric{T}) where T
    w = weights(r)
    nonzero = @. !iszero(w)
    z, w = nodes(r)[nonzero], w[nonzero]
    m = length(w)
    B = diagm( [zero(T); ones(T, m)] )
    E = [zero(T) transpose(w); ones(T, m) diagm(z) ];
    pol = []  # put it into scope
    try
        pol = filter( isfinite, eigvals(E, B) )
    catch
        # generalized eigen not available in extended precision, so:
        λ = filter( z->abs(z)>1e-13, eigvals(E\B) )
        pol = 1 ./ λ
    end
    return pol
end

"""
    residues(r)
    residues(r, p=poles(r))

Return the residues of the rational function `r`. If a vector `p` of poles is given, the
residues are computed at those locations, preserving order.
"""
function residues(r::Barycentric, pol::AbstractVector=poles(r))
    numer = t -> sum( w*y / (t-z) for (z, y, w) in zip(nodes(r), values(r), weights(r)))
    denomdiff = t -> -sum( w / (t-z)^2 for (z, w) in zip(nodes(r), weights(r)))
    res = similar( complex(pol) )
    res .= numer.(pol) ./ denomdiff.(pol)
    return res
end

"""
    roots(r)

Return the roots (zeros) of the rational function `r`.
"""
function roots(r::Barycentric)
    wf = r.w_times_f
    m = length(wf)
    ty = eltype(nodes(r))
    B = diagm( [ty(0); ones(m)] )
    # Thanks to Daan Huybrechs:
    E = [0 transpose(wf); ones(m) diagm(nodes(r))]
    if ty in (Float64,)
        return filter(isfinite, eigvals(E, B))
    else
        # super kludgy; thanks to Daan Huybrechs
        μ = 11 - 17im  # shift since 0 may well be a root
        E -= μ*B
        EB = inv(E)*B
        pp = eigvals(EB)
        large_enough = abs.(pp) .> 1e-10
        return μ .+ 1 ./ pp[large_enough]
    end
end

# Remove any poles on the domain.
function cleanup_poles(f::Approximation, isbad=z->dist(z, f.domain)==0 )
    r = f.fun
    p = filter(!isbad, poles(r))
    res = residues(r, p)
    k = sum(weights(r) .* values(r)) / sum(weights(r))
    pfd(s) = k + sum(a/(s-z) for (a,z) in zip(res, p))
    return pfd
end

"""
    rewind(r, degree)

Rewind a `Barycentric` rational function to a lower degree using stored convergence data.

# Arguments
- `r::Union{Barycentric,Approximation}`: the rational function to rewind
- `degree::Integer`: the degree to rewind to

# Returns
- the rational function of the specified degree (same type as input)

# Examples
```jldoctest
julia> r = aaa(x -> cos(20x), stats=true)
Barycentric function with 25 nodes and values:
    -1.0=>0.408082,  -0.978022=>0.757786,  -0.912088=>0.820908,  …  1.0=>0.408082

julia> rewind(r, 10)
Barycentric function with 11 nodes and values:
    -1.0=>0.408082,  1.0=>0.408082,  -0.466667=>-0.995822,  …  0.898413=>0.636147
```
"""
function rewind(r::Barycentric, degree::Integer)
    m = degree
    return Barycentric(r.stats.nodes[m], r.stats.values[m], r.stats.weights[m])
end

function add_nodes!(r::Barycentric, data, τ, fτ, new_σ, new_f)
    C, L = data
    σ, fσ = r.nodes, r.values
    n, n_new = length(σ), length(new_σ)
    append!(σ, new_σ)
    append!(fσ, new_f)
    append!(r.weights, similar(new_f))
    append!(r.w_times_f, similar(new_f))
    idx_test = CartesianIndices(τ)
    @inbounds @fastmath for i in idx_test, j in n+1:n+n_new
        Δ = τ[i] - σ[j]
        C[i, j] = iszero(Δ) ? 1 / eps() : 1 / Δ
        L[i, j] = (fτ[i] - fσ[j]) * C[i, j]
    end
    # don't recompute the weights yet, since there are still updates needed
    # at new nodes
    return r
end

function update_test_values!(::Type{Barycentric}, numeric_type::Type, num_refine::Integer, max_degree::Integer)
    C = Array{numeric_type}(undef, num_refine, max_degree+1, max_degree+1)
    L = Array{numeric_type}(undef, num_refine, max_degree+1, max_degree+1)
    return C, L
end

function update_test_values!(r::Barycentric, data, τ, fτ, idx_new_test)
    C, L = data
    σ, fσ = r.nodes, r.values
    n = length(σ)
    idx_test = CartesianIndices(τ)

    # update matrices at new test points for all nodes
    @inbounds @fastmath for i in idx_new_test, j in eachindex(σ)
        Δ = τ[i] - σ[j]
        C[i, j] = iszero(Δ) ? 1 / eps() : 1 / Δ
        L[i, j] = (fτ[i] - fσ[j]) * C[i, j]
    end

    # update the weights
    A = reshape(view(L, idx_test, 1:n), :, n)
    _, _, V = svd(A)
    w = V[:, end]
    @. r.weights = w
    @. r.w_times_f = w * fσ

    # evaluate at test points
    # take matrix shape for the active part of C
    Cm = reshape(view(C, idx_test, 1:n), :, n)
    # now use matrix-vector products
    test_values = (Cm*r.w_times_f) ./ (Cm*w)
    return reshape(test_values, size(idx_test))
end
