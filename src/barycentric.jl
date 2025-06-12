
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
    function Barycentric{T}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S},
        wf::AbstractVector{S} = value.*weight
        )  where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight) == length(wf)
        new{T,S}(node, value, weight, wf)
    end
    function Barycentric{T,S}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S},
        wf::AbstractVector{S} = value.*weight
        )  where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight) == length(wf)
        new{T,S}(node, value, weight, wf)
    end
end

"""
    Barycentric(node, value, weight, wf=value .* weight; stats=missing)

Construct a `Barycentric` rational function.

# Arguments
- `node::Vector`: interpolation nodes
- `value::Vector`: values at the interpolation nodes
- `weight::Vector`: barycentric weights

# Keywords
- `wf::Vector = value .* weight`: weights times values

# Returns
- `::Barycentric`: a barycentric rational interpolating function

# Examples
```jldoctest
julia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])
Barycentric function with 3 nodes and values:
    1.0=>1.0,  2.0=>2.0,  3.0=>3.0

julia> r(1.5)
1.5
```
"""
function Barycentric(node, value, weight, wf=value.*weight)
    Barycentric( promote(float(node), float(value), float(weight))..., float(wf))
end

function Barycentric(
    node::Vector{S}, value::Vector{S}, weight::Vector{S}, wf=value.*weight
    ) where {T<:AbstractFloat, S<:RealComplex{T}}
    return Barycentric{T}(node, value, weight, wf)
end

# construct from Loewner matrix
function Barycentric(node::AbstractVector, value::AbstractVector, L::AbstractMatrix)
    _, _, V = svd(L)
    return Barycentric(node, value, V[:, end])
end

Base.copy(r::Barycentric) =
    Barycentric(copy(r.nodes), copy(r.values), copy(r.weights), copy(r.w_times_f))

# convenience accessors and overloads
nodes(r::Barycentric) = r.nodes
Base.values(r::Barycentric) = r.values
weights(r::Barycentric) = r.weights
degrees(r::Barycentric) = (length(r.nodes) - 1, length(r.nodes) - 1)
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

# Evaluate when given Cauchy matrix
function evaluate!(u::AbstractArray, r::Barycentric, C::AbstractMatrix)
    n = length(r.nodes)
    for (i, idx) in enumerate(eachindex(u))
        u[idx] = sum(C[i, j] * r.w_times_f[j] for j in 1:n) / sum(C[i, j] * r.weights[j] for j in 1:n)
    end
    return nothing
end

"""
    poles(r)

Return the poles of the rational function `r`.
"""
function poles(r::Barycentric{T,S}) where {T,S}
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

function residues(r::Barycentric)
    ζ = poles(r)
    res = similar( complex(ζ) )
    z, y, w = nodes(r), values(r), weights(r)
    for (i, t) in pairs(ζ)
        numer = sum( w*y / (t-z) for (z, y, w) in zip(z, y, w))
        denomdiff = -sum( w / (t-z)^2 for (z, w) in zip(z, w))
        res[i] = numer / denomdiff
    end
    return ζ, res
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
        # use generic linear algebra
        _, _, _, _, ⍺, β = schur(complex(E), complex(B))
        filter(isfinite, ⍺ ./ β)
    end
end

# add new nodes to an existing Barycentric function
function add_node(r::Barycentric, C, L, new_σ, new_fσ, τ, fτ, idx_test, idx_new_test)
    σ =  [r.nodes;   new_σ]
    fσ = [r.values; new_fσ]
    n = length(σ)

    # update matrices for all the test points for the new node
    @inbounds @fastmath for i in idx_test
        Δ = τ[i] - σ[n]
        C[i, n] = iszero(Δ) ? 1 / eps() : 1 / Δ
        L[i, n] = (fτ[i] - fσ[n]) * C[i, n]
    end

    # update matrices for new test points for all nodes
    @inbounds @fastmath for i in idx_new_test, j in eachindex(σ)
        Δ = τ[i] - σ[j]
        C[i, j] = iszero(Δ) ? 1 / eps() : 1 / Δ
        L[i, j] = (fτ[i] - fσ[j]) * C[i, j]
    end

    # collect the columns of the Loewner matrix
    if isa(idx_test, CartesianIndices)
        A = reshape(view(L, idx_test, 1:n), :, n)
    else
        # slower route when test point indices are not contiguous
        I = [CartesianIndex((idx, k)) for idx in idx_test, k in 1:n]
        A = reshape(view(L, I), :, n)
    end

    return Barycentric(σ, fσ, A)
end

function _initialize(f, num_ref, max_iter, σ, fσ, path, idx_test)
    τ = path.points
    fτ = Matrix{eltype(fσ)}(undef, size(τ))        # f at test points
    @. fτ[idx_test] = f(τ[idx_test])
    rτ = complex(similar(fτ))    # rational function at test points
    value_type = promote_type(eltype(σ), eltype(fτ))
    C = Array{value_type}(undef, max_iter + 1, num_ref + 1, max_iter + 1)
    L = Array{value_type}(undef, max_iter + 1, num_ref + 1, max_iter + 1)
    @inbounds @fastmath for i in idx_test
        Δ = τ[i] - σ[1]
        C[i, 1] = iszero(Δ) ? 1 / eps() : 1 / Δ
        L[i, 1] = (fτ[i] - fσ[1]) * C[i, 1]
    end
    return τ, fτ, rτ, C, L
end

function approximate(method::Type{Barycentric},
    f::Function, d::ComplexCurveOrPath;
    float_type::Type = promote_type(real_type(d), typeof(float(1))),
    tol::Real = 1000*eps(float_type),
    allowed::Union{Function,Bool} = z -> dist(z, d) > tol,
    max_iter::Int = 150,
    refinement::Int = 3,
    stagnation::Int = 10
    )

    num_ref = 15    # initial number of test points between nodes; decreases to `refinement`
    path = DiscretizedPath(d, [0, 1]; refinement=num_ref, maxpoints=max_iter * refinement)
    σ = [point(d, 0)]
    fσ = f.(σ)         # f at nodes

    # Arrays of test points have one row per node (except the last)
    idx_test = CartesianIndices((1:1, 2:num_ref+1))
    τ, fτ, rτ, C, L = _initialize(f, num_ref, max_iter, σ, fσ, path, idx_test)
    fmax = maximum(abs, view(fτ, idx_test))        # scale of f

    # Initialize rational approximation
    r = Barycentric(σ, fσ, reshape(view(L, idx_test, 1:1), :, 1))
    history = [IterationRecord(r, NaN, missing)]

    # Main iteration
    n = 1       # iteration counter
    while true
        Cmatrix = reshape(view(C, idx_test, 1:n), :, n)
        evaluate!(view(rτ, idx_test), r, Cmatrix)    # r at test points
        err = @. abs(fτ[idx_test] - rτ[idx_test])
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
        new_σ, new_fσ = τ[idx_new], fτ[idx_new]
        # In the initial phase, we throw out the old test points.
        idx_new_test = add_node!(path, idx_new)
        if num_ref > refinement    # initial phase
            num_ref -= 3    # gradually decrease refinement level
            s, _ = collect(path, :nodes)
            path = DiscretizedPath(d, s; refinement=num_ref, maxpoints=max_iter * refinement)
            idx_new_test = idx_test = CartesianIndices((1:n+1, 2:num_ref+1))
            τ, fτ, rτ, C, L = _initialize(f, num_ref, max_iter, σ, fσ, path, idx_test)
            fmax = maximum(abs, view(fτ, idx_test))        # scale of f
        else
            idx_test = CartesianIndices((1:n+1, 2:num_ref+1))
            @. fτ[idx_new_test] = f(τ[idx_new_test])
        end
        r = add_node(r, C, L, new_σ, new_fσ, τ, fτ, idx_test, idx_new_test)
        push!(history, IterationRecord(r, NaN, missing))
        n += 1
    end
    if allowed === true
        allowed = z -> true
    end
    return Approximation(f, d, r, allowed, path, history)
end

function approximate(method::Type{Barycentric},
    y::AbstractVector{T}, z::AbstractVector{S};
    float_type::Type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::AbstractFloat = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = 100,
    stagnation::Int = 10,
    ) where {T<:Number,S<:Number}

    m = length(z)
    fmax = maximum(abs, y)     # scale of f
    values = similar(y)

    _, i₀ = findmin(abs, y)
    idx_test = trues(m)
    idx_test[i₀] = false

    number_type = promote_type(eltype(z), eltype(y))
    C = Array{number_type}(undef, m, max_iter+1)
    L = Array{number_type}(undef, m, max_iter+1)
    @inbounds for i in 1:m
        C[i, 1] = 1 / (z[i] - z[i₀])
        L[i, 1] = (y[i] - y[i₀]) * C[i, 1]
    end
    r = Barycentric([z[i₀]], [y[i₀]], view(L, idx_test, 1:1))
    history = [IterationRecord(r, NaN, missing)]
    n = 1    # iteration counter
    while count(idx_test) > 0
        evaluate!(view(values, idx_test), r, view(C, idx_test, 1:n))    # r at test points
        err = @. abs(y[idx_test] - values[idx_test])
        err_max, idx_max = findmax(err)
        history[n].error = err_max

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping at estimated error $(round(err_max, sigdigits=4)) after $n iterations")
            r = history[n].interpolant
        end
        (status != 0) && break

        # Add new node:
        idx_new = findall(idx_test)[idx_max]
        idx_test[idx_new] = false
        r = add_node(r, C, L, z[idx_new], y[idx_new], z, y, findall(idx_test), [])
        push!(history, IterationRecord(r, NaN, missing))
        n += 1
    end
    return r, history
end
