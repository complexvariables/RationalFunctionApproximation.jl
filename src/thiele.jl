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

function Base.copy(r::Thiele)
    return Thiele(copy(r.nodes), copy(r.values), copy(r.weights))
end

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

function approximate(f::Function, d::Union{ComplexPath,ComplexCurve};
    method::Type{Thiele},
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
    τ = path.points
    idx_test = CartesianIndices((1:1, 2:num_ref+1))
    idx_new_test = idx_test
    fτ = Matrix{eltype(fσ)}(undef, size(τ))        # f at test points
    fτ[idx_test] .= f.(τ[idx_test])
    number_type = Complex{promote_type(eltype(τ), eltype(fτ))}
    values = similar(complex(fτ))

    err = float_type[]    # approximation errors

    # Initialize rational approximation
    r = Thiele(σ, fσ)
    history = IterationRecord{typeof(r),float_type,number_type}[]

    # Main iteration
    n = 1       # iteration counter
    while true
        for i in idx_test
            values[i] = r(τ[i])
        end
        test_values = view(values, idx_test)  # r at test points
        test_actual = view(fτ, idx_test)      # f at test points
        fmax = norm(test_actual, Inf)     # scale of f
        err_max, idx_max = findmax(abs(test_actual[i] - test_values[i]) for i in eachindex(test_actual))
        push!(err, err_max)
        push!(history, IterationRecord(r, err_max, missing))

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping at estimated error $(round(last(err), sigdigits=4)) after $n iterations")
            r = history[status].interpolant
        end
        (status != 0) && break

        ### Refinement
        idx_new = idx_test[idx_max]      # location of worst test point
        add_node!(r, τ[idx_new], fτ[idx_new])
        n += 1

        idx_new_test = add_node!(path, idx_new)
        # In the initial phase, we throw out the old test points.
        if num_ref > refinement    # initial phase
            num_ref -= 3    # gradually decrease refinement level
            s = first(collect(path))
            path = DiscretizedPath(d, s; refinement=num_ref, maxpoints=max_iter * refinement)
            τ = path.points
            idx_test = CartesianIndices((1:n, 2:num_ref+1))
            for i in idx_test
                fτ[i] = f(τ[i])
            end
        else
            # At new test points only, evaluate f.
            for i in idx_new_test
                fτ[i] = f(τ[i])
            end
            idx_test = CartesianIndices((1:n, 2:num_ref+1))
        end

    end
    if allowed === true
        allowed = z -> true
    end
    return Approximation(f, d, r, allowed, path, history)
end


function approximate(y::AbstractVector{T}, z::AbstractVector{S};
    method::Type{Thiele},
    float_type::Type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol::AbstractFloat = 1000*eps(float_type),
    allowed::Union{Function,Bool} = true,
    max_iter::Int = 100,
    stagnation::Int = 16,
    ) where {T<:Number,S<:Number}

    m = length(z)
    n = 1    # iteration counter
    fmax = norm(y, Inf)     # scale of f
    number_type = promote_type(eltype(z), eltype(y))
    err = float_type[]

    _, idx_min = findmin(abs, y)
    values = similar(y)

    r = Thiele([z[idx_min]], [y[idx_min]])
    history = IterationRecord{typeof(r),float_type,number_type}[]
    idx_test = trues(m)
    idx_test[idx_min] = false
    while true
        @. values[idx_test] = r(z[idx_test])  # r at test points
        test_values = view(values, idx_test)  # r at test points
        test_actual = view(y, idx_test)       # f at test points
        err_max, idx_max = findmax(abs(test_actual[i] - test_values[i]) for i in eachindex(test_actual))
        push!(err, err_max)
        push!(history, IterationRecord(r, err_max, missing))

        status = quitting_check(history, stagnation, tol, fmax, max_iter, allowed)
        if status > 0
            @warn("Stopping at estimated error $(round(last(err), sigdigits=4)) after $n iterations")
            r = history[n].interpolant
        end
        (status != 0) && break

        # Add new node:
        idx_new = findall(idx_test)[idx_max]
        add_node!(r, z[idx_new], y[idx_new])
        idx_test[idx_new] = false
        n += 1
    end
    return r, history
end
