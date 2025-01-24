struct Thiele{T,S} <: AbstractRationalInterpolant{T,S}
    nodes::Vector{S}
    values::Vector{S}
    weights::Vector{S}
    stats::Union{Missing,ConvergenceStats{T}}
    function Thiele{T}(
        node::AbstractVector{S},
        value::AbstractVector{S},
        weight::AbstractVector{S};
        stats::Union{Missing,ConvergenceStats{T}} = missing
        )  where {T <: AbstractFloat, S <: RealComplex{T}}
        @assert length(node) == length(value) == length(weight)
        new{T,S}(node, value, weight, stats)
    end
end

nodes(r::Thiele) = r.nodes
Base.values(r::Thiele) = r.values
weights(r::Thiele) = r.weights
stats(r::Thiele) = r.stats

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
    C = diagm(1 => -ones(T, n-2))
    xi = view(r.nodes, 2:n-1)
    D = diagm(0 => view(r.weights, 2:n), -1 => -ones(T, n-2), 1 => -xi)
    z = try
        filter(isfinite, eigvals(D, C))
    catch
        @warn "Problem with generalized eigvals"
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

function Thiele(nodes::AbstractVector, values::AbstractVector, weights::AbstractVector; stats=missing)
    if isempty(nodes) && eltype(nodes) == Any
        nodes = values = weights = Float64[]
    end
    T = eltype(values)
    return Thiele{real_type(T)}(nodes, values, weights; stats)
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
    Δ = Array{numeric_type}(undef, num_refine, max_degree+1, max_degree+1)
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


function t_approximate(f::Function, d::ComplexPath;
    max_degree = 150,
    float_type = promote_type(real_type(d), typeof(float(1))),
    tol = 1000*eps(float_type),
    isbad = z->dist(z, d) < tol,
    refinement = 3,
    lookahead = 10,
    stats = false
    )

    err, nbad = float_type[], Int[]
    iteration = NamedTuple[]

    # Vector of prenodes (except the last, which has no associated test points):
    s = [zero(float_type)]
    σ = point(d, s)            # node points
    if !isclosed(d)
        # Include both endpoints as nodes for open paths:
        push!(σ, point(d, length(d)*one(float_type)))
    end
    if isreal(d)
        # Enforce reality:
        σ = real(σ)
    end
    fσ = f.(σ)           # f at nodes
    n = length(σ)        # number of nodes
    m = length(s)        # number of prenodes having test points to the right

    besterr, bestidx, best = Inf, NaN, nothing
    numref = 16
    test = Matrix{float_type}(undef, numref, max_degree+1)    # parameter values of test points
    τ = similar(σ, numref, max_degree+1)             # test points
    fτ = similar(fσ, numref, max_degree+1)           # f at test points
    Δ = similar(τ, numref, max_degree+1, max_degree+2)    # Differences matrix

    # Update the matrices of test points and f values.
    # Each column belongs to one of the node parameter values and contains points
    # that are equispaced up to the next one.
    function update_test!(test, τ, fτ, s, Δs, δ, rows, cols)
        @inbounds @fastmath for i in rows, j in cols
            test[i, j] = s[j] + δ[i] * Δs[j]
            τ[i, j] = point(d, test[i, j])
            fτ[i, j] = f(τ[i, j])
        end
        return nothing
    end

    # Update the Cauchy and Loewner matrices.
    function update_matrices!(Δ, τ, σ, idx_test, idx_node)
        @inbounds @fastmath for i in idx_test, j in idx_node
            Δ[i, j] = τ[i] - σ[j]
        end
        return nothing
    end

    # Initial refinement
    Δs = [length(d) * one(float_type)]    # prenode spacings
    δ = float_type.((1:numref) / (numref+1))
    update_test!(test, τ, fτ, s, Δs, δ, 1:numref, [1])
    τm = view(τ, 1:numref, 1:m)
    fτm = view(fτ, 1:numref, 1:m)
    update_matrices!(Δ, τm, σ, CartesianIndices(τm), 1:n)
    update_approximant!(f, d, σ, fσ, τ, fτ, Δ, test, 1:numref, 1:m, max_degree, float_type, tol, lookahead, stats)

    φ = copy(fσ)
    for i in eachindex(σ)
        for k in 1:i-1
            φ[i] = (σ[i] - σ[k]) / (φ[i] - φ[k])
        end
    end
    while true
        # Residual at test points:
        Dm = view(Δ, CartesianIndices(τm), 1:m)                  # active part of the delta matrix
        z = φ[n]
        for k in n-1:-1:1
            z = φ[k] .+ Dm[:, :, k] ./ z
        end

        err_max, idx_max = findmax(abs, fτm - z)
        push!(err, err_max)
        any(isnan.(view(fτ, 1:numref, 1:m))) && throw(ArgumentError("Function has NaN value at a test point"))

        # # Poles:
        # zp =  poles(Thiele(σ, fσ, w))
        # I = isbad.(zp)
        zp = []
        push!(iteration, (; weights=copy(φ), vals=copy(fσ), poles=zp, prenodes=copy(s), nodes=copy(σ)))
        push!(nbad, isempty(zp) ? 0 : sum(I))

        # Record the best approximation so far:
        numiter = length(iteration)
        if (last(nbad) == 0) && (last(err) < besterr)
            besterr, bestidx, best = last(err), numiter, last(iteration)
        end

        # Are we done?
        fmax = norm(fτm, Inf)     # scale of f
        if (besterr <= tol*fmax) ||     # goal met
            (n-1 == max_degree) ||          # max degree reached
            ((numiter - bestidx >= lookahead) && (besterr < 1e-2*fmax))    # stagnation
            break
        end

        # Add new node:
        jnew = idx_max[2]           # which node owns the worst test point?
        push!(s, test[idx_max])
        push!(σ, τ[idx_max])
        push!(fσ, fτ[idx_max])
        n += 1

        push!(φ, fτ[idx_max])
        for k in 1:n-1
            φ[n] = (σ[n] - σ[k]) / (φ[n] - φ[k])
        end

        # Update the node spacings:
        Δj = Δs[jnew]
        Δs[jnew] = test[idx_max] - s[jnew]
        push!(Δs, s[jnew] + Δj - test[idx_max])

        # Update test points and matrices:
        m += 1
        if numref > refinement   # initial phase
            numref -= 1
            δ = float_type.((1:numref) / (numref+1))
            # Wipe out the old test points and start over:
            update_test!(test, τ, fτ, s, Δs, δ, 1:numref, 1:m)
            τm = view(τ, 1:numref, 1:m)
            fτm = view(fτ, 1:numref, 1:m)
            update_matrices!(Δ, τm, σ, eachindex(τm), eachindex(σ))
        else                   # steady-state refinement size
            # Replace one column and add a new column of test points:
            update_test!(test, τ, fτ, s, Δs, δ, 1:numref, [jnew, m])
            τm = view(τ, 1:numref, 1:m)
            fτm = view(fτ, 1:numref, 1:m)
            idx = CartesianIndices(τm)
            # New test points at the old nodes:
            update_matrices!(Δ, τm, σ, idx[:, [jnew, m]], 1:n-1)
            # New node at all the test points:
            update_matrices!(Δ, τm, σ, idx, [n])
        end
    end

    # Return the best stuff:
    z, y, w, s = best.nodes, best.vals, best.weights, best.prenodes
    if !isclosed(d)
        push!(s, one(float_type))
    end
    if isreal(z) && isreal(w) && isreal(y)
        z, y, w = real(z), real(y), real(w)
    end
    if stats
        nodes = collect(i.nodes for i in iteration)
        vals = collect(i.vals for i in iteration)
        if isreal(z) && isreal(w) && isreal(y)
            vals = real.(vals)
        end
        weights = collect(i.weights for i in iteration)
        pole = collect(i.poles for i in iteration)
        st = ConvergenceStats(bestidx, err, nbad, nodes, vals, weights, complex.(pole))
        r = Thiele(z, y, w; stats=st)
    else
        r = Thiele(z, y, w)
    end
    # return Approximation(f, d, r, s)
    return r, s
end
