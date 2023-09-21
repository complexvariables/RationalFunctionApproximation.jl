#####
##### Discrete AAA
#####

"""
    aaa(z, y)
    aaa(f)

Adaptively compute a rational interpolant.

# Arguments

## discrete mode
- `z::AbstractVector{<:Number}`: interpolation nodes
- `y::AbstractVector{<:Number}`: values at nodes

## continuous mode
- `f::Function`: function to approximate on the interval [-1,1]

# Keyword arguments
- `max_degree::Integer=150`: maximum numerator/denominator degree to use
- `float_type::Type=Float64`: floating point type to use for the computation
- `tol::Real=1000*eps(float_type)`: tolerance for stopping
- `lookahead::Integer=10`: number of iterations to determines stagnation
- `stats::Bool=false`: return convergence statistics

# Returns
- `r::Barycentric`: the rational interpolant
- `stats::NamedTuple`: convergence statistics, if keyword `stats=true`

See also [`approximate`](@ref) for approximating a function on a region.
"""
function aaa(z::AbstractVector{<:Number}, y::AbstractVector{<:Number};
    max_degree = 150, float_type = Float64, tol = 1000*eps(float_type),
    lookahead = 10, stats = false
    )

    @assert float_type <: AbstractFloat
    T = float_type
    fmax = norm(y, Inf)    # for scaling
    m = length(z)
    iteration = NamedTuple[]
    err = T[]
    besterr, bestidx, best = Inf, NaN, nothing

    # Allocate space for Cauchy matrix, Loewner matrix, and residual
    C = similar(z, (m, m))
    L = similar(z, (m, m))
    R = complex(zeros(size(z)))

    ȳ = sum(y) / m
    s, idx = findmax(abs(y - ȳ) for y in y)
    push!(err, s)

    # The ordering of nodes matters, while the order of test points does not.
    node_index = Int[]
    push!(node_index, idx)
    test_index = Set(1:m)
    delete!(test_index, idx)

    n = 0    # number of poles
    while true
        n += 1
        σ = view(z, node_index)
        fσ = view(y, node_index)
        # Fill in matrices for the latest node
        @inbounds @fastmath for i in test_index
            C[i, n] = 1 / (z[i] - σ[n])
            L[i, n] = (y[i] - fσ[n]) * C[i, n]
        end

        istest = collect(test_index)
        _, _, V = svd( view(L, istest, 1:n) )
        w = V[:, end]    # barycentric weights

        CC = view(C, istest, 1:n)
        num = CC * (w.*fσ)
        den = CC * w
        @. R[istest] = y[istest] - num / den
        push!(err, norm(R, Inf))
        push!(iteration, (; weights=w, active=copy(node_index)))

        if (last(err) < besterr)
            besterr, bestidx, best = last(err), length(iteration), last(iteration)
        end

        # Are we done?
        if (besterr <= tol*fmax) ||
            (n == max_degree + 1) ||
            ((length(iteration) - bestidx >= lookahead) && (besterr < 1e-2*fmax))
            break
        end

        _, j = findmax(abs, R)
        push!(node_index, j)
        delete!(test_index, j)
        R[j] = 0
    end

    idx, w = best.active, best.weights
    r = Barycentric(z[idx], y[idx], w)
    if stats
        return r, (;err, iteration)
    else
        return r
    end
end

#####
##### Adaptive AAA on [-1, 1] only
#####
# refinement in parameter space
function refine(t, N)
    x = sort(t)
    Δx = diff(x)
    d = eltype(x).((1:N) / (N+1))
    return vec( x[1:end-1] .+ (d' .* Δx) )
end

function aaa(
    f::Function;
    max_degree=150, float_type=Float64, tol=1000*eps(float_type),
    refinement=3, lookahead=10, stats=false
    )
    @assert float_type <: AbstractFloat
    T = float_type
    CT = Complex{T}
    # arrays for tracking convergence progress
    err, nbad = T[], Int[]
    nodes, vals, pol, weights = Vector{T}[], Vector{CT}[], Vector{CT}[], Vector{CT}[]

    S = [-one(T), one(T)]                       # initial nodes
    fS = f.(S)
    besterr, bestm = Inf, NaN
    while true                                  # main loop
        m = length(S)
        push!(nodes, copy(S))
        X = refine(S, max(refinement, ceil(16-m)))    # test points
        fX = f.(X)
        push!(vals, copy(fS))
        C = [ 1/(x-s) for x in X, s in S ]
        L = [a-b for a in fX, b in fS] .* C
        _, _, V = svd(L)
        w = V[:,end]
        push!(weights, w)
        R = (C*(w.*fS)) ./ (C*w)                # values of the rational interpolant
        push!(err, norm(fX - R, Inf) )

        zp =  poles(Barycentric(S, fS, w))
        push!(pol, zp)
        I = (imag(zp).==0) .& (abs.(zp).<=1)    # bad poles indicator
        push!(nbad, sum(I))
        # If valid and the best yet, save it:
        if (last(nbad) == 0) && (last(err) < besterr)
            besterr, bestm = last(err), m
        end

        fmax = max( norm(fS, Inf), norm(fX, Inf) )     # scale of f
        # Check stopping:
        if (besterr <= tol*fmax) ||                                # goal met
            (m == max_degree + 1) ||                                   # max degree reached
            ((m - bestm >= lookahead) && (besterr < 1e-2*fmax))    # stagnation
            break
        end

        # We're continuing the iteration, so add the worst test point to the nodes:
        _, j = findmax(abs, fX - R)
        push!(S, X[j])
        push!(fS, fX[j])
    end

    # Use the best result found:
    S, y, w = nodes[bestm-1], vals[bestm-1], weights[bestm-1]
    idx = sortperm(S)
    x, y, w = S[idx], y[idx], w[idx]
    if isreal(w) && isreal(y)
        y, w = real(y), real(w)
    end

    if stats
        if isreal(w) && isreal(y)
            weights = real.(weights)
            vals = real.(vals)
        end
        st = ConvergenceStats(bestm-1, err, nbad, nodes, vals, weights, pol)
        r = Barycentric(x, y, w; stats=st)
    else
        r = Barycentric(x, y, w)
    end

    return r
end
