#####
##### Approximation on a domain
#####

# refinement on a curve
refine(p::ComplexCurve, args...) = refine(Path(p), args...)
function refine(p::ComplexPath, t::AbstractVector, N::Integer=3)
    tt = refine(t, N)
    ττ = point(p, tt)
    if isreal(p)
        ττ = real(ττ)
    end
    # TODO: filter out Inf/Nan?
    return tt, ττ
end

"""
    approximate(f, domain)

Adaptively compute a rational interpolant on a curve, path, or region.

# Arguments
- `f::Function`: function to approximate
- `domain`: curve, path, or region from ComplexRegions

# Keyword arguments
- `max_degree::Integer=150`: maximum numerator/denominator degree to use
- `float_type::Type=Float64`: floating point type to use for the computation
- `tol::Real=1000*eps(float_type)`: relative tolerance for stopping
- `isbad::Function`: function to determine if a pole is bad
- `refinement::Integer=3`: number of test points between adjacent nodes
- `lookahead::Integer=10`: number of iterations to determine stagnation
- `stats::Bool=false`: whether to return convergence statistics with the approximation (slower)

# Returns
- `r::Approximation`: the rational interpolant

See also [`Approximation`](@ref), [`check`](@ref), [`aaa`](@ref).

# Examples
```julia-repl
julia> f(x) = tanh( 40*(x - 0.15) );

julia> r = approximate(f, unit_interval)
Barycentric rational function of type (22,22) on the domain: Path with 1 curve

julia> ( r(0.3), f(0.3) )
(0.9999877116507944, 0.9999877116507956)

julia> check(r);   # accuracy over the domain
[ Info: Max error is 7.09e-14
```
"""
function approximate(f::Function, R::ComplexRegions.AbstractRegion; kw...)
    r = approximate(f, R.boundary; isbad=z->in(z,R), kw...)
    return Approximation(f, R, r.fun, r.prenodes)
end

# Convert curves to paths:
approximate(f::Function, d::ComplexCurve; kw...) = approximate(f, Path(d); kw...)
approximate(f::Function, d::ComplexClosedCurve; kw...) = approximate(f, ClosedPath(d); kw...)

function approximate(f::Function, d::ComplexPath;
    max_degree = 150,
    float_type = typeof(float(1)),
    tol = 1000*eps(float_type),
    isbad = z->dist(z, d) < tol,
    refinement = 3,
    lookahead = 10,
    stats = false
    )

    @assert float_type <: AbstractFloat
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
    C = similar(τ, 17*(max_degree+1), max_degree+2)    # Cauchy matrix
    ty = promote_type(eltype(τ), eltype(fτ))
    L = Matrix{ty}(undef, 17*(max_degree+1), max_degree+2)    # Lowener matrix

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
    function update_matrices!(C, L, τ, fτ, σ, fσ, rows, cols)
        @inbounds @fastmath for i in rows, j in cols
            Δ = τ[i] - σ[j]
            C[i, j] = iszero(δ) ? 1 / eps() : 1 / Δ
            L[i, j] = (fτ[i] - fσ[j]) * C[i, j]
        end
        return nothing
    end

    # Initial refinement
    Δs = [one(float_type)]    # prenode spacings
    δ = float_type.((1:numref) / (numref+1))
    update_test!(test, τ, fτ, s, Δs, δ, 1:numref, [1])
    τm = vec(view(τ, 1:numref, 1:m))
    fτm = vec(view(fτ, 1:numref, 1:m))
    update_matrices!(C, L, τm, fτm, σ, fσ, 1:numref, 1:n)

    while true
        # Barycentric weights:
        numtest = length(τm)
        _, _, V = svd(view(L, 1:numtest, 1:n))
        w = V[:, end]

        # Residual at test points:
        Cm = view(C, 1:numtest, 1:n)                  # active part of the Cauchy matrix
        R = reshape((Cm*(w.*fσ)) ./ (Cm*w), numref, :)
        err_max, idx_max = findmax(abs, view(fτ, 1:numref, 1:m) - R)
        push!(err, err_max)

        # Poles:
        zp =  poles(Barycentric(σ, fσ, w))
        I = isbad.(zp)
        push!(iteration, (; weights=w, vals=copy(fσ), poles=zp, prenodes=copy(s), nodes=copy(σ)))
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
            τm = vec(view(τ, 1:numref, 1:m))
            fτm = vec(view(fτ, 1:numref, 1:m))
            update_matrices!(C, L, τm, fτm, σ, fσ, eachindex(τm), eachindex(σ))
        else                   # steady-state refinement size
            # Replace one column and add a new column of test points:
            update_test!(test, τ, fτ, s, Δs, δ, 1:numref, [jnew, m])
            τm = view(τ, 1:numref, 1:m)
            fτm = view(fτ, 1:numref, 1:m)
            idx = LinearIndices(τm)
            # New test points at the old nodes:
            rows = [idx[:, jnew]; idx[:, m]]
            update_matrices!(C, L, τm, fτm, σ, fσ, rows, 1:n-1)
            # New node at all the test points:
            update_matrices!(C, L, τm, fτm, σ, fσ, vec(idx), [n])
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
        r = Barycentric(z, y, w; stats=st)
    else
        r = Barycentric(z, y, w)
    end
    return Approximation(f, d, r, s)
end

"""
    check(r)

Check the accuracy of a rational approximation `r` on its domain.

# Arguments
- `r::Approximation`: rational approximation

# Returns
- `τ::Vector`: test points
- `err::Vector`: error at test points

See also [`approximate`](@ref), [`aaa`](@ref).
"""
function check(F::Approximation)
    p = F.domain
    if p isa ComplexSCRegion
        p = p.boundary
    end
    t, τ = refine(p, F.prenodes, 30)
    if isreal(F.fun.nodes)
        τ = real(τ)
    end
    err = F.original.(τ) - F.fun.(τ)
    @info f"Max error is {norm(err,Inf):.2e}"
    return τ, err
end
