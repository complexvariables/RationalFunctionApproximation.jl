#####
##### Approximation on a domain
#####

ComplexRegions.dist(z::Number, p::ComplexRegions.AbstractCurve) = minimum(abs(z - point(p, t)) for t in range(0, 1, length=300))

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
    return Approximation(f, R, r.fun, r.prenodes, r.stats)
end

# Convert curves to paths:
approximate(f::Function, d::ComplexCurve; kw...) = approximate(f, Path(d); kw...)
approximate(f::Function, d::ComplexClosedCurve; kw...) = approximate(f, ClosedPath(d); kw...)

# all methods end up here
function approximate(f::Function, d::ComplexPath;
    method = Barycentric,
    max_iter = 100,
    float_type = promote_type(real_type(d), typeof(float(1))),
    tol = 1000*eps(float_type),
    isbad = z -> dist(z, d) < tol,
    refinement = 3,
    lookahead = 16,
    stats = false
    )

    err = float_type[]
    all_poles = fill(complex(float_type[]), max_iter)
    unacceptable = fill(false, max_iter)

    # Vector of prenodes (except the last, which has no associated test points):
    s = [convert(float_type, 321//654), zero(float_type)]
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
    n = 1    # iteration counter
    numref = 14
    test = Matrix{float_type}(undef, numref, max_iter)    # parameter values of test points
    τ = similar(σ, numref, max_iter + 1)             # test points
    fτ = similar(fσ, numref, max_iter + 1)           # f at test points
    values = similar(fτ)

    # Update the matrices of test points and f values.
    # Each column belongs to one of the node parameter values and contains points
    # that are equispaced up to the next one.
    function update_test_points!(test, τ, fτ, s, Δs, δ, rows, cols)
        @inbounds @fastmath for i in rows, j in cols
            test[i, j] = s[j] + δ[i] * Δs[j]
            τ[i, j] = point(d, test[i, j])
            fτ[i, j] = f(τ[i, j])
        end
        return nothing
    end

    # Initial refinement
    Δs = length(d) * [1 - s[1], s[1]]    # prenode spacings
    δ = float_type.((1:numref) / (numref+1))
    update_test_points!(test, τ, fτ, s, Δs, δ, 1:numref, [1])
    number_type = promote_type(eltype(τ), eltype(fτ))
    data = update_test_values!(method, number_type, numref, max_iter)
    r = method(number_type[], number_type[], number_type[])
    idx_new_test = idx_test = CartesianIndices((1:numref, 1:1))
    @views add_nodes!(r, data, τ[idx_test], fτ[idx_test], σ, fσ)
    lengths = Vector{Int}(undef, max_iter)
    all_weights = Matrix{number_type}(undef, max_iter + 2, max_iter + 2)
    while true
        test_values = update_test_values!(values, r, data, τ, fτ, idx_test, idx_new_test)
        test_actual = view(fτ, idx_test)
        fmax = norm(view(fτ, idx_test), Inf)     # scale of f
        any(isnan.(test_actual)) && throw(ArgumentError("Function has NaN value at a test point"))
        err_max, idx_max = findmax(abs, test_actual - test_values)
        push!(err, err_max)
        lengths[n] = L = length(nodes(r))
        all_weights[1:L, L] .= weights(r)

        # Do we compute poles?
        if stats || (last(err) < tol*fmax)
            all_poles[n] = poles(r)
            n, length(all_poles[n])
            bad = isempty(all_poles[n]) ? BitVector() : isbad.(all_poles[n])
            if count(bad) > 0
                unacceptable[n] = true   # confirmed unacceptable
            end
        end

        # Are we done?
        plateau2 = exp(mean(log.(last(err, lookahead))))
        # plateau2 = 0.9*median(last(err, lookahead))
        plateau = 0.9*mean(last(err, lookahead))
        if (!unacceptable[n] && (last(err) <= tol*fmax)) ||     # goal met
            (n == max_iter) ||          # max iterations reached
            ((n > lookahead) && (plateau2 < last(err) < 1e-2*fmax))    # stagnation
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
        if numref > refinement   # initial phase
            numref -= 1
            δ = float_type.((1:numref) / (numref+1))
            # Wipe out the old test points and start over:
            update_test_points!(test, τ, fτ, s, Δs, δ, 1:numref, 1:n)
            idx_new_test = idx_test = CartesianIndices((1:numref, 1:n))
            r = method(number_type[], number_type[], number_type[])
            @views add_nodes!(r, data, τ[idx_test], fτ[idx_test], σ, fσ)
        else                   # steady-state refinement size
            # Replace one column and add a new column of test points:
            @views add_nodes!(r, data, τ[idx_test], fτ[idx_test], [σ[end]], [fσ[end]])
            update_test_points!(test, τ, fτ, s, Δs, δ, 1:numref, [jnew, n])
            idx_test = CartesianIndices((1:numref, 1:n))
            idx_new_test = idx_test[:, [jnew, n]]
        end
    end

    # among the unchecked approximations, check in order of error
    order = sortperm(err, rev=true)
    idx = 0
    while true
        idx = pop!(order)
        L = lengths[idx]
        if unacceptable[idx]
            # continue, unless out of options
            if isempty(order)
                r = method(σ[1:L], fσ[1:L], all_weights[1:L, L])
                break
            end
        else
            r = method(σ[1:L], fσ[1:L], all_weights[1:L, L])
            if isempty(all_poles[idx])
                all_poles[idx] = poles(r)
            end
            if any(isbad, all_poles[idx])
                unacceptable[idx] = true
            else
                break     # we got it
            end
        end
    end

    # Return the best stuff:
    s = first(s, idx)
    if !isclosed(d)
        push!(s, one(float_type))
    end
    return if stats
        all_poles = all_poles[1:n]
        _weights = [all_weights[1:lengths[i], lengths[i]] for i in 1:n]
        # for i in findall(isempty.(all_poles))
        #     ri = method(_nodes[i], _values[i], _weights[i])
        #     all_poles[i] = poles(ri)
        # end
        _isbad = map(x -> BitVector(map(isbad, x)), all_poles)
        len = lengths[n]
        st = ConvergenceStats(idx, err, σ[1:len], fσ[1:len], _weights, complex.(all_poles), _isbad)
        Approximation(f, d, r, s, st)
    else
        Approximation(f, d, r, s, missing)
    end

end

function clean(r::Approximation;
    tol=1000*eps(real(eltype(nodes(r)))),
    isbad = z -> dist(z, r.domain) < tol,
    )
    t = r.fun
    s = r.prenodes
    remove = [true]
    while any(remove)
        z, y = nodes(t), values(t)
        remove = falses(length(z))
        for p in filter(isbad, poles(t))
            _, idx = findmin(abs, z .- p)
            remove[idx] = true
        end
        remove[3] = false  # hack
        # @show findall(remove)
        @infiltrate
        t = Thiele(z[.!remove], y[.!remove])
        # @show length(s), length(t.nodes)
        # deleteat!(s, findall(remove))
    end
    return Approximation(r.original, r.domain, t, s, missing)
end
