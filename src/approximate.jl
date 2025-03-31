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

# Keywords
- `max_degree::Integer=150`: maximum numerator/denominator degree to use
- `float_type::Type`: floating point type to use for the computation¹
- `tol::Real=1000*eps(float_type)`: relative tolerance for stopping
- `allowed::Function`: function to determine if a pole is allowed
- `refinement::Integer=3`: number of test points between adjacent nodes
- `lookahead::Integer=10`: number of iterations to determine stagnation

¹By default, `float_type` is the promotion of `float(1)` and the float type of the domain.

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
function approximate(f::Function, R::ComplexRegions.AbstractRegion, poles; kw...)
    # Given a region as domain, we interpret poles as not being allowed in that region.
    r = approximate(f, R.boundary, z->!in(z,R), kw...)
    return Approximation(f, R, r.fun, r.allowed, r.prenodes, r.history)
end

# Convert curves to paths, to reduce the number of dispatch points.
function approximate(f::Function, d::ComplexCurve, poles; kw...)
     approximate(f, Path(d), poles; kw...)
end
function approximate(f::Function, d::ComplexClosedCurve, poles; kw...)
    approximate(f, ClosedPath(d), poles; kw...)
end

# Update the matrices of test points and f values.
# Each column belongs to one of the node parameter values and contains points
# that are equispaced up to the next one.
function update_test_points!(path, fun, test, τ, fτ, s, Δs, δ, rows, cols)
    @inbounds @fastmath for i in rows, j in cols
        test[i, j] = s[j] + δ[i] * Δs[j]
        τ[i, j] = point(path, test[i, j])
        fτ[i, j] = fun(τ[i, j])
    end
    return nothing
end

# all methods end up here
function approximate(f::Function, d::ComplexPath, allowed::Union{Function,Bool} = z -> dist(z, d) > tol;
    method = Barycentric,
    max_iter = 100,
    float_type = promote_type(real_type(d), typeof(float(1))),
    tol = 1000*eps(float_type),
    refinement = 3,
    lookahead = 16,
    degree = nothing
    )

    if allowed==true || allowed==:all
        allowed = z -> true
    end
    err = float_type[]
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
    n, n_max = 1, 1    # iteration counter, all-time max
    numref = 14
    test = Matrix{float_type}(undef, numref, max_iter)    # parameter values of test points
    τ = similar(σ, numref, max_iter + 1)             # test points
    fτ = similar(fσ, numref, max_iter + 1)           # f at test points
    values = similar(fτ)

    # Initial refinement
    Δs = length(d) * [1 - s[1], s[1]]    # prenode spacings
    δ = float_type.((1:numref) / (numref+1))
    update_test_points!(d, f, test, τ, fτ, s, Δs, δ, 1:numref, [1])
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
        if any(isnan, test_actual)
            throw(ArgumentError("Function has NaN value at a test point"))
        end
        err_max, idx_max = findmax(abs, test_actual - test_values)
        push!(err, err_max)
        lengths[n] = L = length(nodes(r))
        all_weights[1:L, L] .= weights(r)

        # Are we done?
        if (last(err) <= tol*fmax) && all(allowed, poles(r))    # success
            break
        end

        # Do we quit?
        plateau = exp(median(log.(last(err, lookahead))))
        stagnant = (n > lookahead) && (plateau < last(err) < 1e-2*fmax)
        if (n == max_iter) || stagnant
            @warn("May not have converged to desired tolerance")
            # backtrack to last acceptable approximation
            n += 1
            accepted = false
            while !accepted
                if n == 0
                    @error("No acceptable approximation found")
                end
                L = lengths[n]
                r = method(σ[1:L], fσ[1:L], all_weights[1:L, L])
                accepted = all(allowed, poles(r))
            end
            break
        end

        # Add new node:
        jnew = idx_max[2]           # which node owns the worst test point?
        push!(s, test[idx_max])
        push!(σ, τ[idx_max])
        push!(fσ, fτ[idx_max])
        n_max = n += 1

        # Update the node spacings:
        Δj = Δs[jnew]
        Δs[jnew] = test[idx_max] - s[jnew]
        push!(Δs, s[jnew] + Δj - test[idx_max])

        # Update test points and matrices:
        if numref > refinement   # initial phase
            numref -= 1
            δ = float_type.((1:numref) / (numref+1))
            # Wipe out the old test points and start over:
            update_test_points!(d, f, test, τ, fτ, s, Δs, δ, 1:numref, 1:n)
            idx_new_test = idx_test = CartesianIndices((1:numref, 1:n))
            r = method(number_type[], number_type[], number_type[])
            @views add_nodes!(r, data, τ[idx_test], fτ[idx_test], σ, fσ)
        else                   # steady-state refinement size
            # Replace one column and add a new column of test points:
            @views add_nodes!(r, data, τ[idx_test], fτ[idx_test], [σ[end]], [fσ[end]])
            update_test_points!(d, f, test, τ, fτ, s, Δs, δ, 1:numref, [jnew, n])
            idx_test = CartesianIndices((1:numref, 1:n))
            idx_new_test = idx_test[:, [jnew, n]]
        end
    end

    # Return the best stuff:
    s = first(s, lengths[n])
    if !isclosed(d)
        push!(s, one(float_type))
    end
    hist = History{float_type}(σ, fσ, all_weights, lengths[1:n_max], n)
    return Approximation(f, d, r, allowed, s, hist)
end

# function clean(r::Approximation;
#     tol=1000*eps(real(eltype(nodes(r)))),
#     isbad = z -> dist(z, r.domain) < tol,
#     )
#     t = r.fun
#     s = r.prenodes
#     remove = [true]
#     while any(remove)
#         z, y = nodes(t), values(t)
#         remove = falses(length(z))
#         for p in filter(isbad, poles(t))
#             _, idx = findmin(abs, z .- p)
#             remove[idx] = true
#         end
#         remove[3] = false  # hack
#         # @show findall(remove)
#         @infiltrate
#         t = Thiele(z[.!remove], y[.!remove])
#         # @show length(s), length(t.nodes)
#         # deleteat!(s, findall(remove))
#     end
#     return Approximation(r.original, r.domain, t, s, missing)
# end

function get_history(r::Approximation{T,S}) where {T,S}
    hist = r.history
    deg = Int[]
    zp = Vector{complex(S)}[]
    err = T[]
    allowed = BitVector[]
    τ, _ = check(r, quiet=true)
    fτ = r.original.(τ)
    scale = maximum(abs, fτ)
    for (idx, n) in enumerate(hist.len)
        rn = rewind(r, idx)
        push!(deg, degree(rn))
        push!(zp, poles(rn))
        push!(allowed, r.allowed.(zp[end]))
        push!(err, maximum(abs, rn.(τ) - fτ) / scale)
    end
    best = findfirst(length(nodes(r)) .== hist.len)
    return deg, err, zp, allowed, best
end

# Discretize a path so that the distance to neighbors is no more than half the distance to the nearest pole.
function _discretize(d::ComplexPath, poles::AbstractVector)
    function pole_to_nbr_ratio(z, zp, k)
        to_pole = sqrt(minimum(abs2(z[k] - p) for p in zp; init=Inf))
        to_nbr = max(abs(z[k] - z[k-1]), abs(z[k] - z[k+1]))
        return to_pole / to_nbr
    end

    t = range(0, length(d), 100)
    z = point(d, t)
    while any(pole_to_nbr_ratio(z, poles, k) < 2 for k in 2:length(z)-1)
        t, z = refine(d, t, 2)
        if length(t) > 70_000
            @error("Unable to refine the path sufficiently")
            break
        end
    end
    return t, z
end

function approximate(f::Function, d::ComplexPath, allowed::AbstractVector;
    degree = max(4, div(length(allowed), 2))
    )

    t, z = _discretize(d, allowed)
    B = ArnoldiBasis(z, degree)
    C = [1 / (z - zp) for z in z, zp in allowed]
    c = isempty(C) ? B.Q \ f.(z) : [B.Q C] \ f.(z)
    p = ArnoldiPolynomial(c[1:degree+1], B)
    return PartialFractions(p, allowed, c[degree+2:end])
end
