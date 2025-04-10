#####
##### TYPES
#####

##### Approximation

"""
    Approximation (type)

Approximation of a function on a domain.

# Fields
- `original`: the original function
- `domain`: the domain of the approximation
- `fun`: the barycentric representation of the approximation
- `prenodes`: the prenodes of the approximation
- `stats`: convergence statistics
"""
struct Approximation{T,S} <: Function
    original::Function
    domain::Domain
    fun::AbstractRationalInterpolant{T,S}
    allowed::Function
    prenodes::Vector{T}
    test_points::Vector{S}
    history::RFIVector
end

function Approximation(
    f::Function,
    domain::Domain,
    fun::AbstractRationalInterpolant{T,S},
    allowed::Function,
    prenodes::AbstractVector,
    test_points::AbstractVector,
    hist::RFIVector = RFIVector{typeof(fun)}()
    ) where {T,S}
    test = convert(Vector{S}, test_points)
    return Approximation(f, domain, fun, allowed, prenodes, test, hist)
end

(f::Approximation)(z) = f.fun(z)

function Base.show(io::IO, ::MIME"text/plain", f::Approximation)
    print(io, f.fun, " on the domain: ", f.domain)
end

nodes(r::Approximation, args...) = nodes(r.fun, args...)
Base.values(r::Approximation, args...) = values(r.fun, args...)
weights(r::Approximation, args...) = weights(r.fun, args...)
degree(r::Approximation) = degree(r.fun)
poles(F::Approximation) = poles(F.fun)
residues(f::Approximation, args...) = residues(f.fun, args...)
roots(f::Approximation) = roots(f.fun)

#####
##### IMPLEMENTATION
#####

##### Create an approximation on a continuous domain

"""
    approximate(f, domain)

Adaptively compute a rational interpolant on a continuous or discrete domain.

# Arguments
## Continuous domain
- `f::Function`: function to approximate
- `domain`: curve, path, or region from ComplexRegions

## Discrete domain
- `f::Function`: function to approximate
- `z::AbstractVector`: point set on which to approximate

# Keywords
- `max_iter::Integer=150`: maximum number of iterations on node addition
- `float_type::Type`: floating point type to use for the computation¹
- `tol::Real=1000*eps(float_type)`: relative tolerance for stopping
- `allowed::Function`: function to determine if a pole is allowed²
- `refinement::Integer=3`: number of test points between adjacent nodes (continuum only)
- `lookahead::Integer=20`: number of iterations to determine stagnation

¹Default of `float_type` is the promotion of `float(1)` and the float type of the domain.
²Default is to disallow poles on the curve or in the interior of a continuous domain, or to accept all poles on a discrete domain. Use `allowed=true` to allow all poles.

# Returns
- `r::Approximation`: the rational interpolant

See also [`Approximation`](@ref), [`check`](@ref), [`rewind`](@ref).

# Examples
```julia-repl
julia> f = x -> tanh( 40*(x - 0.15) );

julia> r = approximate(f, unit_interval)
Barycentric{Float64, Float64} rational function of type (22, 22) on the domain: Path{Float64} with 1 curve

julia> ( r(0.3), f(0.3) )
(0.9999877116508015, 0.9999877116507956)

julia> check(r);   # accuracy over the domain
[ Info: Max error is 1.58e-13
```
"""
function approximate(f::Function, R::ComplexRegions.AbstractRegion; kw...)
    # ::Function, ::AbstractRegion
    # Given a region as domain, we interpret poles as not being allowed in that region.
    r = approximate(f, R.boundary; allowed=z->!in(z,R), kw...)
    return Approximation(f, R, r.fun, r.allowed, r.prenodes, r.test_nodes, r.history)
end

# ::Function, ::ComplexCurve
approximate(f::Function, d::ComplexCurve; kw...) = approximate(f, Path(d); kw...)

# ::Function, ::ComplexPath
approximate(f::Function, d::ComplexClosedCurve; kw...) = approximate(f, ClosedPath(d); kw...)

# Update the matrices of test points and f values.
function update_test_points!(path, fun, test, τ, fτ, s, Δs, δ, rows, cols)
    # Each column belongs to one of the node parameter values and contains points
    # that are equispaced up to the next one.
    @inbounds @fastmath for i in rows, j in cols
        test[i, j] = s[j] + δ[i] * Δs[j]
        τ[i, j] = point(path, test[i, j])
        fτ[i, j] = fun(τ[i, j])
    end
    return nothing
end

# All continuum methods end up here:
# ::Function, ::ComplexPath
function approximate(f::Function, d::ComplexPath;
    method = Barycentric,
    float_type = promote_type(real_type(d), typeof(float(1))),
    tol = 1000*eps(float_type),
    allowed::Union{Function,Bool} = z -> dist(z, d) > tol,
    max_iter = 150,
    refinement = 3,
    lookahead = 20
    )

    num_ref = 14    # initial number of test points between nodes; decreases to `refinement`
    δ = float_type.((1:num_ref) / (num_ref+1))    # refinement fractions
    if allowed==true
        allowed = z -> true
    end
    s = [zero(float_type)]         # pre-node parameters (except the last)
    Δs = length(d) * [1 - s[1]]    # pre-node spacings
    σ = point(d, s)            # vector of nodes points
    if !isclosed(d)
        # Include both endpoints as nodes for open paths:
        push!(σ, point(d, length(d)*one(float_type)))
    end
    if isreal(d)
        # Enforce reality for a real domain:
        σ = real(σ)
    end
    fσ = f.(σ)         # f at nodes

    # Arrays of test points have one column per node (except the last), num_ref rows
    test = Matrix{float_type}(undef, num_ref, max_iter)    # parameter values of test points
    τ = similar(σ, num_ref, max_iter + 1)             # test points
    fτ = similar(fσ, num_ref, max_iter + 1)           # f at test points
    number_type = promote_type(eltype(τ), eltype(fτ))
    values = similar(complex(fτ))

    # Arrays to track iteration data and progress
    err = float_type[]    # approximation errors
    all_weights = Matrix{number_type}(undef, max_iter+2, max_iter+2)
    lengths = Vector{Int}(undef, max_iter)

    # Initialize test points, data matrices, rational approximation
    update_test_points!(d, f, test, τ, fτ, s, Δs, δ, 1:num_ref, [1])
    data = update_test_values!(method, number_type, num_ref, max_iter)
    r = method(number_type[], number_type[], number_type[])
    idx_new_test = idx_test = CartesianIndices((1:num_ref, 1:1))

    # Add first node
    @views add_nodes!(r, data, τ, fτ, idx_test, σ, fσ)

    # Main iteration
    n, n_max = 1, 1       # iteration counter, all-time max
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

        # Have we succeeded?
        if (last(err) <= tol*fmax) && all(allowed, poles(r))
            break
        end

        # Do we quit?
        plateau = exp(median(log.(last(err, lookahead))))
        stagnant = (n > lookahead) && (plateau < last(err) < 1e-2*fmax)
        if (n == max_iter) || stagnant
            @warn("May not have converged to desired tolerance")
            # Backtrack to last acceptable approximation:
            n += 1
            accepted = false
            while !accepted
                n -= 1
                if n == 0
                    @warn "No acceptable approximation found"
                    n = n_max
                    break
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
        if num_ref > refinement    # initial phase
            num_ref -= 1    # gradually decrease initial refinement level
            δ = float_type.((1:num_ref) / (num_ref+1))
            # Wipe out the old test points and start over:
            update_test_points!(d, f, test, τ, fτ, s, Δs, δ, 1:num_ref, 1:n)
            idx_new_test = idx_test = CartesianIndices((1:num_ref, 1:n))
            r = method(number_type[], number_type[], number_type[])
            @views add_nodes!(r, data, τ, fτ, idx_test, σ, fσ)
        else    # steady-state refinement level
            # Replace one column and add a new column of test points:
            @views add_nodes!(r, data, τ, fτ, idx_test, [σ[end]], [fσ[end]])
            update_test_points!(d, f, test, τ, fτ, s, Δs, δ, 1:num_ref, [jnew, n])
            idx_test = CartesianIndices((1:num_ref, 1:n))
            idx_new_test = idx_test[:, [jnew, n]]
        end
    end

    s = first(s, lengths[n])
    test_points = vec(τ[1:num_ref, 1:n])
    if !isclosed(d)
        # Put the last node back in:
        push!(s, one(float_type))
    end
    hist = RFIVector{typeof(r)}(σ, fσ, all_weights, lengths[1:n_max], n)
    return Approximation(f, d, r, allowed, s, test_points, hist)
end

##### Create an approximation on a discrete domain

# ::Function, ::AbstractVector
function approximate(
    f::Function, z::AbstractVector;
    allowed = z -> true,
    kw...
    )
    r, history, test = approximate(f.(z), z; history=true, kw...)
    return Approximation(f, z, r, allowed, Float64[], test, history)
end

# ::AbstractVector, ::AbstractVector
function approximate(y::AbstractVector{T}, z::AbstractVector{S};
    method = Barycentric,
    float_type = promote_type(real_type(eltype(z)), typeof(float(1))),
    tol = 1000*eps(float_type),
    max_iter = 100,
    lookahead = 16,
    history = false
    ) where {T<:Number,S<:Number}

    m = length(z)
    n = 1    # iteration counter
    fmax = norm(y, Inf)     # scale of f
    number_type = promote_type(eltype(z), eltype(y))
    err = float_type[]

    ȳ = sum(y) / m
    idx_max = argmax(abs(y - ȳ) for y in y)
    σ = [z[idx_max]]
    fσ = [y[idx_max]]
    τ = reshape(z, 1, m)
    fτ = reshape(y, 1, m)
    values = similar(fτ)

    data = update_test_values!(method, number_type, 1, max_iter, m)
    r = method(number_type[], number_type[], number_type[])
    idx_test = CartesianIndex.(1, 1:m)
    deleteat!(idx_test, idx_max)
    idx_new_test = reshape(idx_test, 1, :)
    @views add_nodes!(r, data, τ, fτ, reshape(idx_test, 1, :), σ, fσ)
    lengths = Vector{Int}(undef, max_iter)
    all_weights = Matrix{number_type}(undef, max_iter + 2, max_iter + 2)
    while true
        test_values = update_test_values!(values, r, data, τ, fτ, reshape(idx_test, 1, :), idx_new_test)
        test_actual = transpose(view(fτ, idx_test))
        if any(isnan, test_actual)
            throw(ArgumentError("Function has NaN value at a test point"))
        end
        err_max, idx_max = findmax(abs, test_actual - test_values)
        push!(err, err_max)
        lengths[n] = L = length(nodes(r))
        all_weights[1:L, L] .= weights(r)

        # Are we done?
        if (last(err) <= tol*fmax)    # success
            break
        end

        # Do we quit?
        plateau = exp(median(log.(last(err, lookahead))))
        stagnant = (n > lookahead) && (plateau < last(err) < 1e-2*fmax)
        if (n == max_iter) || stagnant
            @warn("May not have converged to desired tolerance")
            break
        end

        # Add new node:
        push!(σ, τ[idx_test[idx_max[2]]])
        push!(fσ, fτ[idx_test[idx_max[2]]])
        n += 1

        # Update at the test points:
        deleteat!(idx_test, idx_max[2])
        @views add_nodes!(r, data, τ, fτ, idx_test, [σ[end]], [fσ[end]])
        idx_new_test = CartesianIndices((1:1, 1:0))
    end
    if history
        hist = RFIVector{typeof(r)}(σ, fσ, all_weights, lengths[1:n], n)
        return r, hist, vec(τ[idx_test])
    else
        return r
    end
end

##### Helper functions

"""
    rewind(r, index)

Rewind a rational approximation to a state encountered during an iteration.

# Arguments
- `r::Approximation}`: the approximation to rewind
- `index::Integer`: the iteration number to rewind to

# Returns
- the rational function of the specified index (same type as input)

# Examples
```jldoctest
julia> r = approximate(x -> cos(20x), unit_interval)
Barycentric{Float64, Float64} rational interpolant of type (24, 24) on the domain: Path{Float64} with 1 curve

julia> rewind(r, 10)
Barycentric{Float64, Float64} rational interpolant of type (10, 10) on the domain: Path{Float64} with 1 curve
```
"""
function rewind(r::Approximation, idx::Integer)
    if isempty(r.history)
        @error("No convergence history exists.")
    end
    new_rat = r.history[idx]
    return Approximation(r.original, r.domain, new_rat, r.allowed, r.prenodes, r.history)
end

"""
    check(r; quiet=false, prenodes=false)

Check the accuracy of a rational approximation `r` on its domain. Returns the test points and the error at those points.

# Arguments
- `r::Approximation`: rational approximation

# Keywords
- `quiet::Bool=false`: suppress @info output
- `prenodes::Bool=false`: return prenodes of the approximation as well

# Returns
- `τ::Vector`: test points
- `err::Vector`: error at test points

See also [`approximate`](@ref).
"""
function check(F::Approximation; quiet=false, prenodes=false)
    p = F.domain
    if p isa AbstractVector    # discrete domain
        τ = p
        t = collect(eachindex(τ))
    else
        if p isa ComplexSCRegion
            p = p.boundary
        end
        s = sort(F.prenodes)
        if isclosed(p)
            s = [s; length(p)]
        end
        t, τ = refine(p, s, 30, true)
        idx = sortperm(t)
        t = t[idx]
        τ = τ[idx]
    end
    if isreal(nodes(F.fun))
        τ = real(τ)
    end
    err = F.original.(τ) - F.fun.(τ)
    !quiet && @info f"Max error is {norm(err, Inf):.2e}"
    return prenodes ? (t, τ, err) : (τ, err)
end

"""
    get_history(r::Approximation)

Parse the convergence history of a rational approximation.
# Arguments
- `r::Approximation`: the approximation to get the history from
# Returns
- `::Vector`: degrees of the approximations
- `::Vector`: estimated maximum errors of the approximations
- `::Vector{Vector}`: poles of the approximations
- `::Vector{Vector}`: allowed poles of the approximations
- `::Integer`: index of the best approximation

See also [`convergenceplot`](@ref).
"""
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
        rn = hist[idx]
        push!(deg, degree(rn))
        push!(zp, poles(rn))
        push!(allowed, r.allowed.(zp[end]))
        push!(err, maximum(abs, rn.(τ) - fτ) / scale)
    end
    return deg, err, zp, allowed, hist.best
end

# Remove any poles on the domain.  (TODO: Remove?)
function cleanup_poles(f::Approximation, isbad=z->dist(z, f.domain)==0 )
    r = f.fun
    p = filter(!isbad, poles(r))
    res = residues(r, p)
    k = sum(weights(r) .* values(r)) / sum(weights(r))
    pfd(s) = k + sum(a/(s-z) for (a,z) in zip(res, p))
    return pfd
end
