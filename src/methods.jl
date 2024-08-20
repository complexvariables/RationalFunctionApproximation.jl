"""
    r(z)
    evaluate(r, z)

Evaluate the rational function at `z`.
"""

(r::Barycentric)(z) = evaluate(r, z)
(f::Approximation)(z) = evaluate(f.fun, z)

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
poles(F::Approximation) = poles(F.fun)
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
residues(f::Approximation, args...) = residues(f.fun, args...)
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
roots(f::Approximation) = roots(f.fun)
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

"""
    decompose(r)

Return the roots, poles, and residues of the rational function `r`.
"""
function decompose(r::Barycentric)
    p = poles(r)
    return roots(r), p, residues(r,p)
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
