
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

##
# Decompose into partial fractions data
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

residues(F::Approximation, args...) = residues(F.fun, args...)
residues(r::Barycentric) = residues(r, poles(r))
function residues(r::Barycentric, pol)
    numer = t -> sum( w*y / (t-z) for (z, y, w) in zip(nodes(r), values(r), weights(r)))
    denomdiff = t -> -sum( w / (t-z)^2 for (z, w) in zip(nodes(r), weights(r)))
    res = similar( complex(pol) )
    res .= numer.(pol) ./ denomdiff.(pol)
    return res
end

roots(F::Approximation) = roots(F.fun)
function roots(r::Barycentric)
    wf = r.w_times_f
    m = length(wf)
    B = diagm( [0; ones(m)] )
    # Thanks to Daan Huybrechs
    E = [0 transpose(wf); ones(m) diagm(nodes(r))]
    return filter(isfinite, eigvals(E, B))
end

# remove poles on the curve
function cleanup_poles(F::Approximation, isbad=z->dist(z, F.domain)==0 )
    r = F.fun
    k = sum(weights(r) .* values(r)) / sum(weights(r))
    res = residues(r, pol)
    good = @. !isbad(pol)
    data = zip(res[good], pol[good])
    pfd(s) = k + sum(a/(s-z) for (a,z) in data)
    return pfd
end

function decompose(r::Barycentric)
    p = poles(r)
    return roots(r),p,residues(r,p)
end

function diffmat(r::Barycentric{T}) where T
    x, y, w = nodes(r), values(r), weights(r)
    n = length(x)
    entry(i,j) = i==j ? T(0) : w[j] / ( w[i] * (x[i] - x[j]) )
    D = [ entry(i,j) for i in 1:n, j in 1:n ]
    for i in 1:n
        D[i,i] = -sum(D[i,:])
    end
    return D
end

# needs to recompute weights
Base.diff(r::Barycentric) = Barycentric(nodes(r), diffmat(r)*values(r), weights(r))

# needs to recompute weights
function integrate(r::Barycentric)
    x = nodes(r)
    n = length(x)
    j = argmin(x)
    D = diffmat(r)
    idx = [1:j-1; j+1:n]
    y = zeros(eltype(x), n)σ
    y[idx] = D[idx, idx] \ values(r)[idx]
    return Barycentric(x, y, weights(r))
end
