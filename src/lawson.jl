# Applies the Lawson iteration:
# iteratively reweighting the nodes to minimize the maximum error
function lawson(test, ftest, node, fnode, weight, nsteps)
    n, m = length(test), length(node)
    wt = ones(m + n)
    C = [ 1/(x-s) for x in test, s in node ]
    R = (C * (weight .* fnode)) ./ (C * weight)
    F = [ftest; fnode];  fmax = norm(F,Inf)
    R = [R; fnode]
    A = [C; I];                                 # numerator
    B = [(ftest / fmax) .* C; diagm(fnode / fmax)]    # f*denominator
    AB = [A B]
    ⍺, β = [], []    # put into scope
    for _ in 1:nsteps
        _, _, V =  svd( sqrt.(wt) .* AB )
        ⍺ = V[1:m, 2m] * fmax             # numerator weights
        β = -V[m+1:2m, 2m]                # denominator weights
        denom = A*β
        if any(iszero, denom)
            break    # TODO: why does this happen?
        end
        R = (A*⍺) ./ denom
        wt = normalize(wt .* abs.(F - R), Inf)  # iterative reweighting
    end
    return ⍺, β
end

"""
    minimax(r::Barycentric, f::Function, nsteps::Integer=20)
    minimax(r::Approximation, nsteps::Integer=20)

Compute an approximately minimax rational approximation to a function `f` on the nodes of a
given rational function in barycentric form. The returned approximation has the same type
as the first input argument.

The `nsteps` argument controls the number of Lawson iterations. The default value is 20.

# Examples
```julia-repl
julia> f(x) = tanh( 40*(x - 0.15) );

julia> r = approximate(f, unit_interval, max_degree=8);  # least-squares approximation

julia> check(r);
[ Info: Max error is 1.06e-02

julia> r̂ = minimax(r);

julia> check(r̂);
[ Info: Max error is 1.40e-03
```
"""
function minimax(r::Barycentric, f::Function, nsteps::Integer=20)
    node, val, weight = nodes(r), values(r), weights(r)
    test = sort(refine(node, 20))
    ⍺, β = lawson(test, f.(test), node, val, weight, nsteps)
    return Barycentric(node, ⍺ ./ β, β, ⍺)
end

function minimax(r::Approximation, nsteps::Integer=20)
    if !isa(r.fun, Barycentric)
        error("`minimax`` only works with barycentric approximations")
    end
    p = DiscretizedPath(r.path, 20)
    _, test = collect(p, :test)
    ftest = r.original.(test)
    ⍺, β = lawson(test, ftest, nodes(r.fun), values(r.fun), weights(r.fun), nsteps)
    r̂ = Barycentric(nodes(r.fun), ⍺ ./ β, β, ⍺)
    return Approximation(r.original, r.domain, r̂, r.allowed, r.path)
end
