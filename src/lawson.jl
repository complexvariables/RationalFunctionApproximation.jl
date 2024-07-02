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
        _, _, V = svd( sqrt.(wt) .* AB )
        ⍺ = V[1:m, 2m] * fmax             # numerator weights
        β = -V[m+1:2m, 2m]                # denominator weights
        R = (A*⍺) ./ (A*β)
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
"""
function minimax(r::Barycentric, f::Function, nsteps::Integer=20)
    node, val, weight = nodes(r), values(r), weights(r)
    test = sort(refine(node, 20))
    ⍺, β = lawson(test, f.(test), node, val, weight, nsteps)
    return Barycentric(node, ⍺ ./ β, β, ⍺)
end

function minimax(r::Approximation, nsteps::Integer=20)
    s = r.prenodes
    p = isa(r.domain, ComplexSCRegion) ? boundary(r.domain) : r.domain
    _, test = refine(p, s, 20)
    ftest = r.original.(test)
    ⍺, β = lawson(test, ftest, nodes(r.fun), values(r.fun), weights(r.fun), nsteps)
    r̂ = Barycentric(nodes(r.fun), ⍺ ./ β, β, ⍺)
    return Approximation(r.original, r.domain, r̂, r.prenodes)
end
