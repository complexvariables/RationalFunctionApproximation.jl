var documenterSearchIndex = {"docs":
[{"location":"convergence/#Convergence-of-AAA-rational-approximations","page":"Convergence","title":"Convergence of AAA rational approximations","text":"","category":"section"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"using RationalFunctionApproximation, CairoMakie","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"For a function that is analytic on its domain, the AAA algorithm typically converges at least root-exponentially. In order to observe the convergence, we can construct an approximation that preserves the history of its construction. For example, we approximate an oscillatory function over -11 via:","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"f = x -> cos(11x)\nr = approximate(f, unit_interval, stats=true)\nconvergenceplot(r)","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"In the plot above, the markers show the estimated max-norm error of the n-point AAA rational interpolant over the domain as a function of the numerator/denominator degree n-1. The red circles indicate that a pole of the rational interpolant lies on the interval itself. We can verify this by using rewind to recover the degree-7 approximation that was found along the way to r:","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"r7 = rewind(r, 7)\npoles(r7)","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"When the AAA iteration encounters an approximation with such undesired poles, or having less accuracy than a predecessor, the AAA iteration simply disregards that approximation and continues–-unless there have been more than a designated number of consecutive failures, at which the best interpolant ever encountered is returned. That interpolant is indicated by the gold halo in the convergence plot above.","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"When a singularity is very close to the approximation domain, it can cause stagnation and a large number of bad-pole failures:","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"f = x -> tanh(3000*(x - 1/5))\nr = approximate(f, unit_interval, stats=true)\nconvergenceplot(r)","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"This effect is thought to be mainly due to roundoff and conditioning of the problem. If we use more accurate floating-point arithmetic, we can see that the AAA convergence continues steadily  past the previous plateau:","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"using DoubleFloats\nr = approximate(f, unit_interval, float_type=Double64, stats=true)\nconvergenceplot(r)","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"In the extreme case of a function with a singularity on the domain, the convergence can be substantially affected:","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"f = x -> abs(x - 1/8)\nr = approximate(f, unit_interval, stats=true)\nconvergenceplot(r)","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"In such a case, we might get improvement by increasing the number of allowed consecutive failures via the lookahead keyword argument:","category":"page"},{"location":"convergence/","page":"Convergence","title":"Convergence","text":"r = approximate(f, unit_interval, stats=true, lookahead=20)\nconvergenceplot(r)","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"CurrentModule = RationalFunctionApproximation","category":"page"},{"location":"functions/#Functions-and-types","page":"Functions","title":"Functions and types","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [RationalFunctionApproximation]","category":"page"},{"location":"functions/#RationalFunctionApproximation.Approximation","page":"Functions","title":"RationalFunctionApproximation.Approximation","text":"Approximation (type)\n\nApproximation of a function on a domain.\n\nFields\n\noriginal: the original function\ndomain: the domain of the approximation\nfun: the barycentric representation of the approximation\nprenodes: the prenodes of the approximation\n\n\n\n\n\n","category":"type"},{"location":"functions/#RationalFunctionApproximation.Barycentric","page":"Functions","title":"RationalFunctionApproximation.Barycentric","text":"Barycentric (type)\n\nBarycentric representation of a rational function.\n\nFields\n\nnode: the nodes of the rational function\nvalue: the values of the rational function\nweight: the weights of the rational function\nwf: the weighted values of the rational function\nstats: convergence statistics\n\n\n\n\n\n","category":"type"},{"location":"functions/#RationalFunctionApproximation.Barycentric-2","page":"Functions","title":"RationalFunctionApproximation.Barycentric","text":"Barycentric(node, value, weight, wf=value.*weight; stats=missing)\n\nConstruct a Barycentric rational function.\n\nArguments\n\nnode::AbstractVector: interpolation nodes\nvalue::AbstractVector: values at the interpolation nodes\nweight::AbstractVector: barycentric weights\nwf::AbstractVector: weights times values (optional)\nstats::ConvergenceStatistics`: convergence statistics (optional)\n\nExamples\n\njulia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])\nBarycentric function with 3 nodes and values:\n    1.0=>1.0,  2.0=>2.0,  3.0=>3.0\n\njulia> r(1.5)\n1.5\n\n\n\n\n\n","category":"type"},{"location":"functions/#RationalFunctionApproximation.ConvergenceStats","page":"Functions","title":"RationalFunctionApproximation.ConvergenceStats","text":"ConvergenceStats{T}(bestidx, error, nbad, nodes, values, weights, poles)\n\nConvergence statistics for a sequence of rational approximations.\n\nFields\n\nbestidx: the index of the best approximation\nerror: the error of each approximation\nnbad: the number of bad nodes in each approximation\nnodes: the nodes of each approximation\nvalues: the values of each approximation\nweights: the weights of each approximation\npoles: the poles of each approximation\n\nSee also: approximate, Barycentric\n\n\n\n\n\n","category":"type"},{"location":"functions/#Base.values-Tuple{Barycentric}","page":"Functions","title":"Base.values","text":"values(r) returns the nodal values of the rational interpolant r as a vector.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.aaa-Tuple{AbstractVector{<:Number}, AbstractVector{<:Number}}","page":"Functions","title":"RationalFunctionApproximation.aaa","text":"aaa(z, y)\naaa(f)\n\nAdaptively compute a rational interpolant.\n\nArguments\n\ndiscrete mode\n\nz::AbstractVector{<:Number}: interpolation nodes\ny::AbstractVector{<:Number}: values at nodes\n\ncontinuous mode\n\nf::Function: function to approximate on the interval [-1,1]\n\nKeyword arguments\n\nmax_degree::Integer=150: maximum numerator/denominator degree to use\nfloat_type::Type=Float64: floating point type to use for the computation\ntol::Real=1000*eps(float_type): tolerance for stopping\nlookahead::Integer=10: number of iterations to determines stagnation\nstats::Bool=false: return convergence statistics\n\nReturns\n\nr::Barycentric: the rational interpolant\nstats::NamedTuple: convergence statistics, if keyword stats=true\n\nSee also approximate for approximating a function on a region.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.approximate-Tuple{Function, ComplexRegions.AbstractRegion}","page":"Functions","title":"RationalFunctionApproximation.approximate","text":"approximate(f, domain)\n\nAdaptively compute a rational interpolant on a curve, path, or region.\n\nArguments\n\nf::Function: function to approximate\ndomain: curve, path, or region from ComplexRegions\n\nKeyword arguments\n\nmax_degree::Integer=150: maximum numerator/denominator degree to use\nfloat_type::Type=Float64: floating point type to use for the computation\ntol::Real=1000*eps(float_type): relative tolerance for stopping\nisbad::Function: function to determine if a pole is bad\nrefinement::Integer=3: number of test points between adjacent nodes\nlookahead::Integer=10: number of iterations to determine stagnation\nstats::Bool=false: return convergence statistics with the approximation? (slower)\n\nReturns\n\nr::Approximation: the rational interpolant\n\nSee also Approximation, check, aaa.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.check-Tuple{RationalFunctionApproximation.Approximation}","page":"Functions","title":"RationalFunctionApproximation.check","text":"check(r)\n\nCheck the accuracy of a rational approximation r on its domain.\n\nArguments\n\nr::Approximation: rational approximation\n\nReturns\n\nτ::Vector: test points\nerr::Vector: error at test points\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.decompose-Tuple{Barycentric}","page":"Functions","title":"RationalFunctionApproximation.decompose","text":"decompose(r)\n\nReturn the roots, poles, and residues of the rational function r.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.degree-Tuple{Barycentric}","page":"Functions","title":"RationalFunctionApproximation.degree","text":"degree(r) returns the degree of the numerator and denominator of the rational r.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.nodes-Tuple{Barycentric}","page":"Functions","title":"RationalFunctionApproximation.nodes","text":"nodes(r) returns the nodes of the rational interpolant r as a vector.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.poles-Tuple{RationalFunctionApproximation.Approximation}","page":"Functions","title":"RationalFunctionApproximation.poles","text":"poles(r)\n\nReturn the poles of the rational function r.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.residues-Tuple{RationalFunctionApproximation.Approximation, Vararg{Any}}","page":"Functions","title":"RationalFunctionApproximation.residues","text":"residues(r)\nresidues(r, p=poles(r))\n\nReturn the residues of the rational function r. If a vector p of poles is given, the residues are computed at those locations, preserving order.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.rewind-Tuple{Barycentric, Integer}","page":"Functions","title":"RationalFunctionApproximation.rewind","text":"rewind(r, degree)\n\nRewind a Barycentric rational function to a lower degree using stored convergence data.\n\nArguments\n\nr::Union{Barycentric,Approximation}: the rational function to rewind\ndegree::Integer: the degree to rewind to\n\nReturns\n\nthe rational function of the specified degree (same type as input)\n\nExamples\n\njulia> r = aaa(x -> cos(20x), stats=true)\nBarycentric function with 25 nodes and values:\n    -1.0=>0.408082,  -0.978022=>0.757786,  -0.912088=>0.820908,  …  1.0=>0.408082\n\njulia> rewind(r, 10)\nBarycentric function with 11 nodes and values:\n    -1.0=>0.408082,  1.0=>0.408082,  -0.466667=>-0.995822,  …  0.898413=>0.636147\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.roots-Tuple{RationalFunctionApproximation.Approximation}","page":"Functions","title":"RationalFunctionApproximation.roots","text":"roots(r)\n\nReturn the roots (zeros) of the rational function r.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.stats-Tuple{Barycentric}","page":"Functions","title":"RationalFunctionApproximation.stats","text":"stats(r) returns the convergence statistics of the rational interpolant r.\n\n\n\n\n\n","category":"method"},{"location":"functions/#RationalFunctionApproximation.weights-Tuple{Barycentric}","page":"Functions","title":"RationalFunctionApproximation.weights","text":"weights(r) returns the weights of the rational interpolant r as a vector.\n\n\n\n\n\n","category":"method"},{"location":"mode/#Discrete-vs.-continuous-mode","page":"Discrete vs. continuous","title":"Discrete vs. continuous mode","text":"","category":"section"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"The original AAA algorithm (Nakatsukasa, Sète, Trefethen 2018) works with a fixed set of points on the domain of approximation. The aaa method can work with this type of data:","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"using RationalFunctionApproximation, ComplexRegions\nx = -1:0.01:1\nf = x -> tanh(5 * (x - 0.2))\nr = aaa(x, f.(x))","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"As long as there are no singularities as close to the domain as the sample points are to one another, this fully discrete approach should be fine:","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"I = unit_interval\nprintln(\"nearest pole is $(minimum(dist(z, I) for z in poles(r))) away\")\nxx = -1:0.0005:1\nprintln(\"error = $(maximum(abs, @. f(xx) - r(xx) ))\")","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"But if the distance to a singularity is comparable to the sample spacing, the quality of the approximation may suffer:","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"f = x -> tanh(500 * (x - 0.2))\nr = aaa(x, f.(x))\nprintln(\"nearest pole is $(minimum(dist(z, I) for z in poles(r))) away\")\nprintln(\"error = $(maximum(abs, @. f(xx) - r(xx) ))\")","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"In the continuous mode (Driscoll, Nakatsukasa, Trefethen) used by the approximate method, the samples are refined adaptively to try to ensure that the approximation is accurate everywhere:","category":"page"},{"location":"mode/","page":"Discrete vs. continuous","title":"Discrete vs. continuous","text":"r = approximate(f, I)\nprintln(\"error = $(maximum(abs, @. f(xx) - r(xx) ))\")","category":"page"},{"location":"#Rational-function-approximation-in-Julia","page":"Walkthrough","title":"Rational function approximation in Julia","text":"","category":"section"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Documentation for RationalFunctionApproximation.jl.","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"This package uses the continuous form of the AAA algorithm to adaptively compute rational approximations of functions on intervals and other domains in the complex plane.  See AAA rational approximation on a continuum, which is to appear in SISC.","category":"page"},{"location":"#Approximation-on-[-1,-1]","page":"Walkthrough","title":"Approximation on [-1, 1]","text":"","category":"section"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Here's a smooth, gentle function on the interval -1 1:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"using RationalFunctionApproximation, CairoMakie\nconst shg = current_figure\nf = x -> exp(cos(4x) - sin(3x))\nlines(-1..1, f)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"To create a rational function that approximates f well on this domain, we use the continuous form of the AAA algorithm:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"r = aaa(f)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"The result is a type (19,19) rational approximant that can be evaluated like a function:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"f(0.5) - r(0.5)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"We see that this approximation has more than 13 accurate over most of the interval:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"lines(-1..1, x -> f(x)-r(x))","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"The rational approximant interpolates f at greedily selected nodes:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"x = nodes(r)\nscatter!(x, 0*x, markersize = 12, color=:black)\nshg()","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Here's another smooth example, the hyperbolic secant function:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"f = sech\nr = aaa(f)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"We can verify that this is accurate to 14 digits:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"x = range(-1, 1, 1000)\nextrema(f.(x) - r.(x))","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Since the sech function has poles in the complex plane, the rational approximant r will have corresponding poles:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"using DomainColoring\ndomaincolor(r, [-8, 8], abs=true)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"The poles closest to the interval are found to about 10 digits, while more distant ones are less accurate:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"poles(r) / π","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Here's an example with a more interesting structure of poles and zeros:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"f = x -> tanh(10*(x - 0.1)^2)\ndomaincolor(aaa(f), abs=true)","category":"page"},{"location":"#Approximation-on-other-domains","page":"Walkthrough","title":"Approximation on other domains","text":"","category":"section"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"The AAA algorithm can also be used to approximate functions on other domains as defined in the ComplexRegions package. For example, here's a function defined on the unit circle:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"using RationalFunctionApproximation, CairoMakie, DomainColoring\nconst shg = current_figure\nf = z -> (z^3 - 1) / sin(z - 0.9 - 1im)\nr = approximate(f, unit_circle)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"This approximation is accurate to 13 digits, as we can see by plotting the error around the circle:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"errorplot(r)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Here is how the approximation looks in the complex plane (using black crosses to mark the poles):","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"using ComplexRegions, ComplexPlots\ndomaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true)\nlines!(unit_circle, color=:white, linewidth=5)\nscatter!(poles(r), markersize=18, color=:black, marker=:xcross)\nlimits!(-1.5, 1.5, -1.5, 1.5)\nshg()","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Above, you can also see the zeros at roots of unity.","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"This function has infinitely many poles and an essential singularity inside the unit disk:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"f = z -> tan(1 / z^4)\nr = approximate(f, unit_circle)\ndomaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true)\nlines!(unit_circle, color=:white, linewidth=5)\nshg()","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"We can request an approximation that is analytic in a region. In this case, it would not make sense to request one on the unit disk, since the singularities are necessary:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"r = approximate(f, unit_disk)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"In the result above, the approximation is simply a constant function, as the algorithm could do no better. However, if we request analyticity in the region exterior to the circle, everything works out:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"r = approximate(f, exterior(unit_circle))\nz, err = check(r)\nmaximum(abs, err)","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"We are not limited to intervals and circles! ","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"import ComplexRegions.Shapes\nr = approximate(z -> log(0.35 + 0.4im - z), interior(Shapes.cross))\ndomaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true)\nlines!(boundary(r.domain), color=:white, linewidth=5)\nshg()","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"c = Shapes.hypo(5)\nr = approximate(z -> (z+4)^(-3.5), interior(c))\ndomaincolor(r, [-5, 5, -5, 5], abs=true)\nlines!(c, color=:white, linewidth=5)\nshg()","category":"page"},{"location":"#Unbounded-domains","page":"Walkthrough","title":"Unbounded domains","text":"","category":"section"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"It's also possible to approximate on unbounded domains, but this capability is not yet automated. For example, the function","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"f = z -> 1 / sqrt(z - (-1 + 3im))","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"is analytic on the right half of the complex plane. In order to produce an approximation on that domain, we can transplant it to the unit disk via a Möbius transformation phi:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"φ = Mobius( [-1, -1im, 1], [1im, 0, -1im])   # unit circle to imag axis\nz = discretize(unit_circle, ds=.02)\nfig, ax, _ = scatter(z, axis=(autolimitaspect=1, ))\nax.title = \"z\"\nax, _ = scatter(fig[1,2], φ.(z))\nax.title = \"φ(z)\"\nlimits!(-4, 4, -4, 4)\nshg()","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"By composing f with phi, we can approximate within the disk while f is evaluated only on its native domain:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"r = approximate(f ∘ φ, interior(unit_circle))\ndomaincolor(r, [-2, 2, -2, 2], abs=true)\nlines!(unit_circle, color=:white, linewidth=5)\nscatter!(nodes(r.fun), color=:black, markersize=10)\nshg()","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"Above, the black markers show the nodes of the interpolant. We can view the same approximation within the right half-plane by composing r with phi^-1:","category":"page"},{"location":"","page":"Walkthrough","title":"Walkthrough","text":"φ⁻¹ = inv(φ)\ndomaincolor(r ∘ φ⁻¹, [-8, 8, -8, 8], abs=true)\nlines!([(0, 8), (0, -8)], color=:white, linewidth=5)\nscatter!(φ.(nodes(r.fun)), color=:black, markersize=10)\nlimits!(-8, 8, -8, 8)\nshg()","category":"page"},{"location":"minimax/#Minimax-approximation","page":"Minimax","title":"Minimax approximation","text":"","category":"section"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"The AAA algorithm used by aaa and approximate minimizes error in a discrete least-squares sense. Using an iteratively reweighted least squares (IRLS) approach initially due to Lawson, we can approach the classical problem of optimization in the infinity- or max-norm sense instead.","category":"page"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"For example, suppose we limit the degree of the rational interpolant of a smooth function:","category":"page"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"using RationalFunctionApproximation, CairoMakie\nconst shg = current_figure\nf = x -> exp(cos(4x) - sin(3x))\nr = approximate(f, unit_interval, max_degree=10)\nerrorplot(r)","category":"page"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"Now we apply 20 Lawson iterations to approach the minimax approximation:","category":"page"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"r = minimax(r, 20)\nerrorplot(r)","category":"page"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"As you can see above, the error is now nearly equioscillatory over the interval. Moreover, the interpolation nodes appear to have shifted to resemble Chebyshev points of the first kind. If we try minimax approximation on the unit circle, however, equioscillation tends to lead to equally spaced nodes:","category":"page"},{"location":"minimax/","page":"Minimax","title":"Minimax","text":"f = z -> cos(4z) - sin(3z)\nr = approximate(f, unit_circle, max_degree=10)\nr = minimax(r, 20)\nerrorplot(r, use_abs=false)","category":"page"}]
}
