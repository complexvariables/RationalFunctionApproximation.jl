var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = RationalFunctionApproximation","category":"page"},{"location":"#RationalFunctionApproximation","page":"Home","title":"RationalFunctionApproximation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RationalFunctionApproximation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RationalFunctionApproximation]","category":"page"},{"location":"#RationalFunctionApproximation.Approximation","page":"Home","title":"RationalFunctionApproximation.Approximation","text":"Approximation (type)\n\nApproximation of a function on a domain.\n\nFields\n\noriginal: the original function\ndomain: the domain of the approximation\nfun: the barycentric representation of the approximation\nprenodes: the prenodes of the approximation\n\n\n\n\n\n","category":"type"},{"location":"#RationalFunctionApproximation.Barycentric","page":"Home","title":"RationalFunctionApproximation.Barycentric","text":"Barycentric (type)\n\nBarycentric representation of a rational function.\n\nFields\n\nnode: the nodes of the rational function\nvalue: the values of the rational function\nweight: the weights of the rational function\nwf: the weighted values of the rational function\nstats: convergence statistics\n\n\n\n\n\n","category":"type"},{"location":"#RationalFunctionApproximation.Barycentric-2","page":"Home","title":"RationalFunctionApproximation.Barycentric","text":"Barycentric(node, value, weight, wf=value.*weight; stats=missing)\n\nConstruct a Barycentric rational function.\n\nArguments\n\nnode::AbstractVector: interpolation nodes\nvalue::AbstractVector: values at the interpolation nodes\nweight::AbstractVector: barycentric weights\nwf::AbstractVector: weights times values (optional)\nstats::ConvergenceStatistics`: convergence statistics (optional)\n\nExamples\n\njulia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])\nBarycentric function with 3 nodes and values:\n    1.0=>1.0,  2.0=>2.0,  3.0=>3.0\n\njulia> r(1.5)\n1.5\n\n\n\n\n\n","category":"type"},{"location":"#RationalFunctionApproximation.ConvergenceStats","page":"Home","title":"RationalFunctionApproximation.ConvergenceStats","text":"ConvergenceStats{T}(bestidx, error, nbad, nodes, values, weights, poles)\n\nConvergence statistics for a sequence of rational approximations.\n\nFields\n\nbestidx: the index of the best approximation\nerror: the error of each approximation\nnbad: the number of bad nodes in each approximation\nnodes: the nodes of each approximation\nvalues: the values of each approximation\nweights: the weights of each approximation\npoles: the poles of each approximation\n\nSee also: approximate, Barycentric\n\n\n\n\n\n","category":"type"},{"location":"#Base.values-Tuple{Barycentric}","page":"Home","title":"Base.values","text":"values(r::Barycentric) returns the nodal values of the interpolant as a vector.\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.aaa-Tuple{AbstractVector{<:Number}, AbstractVector{<:Number}}","page":"Home","title":"RationalFunctionApproximation.aaa","text":"aaa(z, y)\naaa(f)\n\nAdaptively compute a rational interpolant.\n\nArguments\n\ndiscrete mode\n\nz::AbstractVector{<:Number}: interpolation nodes\ny::AbstractVector{<:Number}: values at nodes\n\ncontinuous mode\n\nf::Function: function to approximate on the interval [-1,1]\n\nKeyword arguments\n\ndegree::Integer=150: maximum numerator/denominator degree to use\nfloat_type::Type=Float64: floating point type to use for the computation\ntol::Real=1000*eps(float_type): tolerance for stopping\nlookahead::Integer=10: number of iterations to determines stagnation\nstats::Bool=false: return convergence statistics\n\nReturns\n\nr::Barycentric: the rational interpolant\nstats::NamedTuple: convergence statistics, if keyword stats=true\n\nSee also approximate for approximating a function on a region.\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.approximate-Tuple{Function, ComplexRegions.AbstractRegion}","page":"Home","title":"RationalFunctionApproximation.approximate","text":"approximate(f, domain)\n\nAdaptively compute a rational interpolant on a curve, path, or region.\n\nArguments\n\nf::Function: function to approximate\ndomain: curve, path, or region from ComplexRegions\n\nKeyword arguments\n\ndegree::Integer=150: maximum numerator/denominator degree to use\nfloat_type::Type=Float64: floating point type to use for the computation\ntol::Real=1000*eps(float_type): relative tolerance for stopping\nisbad::Function: function to determine if a pole is bad\nrefinement::Integer=3: number of test points between adjacent nodes\nlookahead::Integer=10: number of iterations to determine stagnation\nstats::Bool=false: return convergence statistics with the approximation? (slower)\n\nReturns\n\nr::Approximation: the rational interpolant\n\nSee also Approximation, check, aaa.\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.check-Tuple{RationalFunctionApproximation.Approximation}","page":"Home","title":"RationalFunctionApproximation.check","text":"check(F, f)\n\nCheck the accuracy of a rational approximation F to a function f on its domain.\n\nArguments\n\nF::Approximation: rational approximation\n\nReturns\n\nτ::Vector: test points\nerr::Vector: error at test points\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.degree-Tuple{Barycentric}","page":"Home","title":"RationalFunctionApproximation.degree","text":"degree(r::Barycentric) returns the degree of the numerator and denominator.\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.nodes-Tuple{Barycentric}","page":"Home","title":"RationalFunctionApproximation.nodes","text":"nodes(r::Barycentric) returns the nodes of the interpolant as a vector.\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.rewind-Tuple{Barycentric, Integer}","page":"Home","title":"RationalFunctionApproximation.rewind","text":"rewind(r, degree)\n\nRewind a Barycentric rational function to a lower degree using stored convergence data.\n\nArguments\n\nr::Union{Barycentric,Approximation}: the rational function to rewind\ndegree::Integer: the degree to rewind to\n\nReturns\n\nthe rational function of the specified degree (same type as input)\n\nExamples\n\njulia> r = aaa(x -> cos(20x), stats=true)\nBarycentric function with 25 nodes and values:\n    -1.0=>0.408082,  -0.978022=>0.757786,  -0.912088=>0.820908,  …  1.0=>0.408082\n\njulia> rewind(r, 10)\nBarycentric function with 11 nodes and values:\n    -1.0=>0.408082,  1.0=>0.408082,  -0.466667=>-0.995822,  …  0.898413=>0.636147\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.stats-Tuple{Barycentric}","page":"Home","title":"RationalFunctionApproximation.stats","text":"stats(r::Barycentric) returns the convergence statistics of the interpolant.\n\n\n\n\n\n","category":"method"},{"location":"#RationalFunctionApproximation.weights-Tuple{Barycentric}","page":"Home","title":"RationalFunctionApproximation.weights","text":"weights(r::Barycentric) returns the weights of the interpolant as a vector.\n\n\n\n\n\n","category":"method"}]
}