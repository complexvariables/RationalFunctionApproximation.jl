


# Functions and types {#Functions-and-types}
- [`RationalFunctionApproximation.Approximation`](#RationalFunctionApproximation.Approximation)
- [`RationalFunctionApproximation.Barycentric`](#RationalFunctionApproximation.Barycentric)
- [`RationalFunctionApproximation.Barycentric`](#RationalFunctionApproximation.Barycentric)
- [`RationalFunctionApproximation.DiscretizedPath`](#RationalFunctionApproximation.DiscretizedPath-Tuple{Union{ComplexRegions.AbstractCurve,%20ComplexRegions.AbstractPath},%20AbstractVector})
- [`RationalFunctionApproximation.RFIVector`](#RationalFunctionApproximation.RFIVector)
- [`Base.collect`](#Base.collect)
- [`Base.values`](#Base.values-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.Res`](#RationalFunctionApproximation.Res-Tuple{Function,%20Number})
- [`RationalFunctionApproximation.aaa`](#RationalFunctionApproximation.aaa-Tuple{AbstractVector{<:Number},%20AbstractVector{<:Number}})
- [`RationalFunctionApproximation.add_node!`](#RationalFunctionApproximation.add_node!-Tuple{DiscretizedPath,%20Any})
- [`RationalFunctionApproximation.approximate`](#RationalFunctionApproximation.approximate-Tuple{Function,%20ComplexRegions.AbstractRegion})
- [`RationalFunctionApproximation.check`](#RationalFunctionApproximation.check-Tuple{RationalFunctionApproximation.Approximation})
- [`RationalFunctionApproximation.convergenceplot`](#RationalFunctionApproximation.convergenceplot-Tuple{Any})
- [`RationalFunctionApproximation.decompose`](#RationalFunctionApproximation.decompose-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.degree`](#RationalFunctionApproximation.degree-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.degrees`](#RationalFunctionApproximation.degrees-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.errorplot`](#RationalFunctionApproximation.errorplot-Tuple{Any})
- [`RationalFunctionApproximation.get_history`](#RationalFunctionApproximation.get_history-Union{Tuple{RationalFunctionApproximation.Approximation{T,%20S}},%20Tuple{S},%20Tuple{T}}%20where%20{T,%20S})
- [`RationalFunctionApproximation.minimax`](#RationalFunctionApproximation.minimax)
- [`RationalFunctionApproximation.nodes`](#RationalFunctionApproximation.nodes-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.poleplot`](#RationalFunctionApproximation.poleplot-Tuple{Any})
- [`RationalFunctionApproximation.poles`](#RationalFunctionApproximation.poles-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.poles`](#RationalFunctionApproximation.poles-Union{Tuple{Barycentric{T}},%20Tuple{T}}%20where%20T)
- [`RationalFunctionApproximation.residues`](#RationalFunctionApproximation.residues-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.rewind`](#RationalFunctionApproximation.rewind-Tuple{RationalFunctionApproximation.Approximation,%20Integer})
- [`RationalFunctionApproximation.roots`](#RationalFunctionApproximation.roots-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant})
- [`RationalFunctionApproximation.roots`](#RationalFunctionApproximation.roots-Tuple{Barycentric})

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.Approximation' href='#RationalFunctionApproximation.Approximation'><span class="jlbinding">RationalFunctionApproximation.Approximation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Approximation (type)
```


Approximation of a function on a domain.

**Fields**
- `original`: the original function
  
- `domain`: the domain of the approximation
  
- `fun`: the barycentric representation of the approximation
  
- `allowed`: function to determine if a pole is allowed
  
- `prenodes`: the prenodes of the approximation
  
- `test_points`: test points where residual was computed
  
- `history`: all approximations in the iteration
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/approximate.jl#L7-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.Barycentric' href='#RationalFunctionApproximation.Barycentric'><span class="jlbinding">RationalFunctionApproximation.Barycentric</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Barycentric (type)
```


Barycentric representation of a rational function.

**Fields**
- `node`: the nodes of the rational function
  
- `value`: the values of the rational function
  
- `weight`: the weights of the rational function
  
- `wf`: the weighted values of the rational function
  
- `stats`: convergence statistics
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/barycentric.jl#L2-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.Barycentric-2' href='#RationalFunctionApproximation.Barycentric-2'><span class="jlbinding">RationalFunctionApproximation.Barycentric</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Barycentric(node, value, weight, wf=value .* weight; stats=missing)
```


Construct a `Barycentric` rational function.

**Arguments**
- `node::Vector`: interpolation nodes
  
- `value::Vector`: values at the interpolation nodes
  
- `weight::Vector`: barycentric weights
  

**Keywords**
- `wf::Vector = value .* weight`: weights times values
  

**Returns**
- `::Barycentric`: a barycentric rational interpolating function
  

**Examples**

```julia
julia> r = Barycentric([1, 2, 3], [1, 2, 3], [1/2, -1, 1/2])
Barycentric function with 3 nodes and values:
    1.0=>1.0,  2.0=>2.0,  3.0=>3.0

julia> r(1.5)
1.5
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/barycentric.jl#L39-L64" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.DiscretizedPath-Tuple{Union{ComplexRegions.AbstractCurve, ComplexRegions.AbstractPath}, AbstractVector}' href='#RationalFunctionApproximation.DiscretizedPath-Tuple{Union{ComplexRegions.AbstractCurve, ComplexRegions.AbstractPath}, AbstractVector}'><span class="jlbinding">RationalFunctionApproximation.DiscretizedPath</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
DiscretizedPath(path, s::AbstractVector; kwargs...)
DiscretizedPath(path, n::Integer=0; kwargs...)
```


Discretize a path, keeping the option of future making local refinements.

**Arguments**
- `path`: a ComplexCurve or ComplexPath
  
- `s`: a vector of parameter values
  
- `n`: number of points to discretize the path
  

**Keyword arguments**
- `refinement`: number of refinements to make between consecutive points
  
- `maxpoints`: maximum number of points ever allowed
  

See also [`collect`](/functions#Base.collect), [`add_node!`](/functions#RationalFunctionApproximation.add_node!-Tuple{DiscretizedPath,%20Any}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/utils.jl#L34-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.RFIVector' href='#RationalFunctionApproximation.RFIVector'><span class="jlbinding">RationalFunctionApproximation.RFIVector</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Sequence of rational interpolants produced by an iteration.

**Fields**
- `nodes`: vector of interpolation nodes
  
- `values`: vector of interpolation values
  
- `weights`: matrix of all weights (upper triangle)
  
- `len`: the number of nodes for each approximation
  
- `best`: the index of the best approximation
  

See also: [`approximate`](/functions#RationalFunctionApproximation.approximate-Tuple{Function,%20ComplexRegions.AbstractRegion})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L99-L110" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.collect' href='#Base.collect'><span class="jlbinding">Base.collect</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
collect(d::DiscretizedPath, which=:nodes)
```


Collect the points and parameters of a discretized path.

**Arguments**
- `d`: a DiscretizedPath object
  
- `which`: return the nodes if :nodes, test points if :test, or all if :all
  

**Returns**
- Tuple of two vectors: parameter values and points on the path
  

See also [`add_node!`](/functions#RationalFunctionApproximation.add_node!-Tuple{DiscretizedPath,%20Any}), [`DiscretizedPath`](/functions#RationalFunctionApproximation.DiscretizedPath-Tuple{Union{ComplexRegions.AbstractCurve,%20ComplexRegions.AbstractPath},%20AbstractVector}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/utils.jl#L133-L146" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.values-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#Base.values-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">Base.values</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



values(r) returns a vector of the nodal values of the rational interpolant `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.Res-Tuple{Function, Number}' href='#RationalFunctionApproximation.Res-Tuple{Function, Number}'><span class="jlbinding">RationalFunctionApproximation.Res</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
Res(r, z)
```


Returns the residue of the rational function `r` at the point `z`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L34-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.aaa-Tuple{AbstractVector{<:Number}, AbstractVector{<:Number}}' href='#RationalFunctionApproximation.aaa-Tuple{AbstractVector{<:Number}, AbstractVector{<:Number}}'><span class="jlbinding">RationalFunctionApproximation.aaa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
aaa(y, z)
aaa(f)
```


Adaptively compute a rational interpolant.

**Arguments**

**discrete mode**
- `y::AbstractVector{<:Number}`: values at nodes
  
- `z::AbstractVector{<:Number}`: interpolation nodes
  

**continuous mode**
- `f::Function`: function to approximate on the interval [-1,1]
  

**Keyword arguments**
- `max_degree::Integer=150`: maximum numerator/denominator degree to use
  
- `float_type::Type=Float64`: floating point type to use for the computation
  
- `tol::Real=1000*eps(float_type)`: tolerance for stopping
  
- `stagnation::Integer=10`: number of iterations to determines stagnation
  
- `stats::Bool=false`: return convergence statistics
  

**Returns**
- `r::Barycentric`: the rational interpolant
  
- `stats::NamedTuple`: convergence statistics, if keyword `stats=true`
  

**Examples**

```julia
julia> z = 1im * range(-10, 10, 500);

julia> y = @. exp(z);

julia> r = aaa(z, y);

julia> degree(r)   # both numerator and denominator
12

julia> first(nodes(r), 4)
4-element Vector{ComplexF64}:
 0.0 - 6.272545090180361im
 0.0 + 9.43887775551102im
 0.0 - 1.1022044088176353im
 0.0 + 4.909819639278557im

julia> r(1im * π / 2)
-2.637151617496356e-15 + 1.0000000000000002im
```


See also [`approximate`](/functions#RationalFunctionApproximation.approximate-Tuple{Function,%20ComplexRegions.AbstractRegion}) for approximating a function on a curve or region.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/aaa.jl#L8-L57" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.add_node!-Tuple{DiscretizedPath, Any}' href='#RationalFunctionApproximation.add_node!-Tuple{DiscretizedPath, Any}'><span class="jlbinding">RationalFunctionApproximation.add_node!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
add_node!(d::DiscretizedPath, idx)
```


Add a new node to the discretization, and return the indexes of all affected points. The indexes are valid on the `points` and `params` fields.

**Arguments**
- `d`: a DiscretizedPath object
  
- `idx`: a 2-element tuple, vector, or `CartesianIndex` into the `params` field. This identifies
  

the point to be promoted to a node.

**Returns**
- A 2-element vector of `CartesianIndices` into the `params` and `points` fields.
  

See also [`DiscretizedPath`](/functions#RationalFunctionApproximation.DiscretizedPath-Tuple{Union{ComplexRegions.AbstractCurve,%20ComplexRegions.AbstractPath},%20AbstractVector}), [`collect`](/functions#Base.collect).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/utils.jl#L79-L94" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.approximate-Tuple{Function, ComplexRegions.AbstractRegion}' href='#RationalFunctionApproximation.approximate-Tuple{Function, ComplexRegions.AbstractRegion}'><span class="jlbinding">RationalFunctionApproximation.approximate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
approximate(f, domain)
```


Adaptively compute a rational interpolant on a continuous or discrete domain.

**Arguments**

**Continuous domain**
- `f::Function`: function to approximate
  
- `domain`: curve, path, or region from ComplexRegions
  

**Discrete domain**
- `f::Function`: function to approximate
  
- `z::AbstractVector`: point set on which to approximate
  

**Keywords**
- `max_iter::Integer=150`: maximum number of iterations on node addition
  
- `float_type::Type`: floating point type to use for the computation¹
  
- `tol::Real=1000*eps(float_type)`: relative tolerance for stopping
  
- `allowed::Function`: function to determine if a pole is allowed²
  
- `refinement::Integer=3`: number of test points between adjacent nodes (continuum only)
  
- `stagnation::Integer=20`: number of iterations to determine stagnation
  

¹Default of `float_type` is the promotion of `float(1)` and the float type of the domain. ²Default is to disallow poles on the curve or in the interior of a continuous domain, or to accept all poles on a discrete domain. Use `allowed=true` to allow all poles.

**Returns**
- `r::Approximation`: the rational interpolant
  

See also [`Approximation`](/functions#RationalFunctionApproximation.Approximation), [`check`](/functions#RationalFunctionApproximation.check-Tuple{RationalFunctionApproximation.Approximation}), [`rewind`](/functions#RationalFunctionApproximation.rewind-Tuple{RationalFunctionApproximation.Approximation,%20Integer}).

**Examples**

```julia
julia> f = x -> tanh( 40*(x - 0.15) );

julia> r = approximate(f, unit_interval)
Barycentric{Float64, Float64} rational function of type (22, 22) on the domain: Path{Float64} with 1 curve

julia> ( r(0.3), f(0.3) )
(0.9999877116508015, 0.9999877116507956)

julia> check(r);   # accuracy over the domain
[ Info: Max error is 1.58e-13
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/approximate.jl#L65-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.check-Tuple{RationalFunctionApproximation.Approximation}' href='#RationalFunctionApproximation.check-Tuple{RationalFunctionApproximation.Approximation}'><span class="jlbinding">RationalFunctionApproximation.check</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
check(r; quiet=false, prenodes=false)
```


Check the accuracy of a rational approximation `r` on its domain. Returns the test points and the error at those points.

**Arguments**
- `r::Approximation`: rational approximation
  

**Keywords**
- `quiet::Bool=false`: suppress @info output
  
- `prenodes::Bool=false`: return prenodes of the approximation as well
  

**Returns**
- `τ::Vector`: test points
  
- `err::Vector`: error at test points
  

See also [`approximate`](/functions#RationalFunctionApproximation.approximate-Tuple{Function,%20ComplexRegions.AbstractRegion}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/approximate.jl#L455-L472" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.convergenceplot-Tuple{Any}' href='#RationalFunctionApproximation.convergenceplot-Tuple{Any}'><span class="jlbinding">RationalFunctionApproximation.convergenceplot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convergenceplot(r)
```


Plot the convergence history of a rational approximation.

Markers show the maximum error on (the boundary of) the domain as a function of the numerator/denominator degree. A red marker indicates that the approximation has disallowed poles in its domain. A gold halo highlights the best approximation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/RationalFunctionApproximation.jl#L38-L44" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.decompose-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.decompose-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.decompose</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
decompose(r)
```


Return the roots, poles, and residues of the rational interpolant `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L53-L57" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.degree-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.degree-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.degree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



degree(r) returns the degree of the denominator of the rational `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.degrees-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.degrees-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.degrees</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



degrees(r) returns the degrees of the numerator and denominator of the rational `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.errorplot-Tuple{Any}' href='#RationalFunctionApproximation.errorplot-Tuple{Any}'><span class="jlbinding">RationalFunctionApproximation.errorplot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
errorplot(r; use_abs=false)
```


Plot the pointwise error of an `Approximation` on (the boundary of) its domain. If the error is not real, then the real and imaginary parts are plotted separately, unless `use_abs=true`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/RationalFunctionApproximation.jl#L49-L53" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.get_history-Union{Tuple{RationalFunctionApproximation.Approximation{T, S}}, Tuple{S}, Tuple{T}} where {T, S}' href='#RationalFunctionApproximation.get_history-Union{Tuple{RationalFunctionApproximation.Approximation{T, S}}, Tuple{S}, Tuple{T}} where {T, S}'><span class="jlbinding">RationalFunctionApproximation.get_history</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_history(r::Approximation)
```


Parse the convergence history of a rational approximation.

**Arguments**
- `r::Approximation`: the approximation to get the history from
  

**Returns**
- `::Vector`: degrees of the approximations
  
- `::Vector`: estimated maximum errors of the approximations
  
- `::Vector{Vector}`: poles of the approximations
  
- `::Vector{Vector}`: allowed poles of the approximations
  
- `::Integer`: index of the best approximation
  

See also [`convergenceplot`](/functions#RationalFunctionApproximation.convergenceplot-Tuple{Any}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/approximate.jl#L493-L507" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.minimax' href='#RationalFunctionApproximation.minimax'><span class="jlbinding">RationalFunctionApproximation.minimax</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
minimax(r::Barycentric, f::Function, nsteps::Integer=20)
minimax(r::Approximation, nsteps::Integer=20)
```


Compute an approximately minimax rational approximation to a function `f` on the nodes of a given rational function in barycentric form. The returned approximation has the same type as the first input argument.

The `nsteps` argument controls the number of Lawson iterations. The default value is 20.

**Examples**

```julia
julia> f(x) = tanh( 40*(x - 0.15) );

julia> r = approximate(f, unit_interval, max_degree=8);  # least-squares approximation

julia> check(r);
[ Info: Max error is 1.06e-02

julia> r̂ = minimax(r);

julia> check(r̂);
[ Info: Max error is 1.40e-03
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/lawson.jl#L28-L52" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.nodes-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.nodes-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.nodes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



nodes(r) returns a vector of the interpolation nodes of the rational interpolant.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.poleplot-Tuple{Any}' href='#RationalFunctionApproximation.poleplot-Tuple{Any}'><span class="jlbinding">RationalFunctionApproximation.poleplot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
poleplot(r, idx=0)
```


Plot the domain of the approximation `r` and the poles of the rational approximant. If `idx` is nonzero, it should be an index into the convergence history of `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/RationalFunctionApproximation.jl#L58-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.poles-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.poles-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.poles</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



poles(r) returns the poles of the rational interpolant `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.poles-Union{Tuple{Barycentric{T}}, Tuple{T}} where T' href='#RationalFunctionApproximation.poles-Union{Tuple{Barycentric{T}}, Tuple{T}} where T'><span class="jlbinding">RationalFunctionApproximation.poles</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
poles(r)
```


Return the poles of the rational function `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/barycentric.jl#L104-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.residues-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.residues-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.residues</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
residues(r)
```


Returns two vectors of the poles and residues of the rational function `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L27-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.rewind-Tuple{RationalFunctionApproximation.Approximation, Integer}' href='#RationalFunctionApproximation.rewind-Tuple{RationalFunctionApproximation.Approximation, Integer}'><span class="jlbinding">RationalFunctionApproximation.rewind</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rewind(r, index)
```


Rewind a rational approximation to a state encountered during an iteration.

**Arguments**
- `r::Approximation}`: the approximation to rewind
  
- `index::Integer`: the iteration number to rewind to
  

**Returns**
- the rational function of the specified index (same type as input)
  

**Examples**

```julia
julia> r = approximate(x -> cos(20x), unit_interval)
Barycentric{Float64, Float64} rational interpolant of type (24, 24) on the domain: Path{Float64} with 1 curve

julia> rewind(r, 10)
Barycentric{Float64, Float64} rational interpolant of type (10, 10) on the domain: Path{Float64} with 1 curve
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/approximate.jl#L427-L447" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.roots-Tuple{Barycentric}' href='#RationalFunctionApproximation.roots-Tuple{Barycentric}'><span class="jlbinding">RationalFunctionApproximation.roots</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
roots(r)
```


Return the roots (zeros) of the rational function `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/barycentric.jl#L139-L143" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RationalFunctionApproximation.roots-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}' href='#RationalFunctionApproximation.roots-Tuple{RationalFunctionApproximation.AbstractRationalInterpolant}'><span class="jlbinding">RationalFunctionApproximation.roots</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



roots(r) returns the roots of the rational interpolant `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/750aed7fda32202d74afdb177ea52cb8ea0bf9fd/src/abstract.jl#L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

