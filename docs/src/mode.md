# Discrete vs. continuous mode

The original AAA algorithm ([Nakatsukasa, Sète, Trefethen 2018](https://epubs.siam.org/doi/abs/10.1137/16M1106122)) works with a fixed set of points on the domain of approximation. The `aaa` method can work with this type of data:

```@example mode
using RationalFunctionApproximation, ComplexRegions
x = -1:0.01:1
f = x -> tanh(5 * (x - 0.2))
r = aaa(x, f.(x))
```

As long as there are no singularities as close to the domain as the sample points are to one another, this fully discrete approach should be fine:

```@example mode
I = unit_interval
println("nearest pole is $(minimum(dist(z, I) for z in poles(r))) away")
xx = -1:0.0005:1
println("error = $(maximum(abs, @. f(xx) - r(xx) ))")
```

But if the distance to a singularity is comparable to the sample spacing, the quality of the approximation may suffer:

```@example mode
f = x -> tanh(500 * (x - 0.2))
r = aaa(x, f.(x))
println("nearest pole is $(minimum(dist(z, I) for z in poles(r))) away")
println("error = $(maximum(abs, @. f(xx) - r(xx) ))")
``` 

In the continuous mode ([Driscoll, Nakatsukasa, Trefethen](https://arxiv.org/abs/2305.03677)) used by the `approximate` method, the samples are refined adaptively to try to ensure that the approximation is accurate everywhere:

```@example mode
r = approximate(f, I)
println("error = $(maximum(abs, @. f(xx) - r(xx) ))")
```