# Discrete data

The original AAA algorithm ([Nakatsukasa, Sète, Trefethen 2018](https://epubs.siam.org/doi/abs/10.1137/16M1106122)) works with a fixed set of points on the domain of approximation. There is a legacy `aaa` function that can work with this type of data:

```@example mode
using RationalFunctionApproximation, ComplexRegions
x = -1:0.01:1
f = x -> tanh(5 * (x - 0.2))
r = aaa(x, f.(x))
```

However, it's preferable to use the `approximate` function for this purpose, as the result type is more useful within the package. Simply pass the function and, in the form of a vector, the domain.

```@example mode
r = approximate(f, x)
```

As long as there are no singularities as close to the domain as the sample points are to one another, this fully discrete approach should be fine:

```@example mode
I = unit_interval
println("nearest pole is $(minimum(dist(z, I) for z in poles(r))) away")
_, err = check(r);
println("max error on the given domain: ", maximum(abs, err))
```

But if the distance to a singularity is comparable to the sample spacing, the quality of the approximation may suffer. Even worse, the method may not be aware that it has failed.

```@example mode
f = x -> tanh(400 * (x - 0.2))
r = approximate(f, x)
println("nearest pole is $(minimum(dist(z, I) for z in poles(r))) away")
_, err = check(r);
println("max error on the given domain: ", maximum(abs, err))
err = maximum(abs(f(x)- r(x)) for x in range(-1, 1, 3000))
println("max error on finer test points: ", err)
``` 

In the continuous mode, the adaptive sampling of the domain attempts to ensure that the approximation is accurate everywhere.

```@example mode
r = approximate(f, I; tol=1e-12)
err = maximum(abs(f(x)- r(x)) for x in range(-1, 1, 3000))
println("max error on finer test points: ", err)
```