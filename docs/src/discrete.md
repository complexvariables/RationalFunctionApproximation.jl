# Discrete data

For many functions, discretization of the domain is straightforward. But if the function has singularities close to the domain of approximation, a continuum approach may be more robust.

The `approximate` function can take a vector of sample points as a domain. The given function is then evaluated only at those points, and the rational approximation is a fully discrete process that uses only the given data.

```@example mode
using RationalFunctionApproximation, ComplexRegions
x = -1:0.01:1
f = x -> tanh(5 * (x - 0.2))
r = approximate(f, x)
```
You can alternatively provide just the discrete function values yourself. The domain is always given second:

```@example mode
y = f.(x)
r = approximate(y, x)
```

As long as there are no singularities as close to the domain as the sample points are to one another, a basic discretization works well.

```@example mode
I = r.domain
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
