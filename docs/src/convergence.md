# Convergence of AAA rational approximations

```@example convergence
using RationalFunctionApproximation, CairoMakie
```

For a function that is analytic on its domain, the AAA algorithm typically converges at least root-exponentially. In order to observe the convergence, we can construct an approximation that preserves the history of its construction. For example, we approximate an oscillatory function over $[-1,1]$ via:

```@example convergence
f = x -> cos(11x)
r = approximate(f, unit_interval, stats=true)
convergenceplot(r)
```
In the plot above, the markers show the estimated max-norm error of the $n$-point AAA rational interpolant over the domain as a function of the numerator/denominator degree $n-1$. The red circles indicate that a pole of the rational interpolant lies on the interval itself. We can verify this by using `rewind` to recover the degree-7 approximation that was found along the way to `r`:

```@repl convergence
r7 = rewind(r, 7)
poles(r7)
```
When the AAA iteration encounters an approximation with such undesired poles, or having less accuracy than a predecessor, the AAA iteration simply disregards that approximation and continues---unless there have been more than a designated number of consecutive failures, at which the best interpolant ever encountered is returned. That interpolant is indicated by the gold halo in the convergence plot above.

When a singularity is very close to the approximation domain, it can cause stagnation and a large number of bad-pole failures:

```@example convergence
f = x -> tanh(3000*(x - 1//5))
r = approximate(f, unit_interval, stats=true)
convergenceplot(r)
```

This effect is thought to be mainly due to roundoff and conditioning of the problem. If we use more accurate floating-point arithmetic, we can see that the AAA convergence continues steadily past the previous plateau. In the following, we apply `Double64` arithmetic, having used exact rational numbers already in the definition of `f`:

```@example convergence
using DoubleFloats, ComplexRegions
r = approximate(f, Segment{Double64}(-1, 1), stats=true)
convergenceplot(r)
```

In the extreme case of a function with a singularity on the domain, the convergence can be substantially affected:

```@example convergence
f = x -> abs(x - 1/8)
r = approximate(f, unit_interval, stats=true)
convergenceplot(r)
```

In such a case, we might get improvement by increasing the number of allowed consecutive failures via the `lookahead` keyword argument:

```@example convergence
r = approximate(f, unit_interval, stats=true, lookahead=20)
convergenceplot(r)
```