# Minimax approximation

The AAA algorithm used by `aaa` and `approximate` minimizes error in a discrete least-squares sense. Using an iteratively reweighted least squares (IRLS) approach initially due to Lawson, we can approach the classical problem of optimization in the infinity- or max-norm sense instead.

For example, suppose we limit the degree of the rational interpolant of a smooth function:

```@example minimax
using RationalFunctionApproximation, CairoMakie
const shg = current_figure
f = x -> exp(cos(4x) - sin(3x))
r = approximate(f, unit_interval, degree=10)
errorplot(r)
```

Now we apply 20 Lawson iterations to approach the minimax approximation:

```@example minimax
r = minimax(r, 20)
errorplot(r)
```

As you can see above, the error is now nearly equioscillatory over the interval. Moreover, the interpolation nodes appear to have shifted to resemble Chebyshev points of the first kind. If we try minimax approximation on the unit circle, however, equioscillation tends to lead to equally spaced nodes:

```@example minimax
f = z -> cos(4z) - sin(3z)
r = approximate(f, unit_circle, degree=10)
r = minimax(r, 20)
errorplot(r, use_abs=false)
```
