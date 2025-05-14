# Rational function approximation in Julia

*Documentation for [RationalFunctionApproximation.jl](https://github.com/complexvariables/RationalFunctionApproximation.jl).*

This package computes rational approximations of a function or data given in the complex plane.

A rational function is a ratio of two polynomials. Rational functions are capable of very high accuracy and, unlike polynomial interpolation, do not require the interpolation nodes to be distributed in a highly restricted way. They are a good choice for approximating functions with singularities or other complicated behavior. Also unlike polynomials, however, they do not depend linearly on the data, which has historically made them difficult to compute and work with.

Here's a smooth, gentle function on the interval $[-1, 1]$:

```@example interval
using RationalFunctionApproximation, CairoMakie
CairoMakie.update_theme!(size = (600, 400), fontsize=11)
const shg = current_figure

f = x -> exp(cos(4x) - sin(3x))
lines(-1..1, f)
```

To create a rational function that approximates $f$ well on this domain, we make a call to the `approximate` function:

```@example interval
r = approximate(f, unit_interval)
```

The value of `unit_interval` is defined by the package to be the interval $[-1, 1]$. The result `r` is a type (19,19) rational approximant that can be evaluated like a function:

```@example interval
f(0.5) - r(0.5)
```

We see that this approximation is accurate to about 13 places over the interval:

```@example interval
z, err = check(r)
lines(z, err)
```

The rational approximant interpolates $f$ at nodes that were selected iteratively to represent the function well.

```@example interval
x = nodes(r)
scatter!(x, 0*x, markersize = 8, color=:black)
shg()
```

We could choose to approximate over a wider interval:

```@example interval
using ComplexRegions
r = approximate(f, Segment(-2, 4))
```

Note that the degree of the rational function increased to capture the additional complexity.

One interesting feature of a rational function is that it can have poles, or infinite value, at the roots of the denominator polynomial. In this case, the poles hint at where the function is most sharply peaked:

```@example interval
poleplot(r)
```

Let's try approximating the hyperbolic secant function:

```@example interval
r = approximate(sech, unit_interval)
```

The sech function is smooth on the real axis but has poles at the imaginary axis, and the rational approximant locates the ones closest to the domain:

```@example interval
poles(r) / π
```

We can use the [DomainColoring](https://eprovst.github.io/DomainColoring.jl/stable/) package to visualize the rational function in the complex plane. Color is used to show the phase angle of the value, while patterns of lightness show powers of e in the magnitude of the value. The poles stand out as locations of rapid change in both phase and magnitude.

```@example interval
using DomainColoring

domaincolor(r, [-6, 6, -6, 6]; abs=true)
lines!([-1, 1], [0, 0], linewidth=3, color=:white)
shg()
```

The function $\log(x + 0.05i)$ has a branch point at $x = -0.05i$. A rational approximant uses poles to construct a proxy branch cut:

```@example interval
f = x -> log(x + 0.05im)
r = approximate(f, unit_interval)
domaincolor(r, [-1.2, 1.2, -1.2, 1.2]; abs=true)
lines!([-1, 1], [0, 0], linewidth=3, color=:white)
shg()
```

We close with a function having a singularity that lies on the interval: $|x|$. A famous result of Newman in 1964 proved that the best rational approximation of degree $n$ has root-exponential convergence.

```@example interval
r = approximate(abs, unit_interval; tol=1e-12)
convergenceplot(r)
```

(The errors increase for odd degrees above because they are being measured at the test points discovered at the end of the iteration, not the ones during the iteration.) We find that the nodes of the approximant are also distributed (nearly) root-exponentially around the singularity:

```@example interval
z = filter(>(0), nodes(r))
scatter(sort(abs.(z)), axis=(ylabel="| node |", yscale=log10,))
```

## Contents

* [Algorithms](@ref) describes the algorithms available for rational approximation.
* [Approximation on domains](@ref) shows how to approximate functions on different domains.
* [Discrete data](@ref) shows how to approximate data given as points and values rather than as functions.
* [Minimax approximation](@ref) explains the difference between the default approximation and minimax approximation.
* [Usage from Python](@ref) shows how to use the package from Python.
* [Functions and types](@ref) collects the documentation strings of the major package features.