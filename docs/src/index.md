# Rational function approximation in Julia

*Documentation for [RationalFunctionApproximation.jl](https://github.com/complexvariables/RationalFunctionApproximation.jl).*

This package computes rational approximations of a function or data given in the complex plane. For background reading, see [NakatsukasaAAAAlgorithm2018](@cite), [DriscollAAARational2024a](@cite) (or the related arXiv version [DriscollAAARational2023a](@cite)), [CostaAAAleastSquares2023](@cite), and [SalazarCelisNumericalContinued2024](@cite).

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

More typically, however, a function that is well-behaved on the real axis has a singularity structure lurking in the complex plane, and the poles of rational functions provide a unique way to cope with them. For instance, let's try approximating the hyperbolic secant function:

```@example interval
r = approximate(sech, Segment(-4, 4))
```

The sech function is smooth on the real axis but has poles on the imaginary axis at odd multiples of $i\pi/2$. The rational approximant automatically locates the poles closest to the domain:

```@example interval
2 * poles(r) / π
```

We can use the [DomainColoring](https://eprovst.github.io/DomainColoring.jl/stable/) package to visualize the rational function in the complex plane. Color is used to show the phase angle of the value, while dark-to-bright cycles of lightness show powers of e in the magnitude. The poles stand out as locations of rapid change in both phase and magnitude.

```@example interval
using DomainColoring

domaincolor(r, [-6, 6, -6, 6]; abs=true)
lines!(r.domain, linewidth=3, color=:white)
shg()
```

A meromorphic function such as sech has only those isolated poles as singularities, and getting those right is most of the battle. By contrast, the function $\log(x + 0.05i)$ has a branch point at $x = -0.05i$ necessitating a branch cut connecting it to infinity. A rational approximant uses poles to construct a proxy branch cut:

```@example interval
f = x -> log(x + 0.05im)
r = approximate(f, unit_interval)
domaincolor(r, [-1.2, 1.2, -1.2, 1.2]; abs=true)
lines!(r.domain, linewidth=3, color=:white)
shg()
```

We close this quick introduction with approximation of $|x|$, which has a singularity on the interval. A famous result of Newman in 1964 proved that the best rational approximation of degree $n$ has root-exponential convergence.

```@example interval
r = approximate(abs, unit_interval; tol=1e-12)
convergenceplot(r)
```

We find that the nodes of the approximant are also distributed (nearly) root-exponentially around the singularity:

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

## References

```@bibliography
```
