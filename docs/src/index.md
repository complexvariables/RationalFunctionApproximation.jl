# Rational function approximation in Julia

*Documentation for [RationalFunctionApproximation.jl](https://github.com/complexvariables/RationalFunctionApproximation.jl).*

This package computes rational approximations of a function or data given in a domain in the complex plane (including real intervals).

A rational function is a ratio of two polynomials. Rational functions are capable of very high accuracy, usually requiring fewer degrees of freedom than polynomial approximations. They are an especially good choice for approximating functions with singularities near the domain of approximation. Rational functions have applications to finding singularities and roots, computing derivatives, continuation across boundaries, solving PDEs, and more; for a survey of applications, see [NakatsukasaApplicationsAAA2026](@cite). For background on the algorithms in this package, see [NakatsukasaAAAAlgorithm2018](@cite), [DriscollAAARational2024](@cite) (or the related arXiv version [DriscollAAARational2023a](@cite)), [CostaAAAleastSquares2023](@cite), and [SalazarCelisNumericalContinued2024](@cite).

## Basic walkthrough

Here's a smooth, gentle function on the interval $[-1, 1]$:

```@example interval
using CairoMakie
CairoMakie.update_theme!(size = (600, 400), fontsize=11)
const shg = current_figure

f(x) = exp(cos(4x) - sin(3x))
lines(-1..1, f)
```

To create a rational function that approximates $f$, we make a call to the `approximate` function. By default, it uses the real interval $[-1,1]$ as the domain.

```@repl interval
using RationalFunctionApproximation
r = approximate(f)
```

The result `r` is a type (19, 18) rational approximant that can be evaluated like a function:

```@repl interval
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

One important feature of a rational function is that it can have poles, or infinite value, at the roots of the denominator polynomial. In this case, the poles hint at where the function is most sharply peaked:

```@example interval
poleplot(r)
```

As another example, the sech function is smooth on the real axis but has poles on the imaginary axis at odd multiples of $i\pi/2$. An approximation over $[-4,4]$ accurately locates the poles closest to the real axis:

```@example interval
r = approximate(sech, -4..4)
2 * poles(r) / π
```

We can use the [DomainColoring](https://eprovst.github.io/DomainColoring.jl/stable/) package to visualize the rational function in the complex plane. Color is used to show the phase angle of the value, while dark-to-bright cycles of lightness show powers of e in the magnitude. The poles stand out as locations of rapid change in both phase and magnitude.

```@example interval
using DomainColoring

domaincolor(r, 6; abs=true)
lines!(r.domain, linewidth=3, color=:white)
shg()
```

A meromorphic function such as sech has only those isolated poles as singularities, and getting those right is most of the battle. By contrast, the function $\log(x + 0.05i)$ has a branch point at $x = -0.05i$ necessitating a branch cut connecting it to infinity. A rational approximant uses poles to construct a proxy branch cut:

```@example interval
f(x) = log(x + 0.05im)
r = approximate(f)
domaincolor(r, 1.2; abs=true)
lines!(r.domain, linewidth=3, color=:white)
shg()
```

We close this quick introduction with approximation of $|x|$, which has a singularity on the interval. A famous result of Newman in 1964 proved that the best rational approximation of degree $n$ has root-exponential convergence. In order to get the most from the approximation, we need to tell the constructor to be stubborn about declaring the iteration stagnated.

```@example interval
r = approximate(abs; tol=1e-12, stagnation=50)
convergenceplot(r)
```

We find that the nodes of the approximant are also distributed (nearly) root-exponentially around the singularity:

```@example interval
z = filter(>(0), nodes(r))
scatter(sort(z), axis=(xlabel="index", ylabel="node location", xscale=sqrt, yscale=log10,))
```

## Feedback and contributions

If you encounter problems or want to make suggestions, feel free to [file an issue on GitHub](https://github.com/complexvariables/RationalFunctionApproximation.jl/issues).

If you want to contribute to the package, please see [the guidelines](https://github.com/complexvariables/RationalFunctionApproximation.jl/blob/main/CONTRIBUTING.md).

## Contents

* [Installation](@ref) describes the algorithms available for rational approximation.
* [Algorithms](@ref) describes the algorithms available for rational approximation.
* [Approximation on domains](@ref) shows how to approximate functions on different domains.
* [Discrete data](@ref) shows how to approximate data given as points and values rather than as functions.
* [Minimax approximation](@ref) explains the difference between the default approximation and minimax approximation.
* [Usage from Python](@ref) shows how to use the package from Python.
* [Functions and types](@ref) collects the documentation strings of the major package features.

## References

```@bibliography
```
