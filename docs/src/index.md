# Rational function approximation in Julia

*Documentation for [RationalFunctionApproximation.jl](https://github.com/complexvariables/RationalFunctionApproximation.jl).*

This package uses the continuous form of the AAA algorithm to adaptively compute rational approximations of functions on intervals and other domains in the complex plane.  See [AAA rational approximation on a continuum](https://arxiv.org/abs/2305.03677), which is to appear in SISC.

## Approximation on [-1, 1]

Here's a smooth, gentle function on the interval $[-1, 1]$:

```@example interval
using RationalFunctionApproximation, CairoMakie
CairoMakie.update_theme!(size = (400, 250), fontsize=11)
const shg = current_figure

f = x -> exp(cos(4x) - sin(3x))
lines(-1..1, f)
```
 To create a rational function that approximates $f$ well on this domain, we use the continuous form of the AAA algorithm:

```@example interval
r = aaa(f)
```
The result is a type (19,19) rational approximant that can be evaluated like a function:

```@example interval
f(0.5) - r(0.5)
```

We see that this approximation has more than 13 accurate digits over most of the interval:

```@example interval
lines(-1..1, x -> f(x)-r(x))
```

The rational approximant interpolates $f$ at greedily selected nodes:
    
```@example interval
x = nodes(r)
scatter!(x, 0*x, markersize = 8, color=:black)
shg()
```

Here's another smooth example, the hyperbolic secant function:

```@example interval
f = sech
r = aaa(f)
```

We can verify that this is accurate to 14 digits:

```@example interval
x = range(-1, 1, 1000)
extrema(f.(x) - r.(x))
```

Since the sech function has poles in the complex plane, the rational approximant $r$ will have corresponding poles:

```@example interval
using DomainColoring
CairoMakie.update_theme!(size = (360, 360))

domaincolor(r, [-8, 8, -8, 8], abs=true)
```

The poles closest to the interval are found to about 10 digits, while more distant ones are less accurate:

```@example interval
poles(r) / π
```

Here's an example with a more interesting structure of poles and zeros:

```@example interval
f = x -> tanh(10*(x - 0.1)^2)
domaincolor(aaa(f), abs=true)
```

