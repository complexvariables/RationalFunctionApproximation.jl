# RationalFunctionApproximation

Documentation for [RationalFunctionApproximation](https://github.com/complexvariables/RationalFunctionApproximation.jl).

This package uses the continuous form of the AAA algorithm to adaptively compute rational approximations of functions on intervals and other domains in the complex plane.  See [AAA rational approximation on a continuum](https://arxiv.org/abs/2305.03677), which is to appear in SISC.

## Approximation on [-1, 1]

Here's a smooth, gentle function on the interval $[-1, 1]$:

```@example interval
using RationalFunctionApproximation, CairoMakie
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

We see that this approximation has more than 13 accurate over most of the interval:

```@example interval
lines(-1..1, x -> f(x)-r(x))
```

The rational approximant interpolates $f$ at greedily selected nodes:
    
```@example interval
x = nodes(r)
scatter!(x, 0*x, markersize = 12, color=:black)
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
domaincolor(r, [-8, 8], abs=true)
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

## Approximation on other domains

The AAA algorithm can also be used to approximate functions on other domains as defined in the [`ComplexRegions`](https://complexvariables.github.io/ComplexRegions.jl/stable/) package. For example, here's a function defined on the unit circle:

```@example circle
using RationalFunctionApproximation, CairoMakie, DomainColoring
const shg = current_figure
f = z -> (z^3 - 1) / sin(z - 0.9 - 1im)
r = approximate(f, unit_circle)
```

This approximation is accurate to 13 digits, as we can see by plotting the error around the circle:

```@example circle
errorplot(r)
```

Here is how the approximation looks in the complex plane (using black crosses to mark the poles):

```@example circle
using ComplexRegions, ComplexPlots
domaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true)
lines!(unit_circle, color=:white, linewidth=5)
scatter!(poles(r), markersize=18, color=:black, marker=:xcross)
limits!(-1.5, 1.5, -1.5, 1.5)
shg()
```

Above, you can also see the zeros at roots of unity.

This function has infinitely many poles and an essential singularity inside the unit disk:

```@example circle
f = z -> tan(1 / z^4)
r = approximate(f, unit_circle)
domaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true)
lines!(unit_circle, color=:white, linewidth=5)
shg()
```

We can request an approximation that is analytic in a region. In this case, it would not make sense to request one on the unit disk, since the singularities are necessary:

```@example circle
r = approximate(f, unit_disk)
```

In the result above, the approximation is simply a constant function, as the algorithm could do no better. However, if we request analyticity in the region exterior to the circle, everything works out:

```@example circle
r = approximate(f, exterior(unit_circle))
z, err = check(r)
maximum(abs, err)
```

We are not limited to intervals and circles! 

```@example circle
import ComplexRegions.Shapes
r = approximate(z -> log(0.35 + 0.4im - z), interior(Shapes.cross))
domaincolor(r, [-1.5, 1.5, -1.5, 1.5], abs=true)
lines!(boundary(r.domain), color=:white, linewidth=5)
shg()
```

```@example circle
c = Shapes.hypo(5)
r = approximate(z -> (z+4)^(-3.5), interior(c))
domaincolor(r, [-5, 5, -5, 5], abs=true)
lines!(c, color=:white, linewidth=5)
shg()
```

## Unbounded domains

It's also possible to approximate on unbounded domains, but this capability is not yet automated. For example, the function

```@example circle
f = z -> 1 / sqrt(z - (-1 + 3im))
```

is analytic on the right half of the complex plane. In order to produce an approximation on that domain, we can transplant it to the unit disk via a Möbius transformation $\phi$:

```@example circle
φ = Mobius( [-1, -1im, 1], [1im, 0, -1im])   # unit circle to imag axis
z = discretize(unit_circle, ds=.02)
fig, ax, _ = scatter(z, axis=(autolimitaspect=1, ))
ax.title = "z"
ax, _ = scatter(fig[1,2], φ.(z))
ax.title = "φ(z)"
limits!(-4, 4, -4, 4)
shg()
``` 

By composing $f$ with $\phi$, we can approximate within the disk while $f$ is evaluated only on its native domain:

```@example circle
r = approximate(f ∘ φ, interior(unit_circle))
domaincolor(r, [-2, 2, -2, 2], abs=true)
lines!(unit_circle, color=:white, linewidth=5)
scatter!(nodes(r.fun), color=:black, markersize=10)
shg()
```

Above, the black markers show the nodes of the interpolant. We can view the same approximation within the right half-plane by composing $r$ with $\phi^{-1}$:

```@example circle
φ⁻¹ = inv(φ)
domaincolor(r ∘ φ⁻¹, [-8, 8, -8, 8], abs=true)
lines!([(0, 8), (0, -8)], color=:white, linewidth=5)
scatter!(φ.(nodes(r.fun)), color=:black, markersize=10)
limits!(-8, 8, -8, 8)
shg()
```
