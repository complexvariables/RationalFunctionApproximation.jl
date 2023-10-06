---
title: 'RationalFunctionApproximation.jl: A Julia package for approximation by rational functions'
tags:
  - Julia
  - complex analysis
  - function approximation
  - rational approximation
authors:
  - name: Tobin A. Driscoll
    orcid: 0000-0002-1490-2545
    affiliation: "1" 
affiliations:
 - name: Department of Mathematical Sciences, University of Delaware, Newark, DE, USA
   index: 1
date: 5 Oct 2023
bibliography: paper.bib
---

# Summary

In scientific and engineering problems, it is often necessary or useful to replace an exact but lengthy, or even unknown, expression for a varying value by an efficient approximation. Classically, this is known as the problem of *function approximation* in mathematics. In function approximation, the ideal is to find a constructive procedure that achieves high accuracy with as succinct a representation as possible. One class of approximations that can do an outstanding job in this regard are the *rational functions*. A rational function is the ratio of two polynomials, which themselves are the combination of repeated multiplication and addition. For theoretical reasons, it is known that rational functions are excellent for representing many commonly encountered functions. Algorithmic developments within the last five years have opened a way to rapidly find rational approximations that are customized to a given problem. This package provides that capability to the Julia language.

# Statement of need

The `RationalFunctionApproximation.jl` package provides the automatic approximation of a function of a single variable by means of a rational function. Function approximation is a useful tool in scientific computation, as evidenced by the success of Chebfun in Matlab [@DriscollChebfunGuide2014] and ApproxFun in Julia [@OlverPracticalFramework2014], both of which are based on Chebyshev and trigonometric polynomial approximations.

Rational functions are a generalization of polynomials that are well-known to improve on polynomials' approximation power for functions that have singularities located near the approximation domain. The AAA algorithm [@NakatsukasaAAAAlgorithm2018a] opened the door to fast adaptive construction of rational interpolants of a function. A number of interesting applications rapidly followed, e. g., @BudisaRationalApproximation2022, @CostaAAAleastSquares2023, @DengExponentialAsymptotics2023a, @DereviankoESPRITESPIRA2023, @GopalSolvingLaplace2019a. The [`BaryRational`](https://github.com/macd/BaryRational.jl) package implements the original AAA algorithm in Julia.

More recently, a new variation of AAA [@DriscollAAARational2023] treats the domain of the given function as a continuum, rather than requiring it to be discretized. This can have significant advantages when singularities are very close to the domain, since exponential spacing of nodes is required.  `RationalFunctionApproximation.jl` implements this new algorithm in Julia, incorporating the [`ComplexRegions.jl`](https://github.com/complexvariables/ComplexRegions.jl) package[@DriscollComplexRegionsJl2019a] to provide general curves and regions as domains in the complex plane.

# Usage examples

We approximate the function 

$$
f(x) = \frac{0.01}{\sin(x^2 + 0.2x + 0.02)},
$$

which has an off-center hump over the interval $[-1,1]$, using the following:

```julia
julia> using RationalFunctionApproximation
julia> f = x -> 1/10^2 / sin(x^2 + 2x/10 + 2/10^2);
julia> r = approximate(f, unit_interval)
Barycentric rational function of type (10,10) on the domain: Path with 1 curve
```

The result is a rational function of type $(10,10)$, meaning that the numerator and denominator are polynomials of degree 10. The approximation interpolates the original function at 11 greedily selected points retrievable as `nodes(r)`, and it is accurate to over 13 digits over the domain, as shown in \autoref{fig:error_sin}.

![Error in a rational approximation.\label{fig:error_sin}](error_sin.pdf)

The approximation can be visualized in the complex plane using the [`DomainColoring`](https://eprovst.github.io/DomainColoring.jl/) package, which uses hue to indicate phase and brightness to show contours of log magnitude. The result is shown in \autoref{fig:complex}.

![Rational approximation in the complex plane.\label{fig:complex}](complex.pdf)

As one sees in the figure, there are poles in the approximation that are extremely close to the exact poles of $f$ at $x = -0.1 \pm 0.1i$:

```julia
julia> p = poles(r)
10-element Vector{ComplexF64}:
   -2.515988566699771 + 0.0im
  -1.8696464689420669 + 0.0im
  -0.1012939445069749 - 1.784487403808564im
 -0.10129394450697489 + 1.784487403808564im
 -0.09999999999999569 - 0.09999999999998077im
 -0.09999999999999569 + 0.09999999999998077im
 -0.08503139961597209 - 2.1919510617575773im
 -0.08503139961597207 + 2.1919510617575777im
   1.6697943872716796 + 0.0im
   2.3073258267815193 + 0.0im
```

We can also find the residues at those poles:

```julia
julia> residues(r, p[5:6])
 -4.686750511456871e-15 + 0.049999999999976126im
 -4.686750511456871e-15 - 0.049999999999976126im
```

As seen in \autoref{fig:error_sin}, the rational approximation is not quite uniformly accurate over the domain, because the AAA algorithm uses a least-squares criterion. However, there is an iteration due to Lawson [@LawsonContributionsTheory1961] that is employed by `minimax(r)` to more nearly recover the rational approximation that is optimal in the infinity norm. The resulting error, shown in \autoref{fig:error_minimax}, is roughly equioscillatory.

![Error in a nearly-minimax rational approximation.\label{fig:error_minimax}](error_minimax.pdf)

We can see the usefulness of the improved, continuous AAA algorithm [@DriscollAAARational2023] by approximating a function that has a singularity very close to its domain, such as

$$
g(x) = \tanh \bigl(5000(x-0.2)\bigr).
$$

If the proximity of the pole is comparable to or less than the sample spacing, the approximation by the original, discrete AAA method (which is implemented in the package) will be poor:

```julia
julia> g = x -> tanh(5000 * (x - 0.2));
julia> x = -1:0.001:1;
julia> r = aaa(x, g.(x));
julia> xx = -1:.0001:1;
julia> println("error = $(maximum(abs, @. g(xx) - r(xx) ))");
error = 0.1401166557529342
```

The error could be much reduced by using a finer grid of samples. But the continuous AAA algorithm finds those samples automatically:

```julia
julia> r = approximate(g, unit_interval);
julia> println("error = $(maximum(abs, @. g(xx) - r(xx) ))");
error = 4.0453862482081604e-11
```

The domain of an approximation can be any region defined using the `ComplexRegions` package. Here, we approximate the function $\tan(1/z^4)$ by a rational function that is analytic in the exterior of the unit circle, resulting in \autoref{fig:tan}.

```julia
r = approximate(z -> tanh(1/z^4), exterior(unit_circle));
```

![Approximation of $\tan(1/z^4)$.\label{fig:tan}](tan.pdf)

Finally, here is the approximation of a function with a branch cut over the interior of a cross-shaped region:

```julia
julia> import ComplexRegions.Shapes, ComplexPlots
julia> r = approximate(z -> log(0.35 + 0.4im - z), interior(Shapes.cross))
Barycentric rational function of type (5,5) on the domain: Region interior to Polygon with 12 sides
```

The result is shown in \autoref{fig:cross}.

![Approximation of a log function.\label{fig:cross}](cross.pdf)

# Acknowledgements

I grartefully acknowledge the help and support of my SISC paper coauthors, Yuji Nakatsukasa and Nick Trefethen. I also thank Daan Huybrechs for a valuable Julia-specific suggestion.

# References
