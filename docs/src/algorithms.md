# Algorithms

```@example convergence
using RationalFunctionApproximation, CairoMakie
```

There are three algorithms for rational approximation offered in this version of the package. 

## AAA

The most robust method is the [continuum variant of the AAA algorithm](https://doi.org/10.1137/23M1570508) (see also the [arXiv version](https://arxiv.org/abs/2305.03677)). Briefly speaking, the AAA algorithm maintains a set of interpolation nodes and a set of test points on the boundary of the domain. The node set is grown iteratively by finding the best current approximation in a least-squares sense on the test nodes, then greedily adding the worst-case test point to the node set.

The `convergenceplot` function shows the errors of the approximants found during the AAA iteration.

```@example convergence
f = x -> cos(11x)
r = approximate(f, unit_interval)
convergenceplot(r)
```

(The plots in this documentation are made using `CairoMakie`, but the same functions are made available for `Plots`.)  In the plot above, the markers show the estimated max-norm error of the AAA rational interpolant over the domain as a function of the iteration counter. Each iteration adds one degree to both the numerator and denominator of the rational approximation. The gold halo indicates the final approximation chosen by the algorithm. The red dots indicate that a pole of the rational interpolant lies in the approximation domain. We can verify this fact by using `rewind` to recover the approximation from iteration 7:

```@repl convergence
r7 = rewind(r, 7)
poles(r7)
```

In most cases, the poles move out of the approximation domain as the iteration proceeds to better accuracy. However, it is possible for the iteration to stagnate if the original function has a singularity very close to the domain.

```@example convergence
f = x -> tanh(3000*(x - 1//4))
r = approximate(f, unit_interval)
convergenceplot(r)
```

This effect is thought to be mainly due to roundoff and conditioning of the problem. In this case, if we use more accurate floating-point arithmetic, we can see that the AAA convergence continues steadily past the previous plateau. In the following, we apply `Double64` arithmetic, having used exact rational numbers already in the definition of `f`:

```@example convergence
using DoubleFloats, ComplexRegions
r = approximate(f, Segment{Double64}(-1, 1))
convergenceplot(r)
```

In the extreme case of a function with a singularity on the domain, the convergence can be substantially affected:

```@example convergence
f = x -> abs(x - 1/8)
r = approximate(f, unit_interval)
convergenceplot(r)
```

In such a case, we might get improvement by increasing the number of allowed consecutive failures via the `stagnation` keyword argument:

```@example convergence
r = approximate(f, unit_interval, stagnation=50)
convergenceplot(r)
```

## Prescribed poles

While finding a rational interpolant is a nonlinear problem, we can compute a linear variant if we prescribe the poles of the rational function. Specifically, given the poles $\zeta_1,\ldots, \zeta_n$, we can find a polynomial $p$ and residues $w_k$ such that

```math
f(z) \approx p(z) + \sum_{k=1}^n \frac{w_k}{z - \zeta_k}. 
```

When posed on a discrete set of test points, this is a linear least-squares problem. In order to represent the polynomial stably, an Arnoldi iteration is used to find a well-conditioned basis for the test points. 

There is no iteration on the degree of the polynomial or rational parts of the approximant. Instead, the discretization of the boundary of the domain is refined iteratively until either the max-norm error is below a specified threshold or has stopped improving.

```@example convergence
f = x -> tanh(x)
ζ = 1im * π * [-1/2, 1/2, -3/2, 3/2]
r = approximate(f, Segment(-2, 2), ζ)
```

```@example convergence
max_err(r) = println("Max error: ", maximum(abs, check(r, quiet=true)[2]))
max_err(r);
```

To get greater accuracy, we can increase the degree of the polynomial part.

```@example convergence
r = approximate(f, Segment(-2, 2), ζ; degree=20)
max_err(r);
```

Note that the residues (in the exact function, all equal to one) may be accurate at the poles closest to the domain, but much less so elsewhere.

```@example convergence
Pair.(residues(r)...)
```

Suppose now we approximate $|x|$ using AAA. We can extract the poles of the result.

```@example convergence
r = approximate(abs, unit_interval, tol=1e-9)
ζ = poles(r)
```

To what extent might these poles be suitable for a different function that has the same singularity?

```@example convergence
s = approximate(x -> exp(abs(x)), unit_interval, ζ; degree=20)
max_err(r);
```

## Thiele continued fractions (experimental)

The AAA algorithm is remarkably robust, but a potential downside is that constructing an approximant of degree $n$ requires an SVD that takes $O(n^3)$ time. Thus, for an iteration up to degree $n$, the work is $O(n^4)$. In many cases, the value of $n$ is small enough not to cause a concern, but there is some motivation to find a more efficient algorithm.

One possibility is to use a continued fraction representation of the rational approximant. The Thiele representation, dating back to the 1800s, uses inverse divided differences to represent a rational interpolant, requiring just $O(n)$ time to add a node to an approximation of degree $n$. Divided differences are notoriously unstable, but [work by Salazar in 2024](https://doi.org/10.1007/s11253-024-02344-5) (or [the arXiv version](http://arxiv.org/abs/2109.10529)) indicates that a greedy adaptive approach can reduce or eliminate the instability.  

To try greedy Thiele, use `method=Thiele` as an argument to `approximate`. 

```@example convergence
f = x -> cos(41x - 5) * exp(-10x^2)
r = approximate(f, unit_interval; method=Thiele)
convergenceplot(r)
```

The $x$-axis of the convergence plot shows the degree of the denominator polynomial. Because the Thiele method alternates between interpolants of type $(n, n)$ and $(n+1, n)$, there are two dots above for each degree.

The primary appeal of the greedy Thiele method is that adding a node to an interpolant of degree $n$ takes $O(n)$ time, compared to $O(n^3)$ for the AAA method. Some experiments suggest that the Thiele method is faster in practice, though a systematic comparison is still needed.

```@example convergence
f = x -> abs(x-0.5)
@elapsed r = approximate(f, unit_interval; method=Barycentric, max_iter=150, allowed=true)
max_err(r);
```

```@example convergence
@elapsed r = approximate(f, unit_interval; method=Thiele, max_iter=300, allowed=true)
max_err(r);
```
