# RationalFunctionApproximation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/dev/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15446359.svg)](https://doi.org/10.5281/zenodo.15446359)
[![codecov](https://codecov.io/gh/complexvariables/RationalFunctionApproximation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/complexvariables/RationalFunctionApproximation.jl)

This package uses the continuous form of the AAA algorithm to adaptively compute rational approximations of functions on intervals and other domains in the complex plane.  For the mathematical description, see [AAA rational approximation on a continuum](https://arxiv.org/abs/2305.03677), which is to appear in SISC.

The [documentation](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/) includes a walkthrough showing off the main capabilities of the package.

Example: Approximation of $f(z) = (z^3 - 1) / \sin(z - 0.9 - i)$ to more than 13 digits on the unit circle:

![](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/index-099466ee.png)

## Related work

* The [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/) package provides rational functions, but not in a way related to function approximation.
* The [BaryRational](https://juliahub.com/ui/Packages/General/BaryRational) package implements the original (fully discrete) version of the AAA algorithm, as well as Floater–Hormann rational interpolation.
* The [ApproxFun](https://juliaapproximation.github.io/ApproxFun.jl/stable) package provides 1D and multidimensional function approximation using Chebyshev polynomials and Fourier series. It also has extensive functionality for manipulating the approximations and for solving differential equations.
* There is an [ApproxFunRational](https://github.com/tomtrogdon/ApproxFunRational.jl) package, but it is undocumented.