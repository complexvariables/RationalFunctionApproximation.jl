# RationalFunctionApproximation.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/dev/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15446359.svg)](https://doi.org/10.5281/zenodo.15446359)
[![codecov](https://codecov.io/gh/complexvariables/RationalFunctionApproximation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/complexvariables/RationalFunctionApproximation.jl)

![logo](logo-sm.png)

This Julia package adaptively computes rational approximations (i.e., ratios of polynomials) for functions on intervals and other domains in the complex plane.

The [documentation](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/) includes a walkthrough showing off the main capabilities of the package.

## Related work

* The [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/) package provides rational functions, but not in a way related to function approximation.
* The [BaryRational](https://juliahub.com/ui/Packages/General/BaryRational) package implements the original (fully discrete) version of the AAA algorithm, as well as Floater–Hormann rational interpolation.
* The [ApproxFun](https://juliaapproximation.github.io/ApproxFun.jl/stable) package provides 1D and multidimensional function approximation using Chebyshev polynomials and Fourier series. It also has extensive functionality for manipulating the approximations and for solving differential equations.
* There is an [ApproxFunRational](https://github.com/tomtrogdon/ApproxFunRational.jl) package, but it is undocumented.