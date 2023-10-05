# RationalFunctionApproximation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://complexvariables.github.io/RationalFunctionApproximation.jl/dev/)

This package uses the continuous form of the AAA algorithm to adaptively compute rational approximations of functions on intervals and other domains in the complex plane.  For the mathematical description, see [AAA rational approximation on a continuum](https://arxiv.org/abs/2305.03677), which is to appear in SISC.

The [documentation](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/) includes a walkthrough showing off the main capabilities of the package.

Example: Approximation of $f(z) = (z^3 - 1) / \sin(z - 0.9 - i)$ to more than 13 digits on the unit circle:

![](https://complexvariables.github.io/RationalFunctionApproximation.jl/stable/index-099466ee.png)
