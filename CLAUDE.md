# RationalFunctionApproximation.jl — AI Agent Reference

## Purpose

Julia package for adaptively computing rational approximations (ratios of polynomials) to functions on real intervals, complex curves, and discrete point sets. Rational functions often achieve superior accuracy over polynomials for functions with singularities, poles, or steep variations.

Described briefly in preprint by Driscoll, doi:10.48550/arXiv.2512.06140

Implements the AAA algorithm and related variants. Key references:
- Nakatsukasa, Sète, Trefethen (2018) SISC — original AAA (doi:10.1137/16m1106122)
- Driscoll, Nakatsukasa, Trefethen (2024) SISC — AAA on a continuum (doi:10.1137/23M1570508)
- Brubeck, Nakatsukasa, Trefethen (2021) SIAM Review — Vandermonde with Arnoldi (doi:10.1137/19M130100X)

Implements Thiele continued fractions. Key references:
- Salazar Celis (2024), Ukr. Math J. (doi:10.1007/s11253-024-02344-5)
- Driscoll and Zhou (2025 preprint, doi:10.48550/arXiv.2510.07295)
- Driscoll and Zhou (2026 preprint, doi:10.48550/arXiv.2601.10667)

---

## Source File Map

```
src/
  RationalFunctionApproximation.jl   — Main module, exports
  abstract-rational.jl               — Abstract type hierarchy
  aaa.jl                             — Legacy AAA (discrete) implementation
  barycentric.jl                     — Barycentric representation + AAA continuum algorithm
  thiele.jl                          — Thiele continued fraction representation
  parfrac.jl                         — Partial fractions + ArnoldiBasis
  approximation.jl                   — approximate() function, ContinuumApproximation, DiscreteApproximation
  discretized_path.jl                — Path discretization for complex domains
  lawson.jl                          — Minimax approximation via Lawson iteration
  utils.jl                           — Domain helpers and refinement utilities

test/
  runtests.jl, RFATests.jl           — Runner and main test module (ReTest framework)
  aaa.jl, discrete.jl, real_interval.jl, imag_interval.jl
  circle.jl, discretized-path.jl, operations.jl, parfrac.jl, custom.jl

ext/
  RFAForwardDiffExt.jl               — ForwardDiff autodiff
  RFAMakieExt.jl, RFAPlotsExt.jl    — Visualization
  RFAPythonCallExt.jl                — Python interop
  RFAZygoteExt.jl                    — Zygote autodiff
```

---

## Key Types

### Abstract hierarchy
- `AbstractRationalFunction{S}` — base; `S` = value type
- `AbstractRationalInterpolant{T,S}` — adds `T` = float precision

### Concrete rational function types

**`Barycentric{T,S}`** (barycentric.jl:13) — default method; alias `AAA = Barycentric`  
Fields: `nodes`, `values`, `weights`, `w_times_f`

**`Thiele{T,S}`** (thiele.jl) — continued fraction representation; alias `TCF = Thiele`  
Fields: `nodes`, `values`, `weights`  
Two evaluation strategies: classic and numerically stable `onediv`.

**`PartialFractions{S}`** (parfrac.jl:158) — for approximation with prescribed poles  
Fields: `polynomial::ArnoldiPolynomial`, `poles`, `residues`

### Supporting types

- **`ArnoldiBasis{T}`** (parfrac.jl:11) — orthogonal polynomial basis; fields: `nodes`, `Q`, `H`
- **`ArnoldiPolynomial{T}`** (parfrac.jl:63) — polynomial in Arnoldi basis; fields: `coeff`, `basis`
- **`DiscretizedPath{T,F}`** (discretized_path.jl:9) — complex curve discretization with adaptive refinement
- **`ContinuumApproximation{T,S,R}`** (approximation.jl:37) — continuous domain wrapper; fields: `original`, `domain`, `fun`, `allowed`, `path`, `history`
- **`DiscreteApproximation{T,S,R}`** (approximation.jl:75) — discrete point set wrapper; fields: `data`, `domain`, `fun`, `test_index`, `allowed`, `history`
- **`IterationRecord{R,S,T}`** (approximation.jl:5) — convergence history entry; fields: `interpolant`, `error`, `poles`

---

## Public API

### Approximation construction
- `approximate(f, domain; method, max_iter, tol, allowed, refinement, stagnation)` — main entry point
- `approximate(f, domain, poles)` — least-squares with prescribed poles
- `aaa(y, z; kwargs...)` — legacy discrete AAA (deprecated)

### Rational function queries
- `nodes(r)`, `values(r)`, `weights(r)` — interpolation data
- `degrees(r)` → `(num_degree, den_degree)`, `degree(r)` → denominator degree
- `poles(r)`, `residues(r)`, `roots(r)` — singularity/zero data
- `Res(f, z)` — residue via contour integration

### Evaluation & derivatives
- `r(z)` or `evaluate(r, z)` — pointwise evaluation
- `derivative(r, order=1)` — compute derivatives
- `check(r; kwargs...)` — verify approximation accuracy on domain

### Approximation management
- `get_function(r)`, `domain(r)` — extract components
- `rewind(r, index)` — revert to earlier iteration
- `get_history(r)` — convergence history
- `test_points(r)` — test point locations

### Optimization
- `minimax(r, f, nsteps=20)` — Lawson iteration for minimax refinement

### Visualization (optional extensions)
- `convergenceplot(r)`, `errorplot(r)`, `poleplot(r)`

### Domain constants
- `unit_interval` = Segment(-1, 1)
- `unit_circle` = Circle(0, 1)
- `unit_disk` = disk(0, 1)

---

## Algorithms

### AAA 
- Barycentric representation; weights via Loewner matrix SVD
- Iteratively selects new nodes at maximum error locations
- Stagnation detection for early stopping
- Works on continuous domains and discrete point sets

#### Minimax via Lawson Iteration
- Iterative node reweighting to minimize maximum error
- Post-processing step applied to an existing Barycentric approximant

### Thiele Continued Fraction
- Alternative to barycentric; better conditioning for some function classes
- Poles via generalized eigenvalue problem (thiele.jl:180)
- Efficient derivatives via recurrence relations

### Partial Fractions with Prescribed Poles
- Linear least-squares with Arnoldi orthogonalization
- Used when singularity locations are known in advance

---

## Dependencies

### Required
- ComplexRegions (0.3) — complex curve/path abstractions
- ComplexValues (0.3) — complex number utilities
- GenericLinearAlgebra (0.3, 0.4), GenericSchur (0.5) — generic linear algebra
- PrecompileTools, PyFormattedStrings, Logging, Printf, Statistics, LinearAlgebra

### Optional (weak dependencies / extensions)
- ForwardDiff, Zygote — autodiff
- Makie, Plots — visualization
- PythonCall — Python interop

### Test only
- ReTest, DoubleFloats

---

## Design Decisions

1. **Two-parameter type system**: `T` for float precision, `S` for value type — supports generic arithmetic.
2. **Barycentric as default**: most efficient/stable; aliased `AAA` for historical compatibility.
3. **Continuum vs. Discrete split**: `ContinuumApproximation` and `DiscreteApproximation` reflect fundamentally different strategies.
4. **Adaptive path discretization**: `DiscretizedPath` stores multiple refinement levels in matrix form.
5. **`allowed` parameter**: generic function to filter pole locations; enables multiply-connected domains.
6. **Convergence history**: optional recording enables `rewind()` and convergence plots.
7. **Extension architecture**: plotting and autodiff are optional — no hard dependencies.
8. **Precompilation workload**: uses `@compile_workload` for fast time-to-first-approximation.
