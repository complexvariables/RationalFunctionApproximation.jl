# Contributing to RationalFunctionApproximation.jl

Thanks for your interest in contributing! Contributions of all kinds are welcome: bug reports, documentation
fixes, new features, performance improvements, and tests.

## Ways to contribute

- **Report a bug** by opening an [issue](https://github.com/complexvariables/RationalFunctionApproximation.jl/issues).
  Please include a minimal, reproducible example, the output you got, the output
  you expected, and your Julia and package versions (`versioninfo()` and
  `] status RationalFunctionApproximation`).
- **Suggest a feature or enhancement** by opening an issue describing the use case.
- **Improve documentation**, including the docstrings, the walkthrough, and this file.
- **Submit code** via a pull request (see below).

For larger changes, please open an issue to discuss the design before investing
significant effort.

## Development setup

1. Fork and clone the repository.
2. Start Julia in the project and instantiate the environment:

   ```julia
   julia> ] activate .
   julia> ] instantiate
   ```

3. For interactive development, [Revise.jl](https://github.com/timholy/Revise.jl)
   is recommended so that edits are picked up without restarting Julia:

   ```julia
   julia> using Revise
   julia> using RationalFunctionApproximation
   ```

## Running the tests

The test suite uses [ReTest](https://github.com/JuliaTesting/ReTest.jl). To run
the full suite:

```julia
julia> ] test
```

or from the shell:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

ReTest also lets you run a subset of tests interactively. From the `test`
environment:

```julia
julia> include("test/RFATests.jl")
julia> RFATests.runtests()                 # run everything
julia> RFATests.runtests("real_interval")  # run tests matching a pattern
```

Please add or update tests for any change in behavior. New tests generally
belong in the relevant file under `test/` (e.g. `real_interval.jl`,
`circle.jl`, `parfrac.jl`) and should be wired into `test/RFATests.jl` if you
add a new file.

## Building the documentation

The documentation is built with [Documenter.jl](https://documenter.juliadocs.org/).
To build it locally:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site appears in `docs/build`.

## Pull request guidelines

- Branch from `main` and keep each PR focused on a single logical change.
- Match the existing code style: this package uses a two-parameter type system
  (`T` for float precision, `S` for value type) and supports generic arithmetic.
  See `CLAUDE.md` for an overview of the source layout and key types.
- Add docstrings to new public functions and types, and update the relevant
  documentation pages.
- Ensure the full test suite passes locally before opening the PR.
- Write a clear PR description explaining the motivation and approach. Reference
  any related issues.
- CI (tests, documentation, and benchmarks) runs automatically on PRs; please
  address any failures.

## Reporting security issues

If you discover a security-sensitive issue, please contact the maintainer
directly (see the email in `Project.toml`) rather than opening a public issue.

## License and attribution

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers this project.

If your work builds on or contributes to the methods implemented here, the
relevant references are listed in `README.md`, `CITATION.bib`, and `CLAUDE.md`.

## Code of conduct

Please be respectful and constructive in all interactions. We aim to maintain a
welcoming and collaborative community.
