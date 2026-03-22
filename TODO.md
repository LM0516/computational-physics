# TODO

## Correctness & Algorithm Quality
- Audit all algorithm implementations against reference textbooks for correctness
- Fix the typo in `nonlinear_equatioins.jl` (rename to `nonlinear_equations.jl` and update all references in `ComputationalPhysics.jl` and exercises)
- Remove `rk4_old` export from `ComputationalPhysics.jl` (dead code)

## Testing
- Populate `test/runtests.jl` with unit tests using Julia's `Test` stdlib
- Add tests for each module: `linear_systems`, `interpolations`, `nonlinear_equations`, `numerical_integration`, `ode`
- Add property-based / convergence tests (e.g. verify RK4 achieves O(h⁴) on a known ODE)
- Integrate test suite in CI (GitHub Actions with `julia --project=. test/runtests.jl`)

## Performance & Analysis
- Measure and document computational complexity for each major algorithm
- Benchmark algorithms with `BenchmarkTools.jl` and add a `benchmarks/` directory
- Profile memory allocations in the ODE solvers (avoid allocating inside loops)

## Package & Code Quality
- Add `FFT` module (currently `fft.jl` is commented out in `ComputationalPhysics.jl`)
