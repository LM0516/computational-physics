# Computational Physics

A comprehensive collection of numerical methods and exercises written in **Julia**, developed as part of the Computational Physics course at **Università di Milano Bicocca (UniMiB)**. All algorithms are implemented from scratch, accompanied by detailed worked exercises and a comprehensive LaTeX report.

## Topics

| # | Topic | Key Algorithms |
|---|-------|----------------|
| 1 | **Error Analysis** | Machine epsilon, roundoff errors, stable quadratic formula |
| 2 | **Linear Systems** | Forward/backward substitution, LU decomposition (with/without pivoting), QR (MGS), condition numbers, least squares |
| 3 | **Interpolation** | Barycentric Lagrange interpolation, Chebyshev nodes, χ² fit quality |
| 4 | **Roots of Nonlinear Equations** | Bisection, Newton, secant, inverse quadratic interpolation, convergence analysis |
| 5 | **Numerical Integration** | Composite trapezoidal & Simpson, Gauss–Legendre, Fejér, Clenshaw–Curtis, double-exponential |
| 6 | **Ordinary Differential Equations** | Euler, improved Euler (IE2), Runge–Kutta (RK4), Schrödinger equation solver, three-body problem |

## Repository Structure

```
computational-physics/
├── src/                              # Julia package source (ComputationalPhysics module)
│   ├── ComputationalPhysics.jl       # Main module — re-exports all submodules
│   ├── error_analysis.jl             # Variance estimators, Maclaurin series, condition number
│   ├── linear_systems.jl             # LU/PLU, QR, substitution, norms, least squares
│   ├── interpolations.jl             # Barycentric Lagrange, χ², p-value, fit quality
│   ├── nonlinear_equatioins.jl       # Bisection, Newton, secant, IQI, convergence tools
│   ├── numerical_integration.jl      # Trapezoidal, Simpson, GL, Fejér, CC, DE quadrature
│   ├── ode.jl                        # Euler, IE2, RK4 for scalar and vector-valued ODEs
│   ├── schrodinger_equation.jl       # Sparse Hamiltonian, 1D solver, reflection/transmission
│   ├── visualizations.jl             # Consistent plot style and save helpers
│   └── helpers.jl                    # Utility functions (e.g., make_log_safe)
├── exercises/                        # One subdirectory per course topic
│   ├── 01-error-analysis/            # Machine epsilon, roundoff errors, condition numbers
│   ├── 02-linear-systems/            # Direct and iterative solvers, conditioning
│   ├── 03-interpolations/            # Lagrange, least-squares, error analysis
│   ├── 04-roots-of-nonlinear-equations/  # Root-finding algorithms, convergence
│   ├── 05-numerical-integration/     # Quadrature methods, adaptive integration
│   └── 06-ordinary-differential-equations/  # ODE solvers, Schrödinger equation
├── notebooks/                        # Jupyter/Pluto notebooks for exploration
│   ├── 3body_problem.ipynb           # Three-body problem simulation
│   ├── julia_test.ipynb
│   └── ...
├── test/                             # Test suite
│   ├── runtests.jl                   # Standard Julia test entry point
│   ├── submodules/                   # Individual module tests
│   └── test_all.sh                   # Convenience test runner
├── report/                           # LaTeX report with full write-ups
│   ├── main.tex                      # Main report document
│   ├── preamble.sty                  # Custom LaTeX style
│   ├── chapters/                     # Per-topic LaTeX chapters
│   └── figures/                      # Generated plots and diagrams
├── output/                           # Generated figures and GIFs
├── Project.toml                      # Julia package manifest
└── Manifest.toml
```

## Modules

All implementations live in `src/` and are exposed through the `ComputationalPhysics` package. Include the package in any exercise with `using ComputationalPhysics`.

| Module file | Exported functions |
|---|---|
| `error_analysis.jl` | `var_double_pass32/64`, `var_single_pass32/64`, `mclaurin_series`, `kappa_f` |
| `linear_systems.jl` | `forward_substitution`, `backward_substitution`, `lu_decomposition`, `lu_decomposition_pivoting`, `plu_factorization`, `determinant_plu`, `matrix_condition_number`, `solve_least_squares`, `qr_mgs`, `pure_qr_algorithm`, `matrix_infinity_norm` |
| `interpolations.jl` | `barycentric_lagrange`, `chi_square`, `p_value`, `fit_goodness` |
| `nonlinear_equations.jl` | `bisection`, `newton_method`, `secant_method`, `inverse_quadratic_interpolation`, `convergence`, `plot_convergence_analysis`, `plot_global_convergence` |
| `numerical_integration.jl` | `composite_trapezoidal`, `composite_simpson`, `gauss_legendre_quadrature`, `fajer_rule`, `clenshaw_curtis_rule`, `double_exponential_quadrature` |
| `ode.jl` | `euler_method`, `ie2`, `rk4` |
| `schrodinger_equation.jl` | `t_independent_hamiltonian`, `schrodinger_solver1D`, `plot_snapshots`, `compute_RT` |
| `visualizations.jl` | `plot_init`, `plot_generic`, `scatter_generic`, `plot_add!`, `scatter_add!`, `multi_plot`, `plot_comparison`, `plot_convergence`, `save_plot`, `save_gif` |
| `helpers.jl` | `make_log_safe` and other utility functions |

## Getting Started

### Prerequisites

- [Julia](https://julialang.org/downloads/) ≥ 1.10
- Dependencies are managed by `Project.toml` — no manual installation needed.

### Install dependencies

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

### Running an Exercise

Each exercise is a standalone Julia script that loads the package via a relative path. To run one:

```bash
julia --project=. exercises/01-error-analysis/1-2-representation-of-arithmetic-and-roundoff-errors/es1.jl
```

Exercise subdirectories typically contain multiple exercises (e.g., `es1.jl`, `es2.jl`) related to the same topic.

### Running the Test Suite

```bash
julia --project=. test/runtests.jl
```

or using the shell script:

```bash
bash test/test_all.sh
```

### Building the Report

The comprehensive LaTeX report can be built using `latexmk`:

```bash
cd report && latexmk -pdf main.tex
```

The report includes full derivations, numerical results, and all generated figures organized by topic.

## Author

**Lorenzo Minuz** — [GitHub](https://github.com/LM0516)

## License

This project is provided for educational purposes under the MIT License.
