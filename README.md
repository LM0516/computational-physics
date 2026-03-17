# Computational Physics

A collection of numerical methods and exercises written in **Julia**, developed as part of a Computational Physics course. Each topic includes from-scratch implementations of core algorithms, worked exercises, and a companion LaTeX report.

## Topics

| # | Topic | Key Algorithms |
|---|-------|---------------|
| 1 | **Error Analysis** | Machine epsilon, roundoff errors, quadratic equation stability |
| 2 | **Linear Systems** | Forward/backward substitution, LU decomposition (with/without pivoting), QR decomposition, condition numbers, least squares |
| 3 | **Interpolations** | Polynomial interpolation, Lagrange–Waring (barycentric formula), Chebyshev nodes, interpolation error analysis |
| 4 | **Roots of Nonlinear Equations** | Bisection method, Newton's method, secant method, inverse quadratic interpolation, convergence analysis |
| 5 | **Numerical Integration** | Composite trapezoidal & Simpson rules, Gauss–Legendre quadrature, Clenshaw–Curtis rule, double-exponential quadrature |
| 6 | **Ordinary Differential Equations** | Euler method, improved Euler (IE2), Runge–Kutta (RK4), Schrödinger equation solver, three-body problem |

## Repository Structure

```
├── 01 - Error Analysis/          # Exercises on floating-point arithmetic
├── 02 - Linear Systems/          # LU, QR, least squares exercises
├── 03 - Interpolations/          # Lagrange interpolation exercises
├── 04 - Roots of Nonlinear Equations/  # Root-finding exercises
├── 05 - Numerical Integration/   # Quadrature exercises
├── 06 - Ordinary Differential Equations/  # ODE solver exercises
├── modules/                      # Reusable Julia modules
│   ├── linear_systems.jl         # Linear algebra routines
│   ├── interpolations.jl         # Barycentric Lagrange interpolation
│   ├── nonlinear_equatioins.jl   # Root-finding methods & convergence tools
│   ├── numerical_integration.jl  # Quadrature rules
│   ├── ode.jl                    # ODE solvers (Euler, IE2, RK4)
│   ├── schrodinger_equation.jl   # 1D Schrödinger equation solver
│   ├── ConsistentPlots.jl        # Uniform plotting style across all graphs
│   └── ...
├── report/                       # LaTeX report with full solutions
│   ├── main.tex
│   └── chapters/
├── tests/                        # Test suite & notebooks
│   ├── test_suite.jl
│   ├── test_all.sh
│   └── 3body_problem.ipynb
└── gif/                          # Animated visualizations
```

## Modules

All reusable implementations live in the `modules/` directory and are structured as Julia modules that can be included in any exercise:

- **`LinearSystems`** — Forward/backward substitution, LU decomposition (no pivoting, row pivoting), PLU factorization, determinant via PLU, matrix norms, condition number, least squares.
- **`Interpolations`** — Barycentric Lagrange interpolation with support for equispaced and Chebyshev nodes.
- **`NonlinearEquations`** — Bisection, Newton, secant, and inverse quadratic interpolation methods. Includes convergence order analysis and plotting utilities.
- **`NumericalIntegration`** — Composite trapezoidal & Simpson rules, Gauss–Legendre quadrature, Fejér rule, Clenshaw–Curtis rule, double-exponential quadrature.
- **`ODE`** — Euler method, improved Euler (IE2), and 4th-order Runge–Kutta (RK4) for scalar and vector-valued ODEs.
- **`SchrodingerEquation`** — Time-independent Hamiltonian construction (sparse matrices), 1D Schrödinger solver via RK4, reflection/transmission coefficient computation, wavefunction visualization.
- **`ConsistentPlots`** — Shared plotting defaults and helper functions for uniform graph style across all exercises.

## Getting Started

### Prerequisites

- [Julia](https://julialang.org/downloads/) (≥ 1.6)
- Julia packages: `LinearAlgebra`, `Plots`, `LaTeXStrings`, `SparseArrays`, `Measures`, `Distributions`

### Running an Exercise

Each exercise is a standalone Julia script. To run one:

```bash
julia "01 - Error Analysis/1.2 - Representation of Arithmetic and Roundoff Errors/es1.jl"
```

### Running the Test Suite

```bash
julia tests/test_all.sh
```

### Building the Report

The LaTeX report is located in the `report/` directory. Build it with:

```bash
cd report
latexmk -pdf main.tex
```

## TODO

- Implement a script to benchmark all the exercise (Use Benchmarks)
- Implement scripts to test the complexity of every algorithm
- Write an introduction for every section of every major chapter
- Check if it is worth using Make to compile LaTeX

## Author

**Lorenzo Minuz**

## License

This project is provided for educational purposes.
