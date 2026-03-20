module ComputationalPhysics

using LinearAlgebra
using LaTeXStrings
using Plots
using Distributions
using Measures
using SparseArrays
using ForwardDiff

# Module for consistent plotting
include("visualizations.jl")
export plot_init, plot_generic, scatter_generic, plot_add!, scatter_add!, multi_plot
export plot_comparison, plot_convergence, save_plot, save_gif

# Main modules 
include("error_analysis.jl")
export var_double_pass32, var_double_pass64
export var_single_pass32, var_single_pass64
export mclaurin_series, kappa_f

include("linear_systems.jl")
export forward_substitution, forward_substitution_multiple, backward_substitution, backward_substitution_multiple
export lu_decomposition, lu_decomposition_pivoting
export plu_factorization, determinant_plu
export matrix_condition_number, solve_least_squares
export qr_mgs, test_qr, qr_mgs_augmented, pure_qr_algorithm, matrix_infinity_norm

include("interpolations.jl")
export barycentric_lagrange, chi_square, p_value, fit_goodness

include("nonlinear_equatioins.jl")
export bisection, bisection_plot, convergence, plot_convergence_analysis, newton_method, secant_method, inverse_quadratic_interpolation, plot_global_convergence

include("numerical_integration.jl")
export composite_trapezoidal, composite_simpson, gauss_legendre_quadrature, fajer_rule, clenshaw_curtis_rule, double_exponential_quadrature

include("ode.jl")
export euler_method, ie2, rk4_old, rk4

# Extras
include("schrodinger_equation.jl")
export t_independent_hamiltonian, schrodinger_solver1D, plot_snapshots, compute_RT

# include("fft.jl")
include("helpers.jl")
export make_log_safe

end # module ComputationalPhysics
