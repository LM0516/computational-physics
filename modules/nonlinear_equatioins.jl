module NonlinearEquations
include("../modules/linear_systemsV2.jl")
using LinearAlgebra
using LaTeXStrings
using Plots
using .LinearSystemsV2
export bisection, bisection_plot, convergence, plot_convergence_analysis, newton_method, secant_method, inverse_quadratic_interpolation

"""
    bisection(f::Function, a, b; tol=1e-14)

Finds a root of the function `f` within the interval `[a, b]` using the bisection method.

# Arguments
- `f::Function`: The function for which to find a root.
- `a`: The lower bound of the interval.
- `b`: The upper bound of the interval.
- `tol::Float64`: (Optional) The desired tolerance for the interval width. Defaults to `1e-14`.

# Returns
- `Float64`: The approximate root of `f`.
"""
function bisection(f::Function, a, b; tol=1e-14)
    while b - a > tol
        mk = 1 / 2 * (a + b)
        prod = f(mk) * f(a)

        # Check the product of the mean values
        a = prod > 0 ? mk : a
        b = prod < 0 ? mk : b
    end
    sol = 1 / 2 * (a + b)
    return sol
end

"""
    bisection_plot(f::Function, a::Float64, b::Float64, f_equation::String)

Finds a root of the function `f` within the interval `[a, b]` using the bisection method and plots the convergence.

# Arguments
- `f::Function`: The function for which to find a root.
- `a::Float64`: The lower bound of the interval.
- `b::Float64`: The upper bound of the interval.
- `f_equation::String`: A string representation of the function `f` for plotting purposes (e.g., "x^2 - 2").

# Returns
- `Tuple{Vector{Float64}, Float64}`: A tuple containing:
    - `x_history::Vector{Float64}`: A vector of approximate roots at each iteration.
    - `sol::Float64`: The final approximate root.
"""
function bisection_plot(f::Function, a::Float64, b::Float64, f_equation::String)
    count = 0
    tol = 1e-14
    epsilon = 1.0

    # Crate an empty array for the approximation history
    x_history = Float64[]

    # Plot the function
    p = plot(f, a, b, label=latexstring(f_equation), xlabel=L"x", ylabel=L"y", legend=:topleft)

    # Plot the zeros
    hline!([0.0], label="")

    # Execute the bisection method
    while epsilon > tol
        # Show the point for every iteration
        scatter!(p, [a, b], [f(a), f(b)], label="", ms=2)

        mk = 1 / 2 * (a + b)
        prod = f(mk) * f(a)

        # Save the current approximation
        push!(x_history, mk)

        # Check the product of the mean values
        a = prod > 0 ? mk : a
        b = prod < 0 ? mk : b

        count += 1
        epsilon = b - a
    end
    sol = 1 / 2 * (a + b)

    scatter!(p, [sol], [f(sol)], label="Exact solution: x = $sol with $count iterations", ms=4)
    println()
    display(p)
    readline()

    return x_history, sol
end

"""
    convergence(x::Vector{Float64}, r::Float64; skip::Int=10)

Analyzes the order of convergence and asymptotic error constant for a sequence of approximations.

# Arguments
- `x::Vector{Float64}`: A vector of approximations obtained from an iterative method.
- `r::Float64`: The true root (exact solution).
- `skip::Int`: (Optional) Number of initial iterations to skip for the analysis. Defaults to `10`.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing:
    - `q::Float64`: The order of convergence.
    - `C::Float64`: The asymptotic error constant.
"""
function convergence(x, r; skip::Int=10)
    # Number of initial iterations to skip for the analysis
    # Use only the tail of the data
    real_skip = (length(x) > skip + 2) ? skip : 1
    x_tail = x[real_skip:end]

    epsilon = @. abs(x_tail - r)
    # Filter out exact zeros
    epsilon = epsilon[epsilon.>0]
    log_eps = log10.(epsilon)
    n = length(log_eps) - 1

    log_eps_n = log_eps[1:end-1] # ε[n1]
    log_eps_n1 = log_eps[2:end]  # ε[n+1]

    M = [ones(n) log_eps_n] # Design matrix
    b = log_eps_n1 # Target vector

    # Least Squares method for interpolation
    log_C, q = solve_least_squares(M, b)
    C = 10^log_C

    println("Order of convergence q = $q")
    println("Constant C = $C")

    return q, C
end

"""
    plot_convergence_analysis(x::Vector{Float64}, r::Float64; skip_initial::Int=0)

Generates a plot to visually analyze the convergence of an iterative method.
It plots `log10(|ε_{n+1}|)` against `log10(|ε_n|)` and fits a line to determine
the order of convergence and the asymptotic error constant.

# Arguments
- `x::Vector{Float64}`: A vector of approximations obtained from an iterative method.
- `r::Float64`: The true root (exact solution).
- `skip_initial::Int`: (Optional) Number of initial iterations to skip for the plot and analysis. Defaults to `0`.

# Returns
- `Plots.Plot`: The generated convergence analysis plot.
"""
function plot_convergence_analysis(x, r; skip_initial::Int=0)
    # Use only the tail of the data
    x_tail = x[skip_initial+1:end]

    epsilon = @. abs(x_tail - r)

    log_eps = log10.(epsilon)

    log_eps_n = log_eps[1:end-1]
    log_eps_n1 = log_eps[2:end]

    n = length(log_eps_n)
    M = [ones(n) log_eps_n]
    b = log_eps_n1

    # solve_least_squares 
    log10_C, q = solve_least_squares(M, b)
    C = 10^(log10_C)

    println("Analysis based on the last ", length(x_tail), " iterations.")
    println("Order of convergence q ≈ $q")
    println("Asymptotic error constant C ≈ $C")

    conv_plot = scatter(
        log_eps_n,
        log_eps_n1,
        label=L"Data ($\log_{10}|\epsilon_{n+1}|) vs (\log_{10}|\epsilon_n|))$",
        xlabel=L"\log_{10}|\epsilon_n|",
        ylabel=L"\log_{10}|\epsilon_{n+1}|",
        title="Convergence Analysis (Skipping first $skip_initial iterations)",
        legend=:topleft,
        markershape=:circle,
        markerstrokewidth=0
    )

    # Calculate the y-values for the fitted line
    y_fit = @. log10_C + q * log_eps_n

    plot!(
        conv_plot,
        log_eps_n,
        y_fit,
        label="Least-Squares Fit (slope q ≈ $(round(q, digits=3)))",
        color=:red,
        linewidth=2
    )

    display(conv_plot)
    readline()
    return conv_plot
end

"""
    newton_method(f::Function, df::Function, a::Real, b::Real; x_init::Real=b, xtol::Float64=1e-14, ftol::Float64=1e-14, max_iter::Int=100, bracketing::Bool=true)

Finds a root of the function `f` using Newton's method.

# Arguments
- `f::Function`: The function for which to find a root.
- `df::Function`: The derivative of the function `f`.
- `a::Real`: The lower bound of the bracketing interval (used if `bracketing` is true).
- `b::Real`: The upper bound of the bracketing interval (used if `bracketing` is true).
- `x_init::Real`: (Optional) The initial guess for the root. Defaults to `b`.
- `xtol::Float64`: (Optional) The desired tolerance for the change in `x`. Defaults to `1e-14`.
- `ftol::Float64`: (Optional) The desired tolerance for the function value at the root. Defaults to `1e-14`.
- `max_iter::Int`: (Optional) The maximum number of iterations. Defaults to `100`.
- `bracketing::Bool`: (Optional) If `true`, falls back to bisection if the Newton's method diverges outside `[a, b]`. Defaults to `true`.

# Returns
- `Tuple{Float64, Vector{Float64}}`: A tuple containing:
    - `Float64`: The approximate root of `f`.
    - `Vector{Float64}`: A vector of approximations generated during the iterations.
"""
function newton_method(f::Function, df::Function, a::Real, b::Real; x_init::Real=b, xtol::Float64=1e-13, ftol::Float64=1e-13, max_iter::Int=100, bracketing::Bool=true)
    # Initial 2 x values
    x1 = float(x_init == b ? b : x_init)
    x2 = x1 - f(x1) / df(x1)

    iter = 0
    x = zeros(max_iter)
    x[1], x[2] = x1, x2

    while ((abs(x1 - x2) >= xtol) || (abs(f(x1)) >= ftol)) && iter < max_iter
        x_new = x2 - f(x2) / df(x2)

        # Save the x values for the convergence analysis
        x[iter+1] = x_new

        x1 = x2
        x2 = x_new
        iter += 1
    end

    if iter == max_iter
        @warn "Maximum number of iterations reached for the root $x2"
        @warn "The function in this root is $(f(x2))"
    end

    # Bracketing 
    if x2 > b && f(a) * f(b) < 0 && bracketing == true
        println("Using bracketing...")
        x2 = bisection(f, a, b)
    end

    return x2, x
end

"""
    secant_method(f::Function, x_0::Real, x_1::Real; xtol::Float64=1e-14, ftol::Float64=1e-14, max_iter::Int=100, bracketing::Bool=true)

Finds a root of the function `f` using the secant method.

# Arguments
- `f::Function`: The function for which to find a root.
- `x_0::Real`: The first initial guess.
- `x_1::Real`: The second initial guess.
- `xtol::Float64`: (Optional) The desired tolerance for the change in `x`. Defaults to `1e-14`.
- `ftol::Float64`: (Optional) The desired tolerance for the function value at the root. Defaults to `1e-14`.
- `max_iter::Int`: (Optional) The maximum number of iterations. Defaults to `100`.
- `bracketing::Bool`: (Optional) If `true`, falls back to bisection if the secant method diverges outside the initial interval. Defaults to `true`.

# Returns
- `Tuple{Float64, Vector{Float64}}`: A tuple containing:
    - `Float64`: The approximate root of `f`.
    - `Vector{Float64}`: A vector of approximations generated during the iterations.
"""
function secant_method(f::Function, x_0::Real, x_1::Real; xtol::Float64=1e-14, ftol::Float64=1e-14, max_iter::Int=100, bracketing::Bool=true)
    x_prev = float(x_0)
    x_curr = float(x_1)
    iter = 0
    x = [x_prev, x_curr]

    while ((abs(x_curr - x_prev) >= xtol) || (abs(f(x_curr)) >= ftol)) && iter < max_iter
        f_x_curr = f(x_curr)
        f_x_prev = f(x_prev)

        if f_x_curr == f_x_prev
            # Avoid division by zero, return current approximation
            break
        end

        x_next = x_curr - (f_x_curr * (x_curr - x_prev)) / (f_x_curr - f_x_prev)

        # Save the x values for the convergence analysis
        push!(x, x_next)

        x_prev = x_curr
        x_curr = x_next
        iter += 1
    end

    # Bracketing (assuming x_0 and x_1 are the original interval bounds)
    if (x_curr > x_1 || x_curr < x_0) && f(x_0) * f(x_1) < 0 && bracketing == true
        println("Using bracketing...")
        x_curr = bisection(f, x_0, x_1)
    end

    return x_curr, x
end

"""
    inverse_quadratic_interpolation(f::Function, x0::Real, x1::Real, x2::Real; xtol::Float64=1e-14, ftol::Float64=1e-14, max_iter::Int=100)

Finds a root of the function `f` using the inverse quadratic interpolation method.
This method requires three initial guesses.

# Arguments
- `f::Function`: The function for which to find a root.
- `x0::Real`: The first initial guess.
- `x1::Real`: The second initial guess.
- `x2::Real`: The third initial guess.
- `xtol::Float64`: (Optional) The desired tolerance for the change in `x` or interval width. Defaults to `1e-14`.
- `ftol::Float64`: (Optional) The desired tolerance for the function value at the root. Defaults to `1e-14`.
- `max_iter::Int`: (Optional) The maximum number of iterations. Defaults to `100`.

# Returns
- `Tuple{Float64, Vector{Float64}}`: A tuple containing:
    - `Float64`: The approximate root of `f`.
    - `Vector{Float64}`: A vector of approximations generated during the iterations.
"""
function inverse_quadratic_interpolation(f::Function, x0::Real, x1::Real, x2::Real; xtol::Float64=1e-14, ftol::Float64=1e-14, max_iter::Int=100)
    a, b, c = float(x0), float(x1), float(x2)
    iter = 0
    x_history = [a, b, c]

    while iter < max_iter
        fa, fb, fc = f(a), f(b), f(c)

        if abs(fc) < ftol || abs(c - b) < xtol
            return c, x_history
        end

        # Lagrange interpolation formula for inverse function x(y) at y=0
        L0 = (fb * fc) / ((fa - fb) * (fa - fc))
        L1 = (fa * fc) / ((fb - fa) * (fb - fc))
        L2 = (fa * fb) / ((fc - fa) * (fc - fb))

        x_new = L0 * a + L1 * b + L2 * c

        push!(x_history, x_new)

        # Update points: shift old ones out
        a = b
        b = c
        c = x_new

        iter += 1
    end

    return c, x_history
end

end
