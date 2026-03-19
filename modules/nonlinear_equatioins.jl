module NonlinearEquations
include("../modules/linear_systemsV2.jl")
include("../modules/ConsistentPlots.jl")
using LinearAlgebra
using LaTeXStrings
using Plots
using .LinearSystemsV2
using .ConsistentPlots
export bisection, bisection_plot, convergence, plot_convergence_analysis, newton_method, secant_method, inverse_quadratic_interpolation, plot_global_convergence

"""
Finds a root of the function `f` within the interval `[a, b]` using the bisection method.
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
Finds a root of the function `f` within the interval `[a, b]` using the bisection method and plots the convergence.
"""
function bisection_plot(f::Function, a::Float64, b::Float64, f_equation::String; save_dir::String="test")
    initialize_style()
    count = 0
    tol = 1e-14
    epsilon = 1.0

    # Crate an empty array for the approximation history
    x_history = Float64[]

    # Plot the function
    p = plot(f, a, b, label=latexstring(f_equation), xlabel=L"x", ylabel=L"y")

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

    label = "Exact solution: x = $sol with $count iterations"

    scatter!(p, [sol], [f(sol)], label=label, ms=4)

    save_graph(p, label, save_dir)
    display(p)
    readline()

    return x_history, sol
end

"""
Analyzes the order of convergence and asymptotic error constant for a sequence of approximations.
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
Generates a plot to visually analyze the convergence of an iterative method.
It plots `log10(|ε_{n+1}|)` against `log10(|ε_n|)` and fits a line to determine
the order of convergence and the asymptotic error constant.
"""
function plot_convergence_analysis(x, r; skip_initial::Int=0, save_dir::String="test")
    initialize_style()
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
    )

    save_graph(conv_plot, "conv-analysis-with-$skip_initial-skipped", save_dir)
    display(conv_plot)
    readline()
    return conv_plot
end

"""
Finds a root of the function `f` using Newton's method.
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
Finds a root of the function `f` using the secant method.
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
Finds a root of the function `f` using the inverse quadratic interpolation method.
This method requires three initial guesses.
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

"""
Plots a_n = -log10|x_n - r| as a function of n, directly addressing the exercise hint.
It calculates the global asymptotic error constant C from the slope.
"""
function plot_global_convergence(x, r; save_dir::String="test")
    initialize_style()
    
    # Calculate absolute errors
    epsilon = @. abs(x - r)
    
    # Filter out exact zeros (usually the very last iteration) to avoid log10(0)
    valid_indices = epsilon .> 0
    epsilon = epsilon[valid_indices]
    
    # Calculate a_n = -log10(epsilon) as requested by the prompt
    a_n = @. -log10(epsilon)
    
    # Create iteration array n (x-axis)
    n_iter = collect(1:length(a_n))
    
    # Fit a line to find the slope
    # Since epsilon_n ≈ C^n * epsilon_0, then a_n ≈ -n*log10(C) - log10(epsilon_0)
    # So the slope of a_n vs n is -log10(C)
    M = [ones(length(n_iter)) n_iter]
    intercept, slope = solve_least_squares(M, a_n)
    
    # Calculate C from the slope
    C_global = 10^(-slope)
    
    println("Global Analysis based on n vs a_n:")
    println("Slope ≈ $slope")
    println("Global Asymptotic error constant C ≈ $C_global")
    
    global_plot = scatter(
        n_iter, 
        a_n,
        label=L"Data $(n, a_n)$",
        xlabel=L"Iteration $n$",
        ylabel=L"a_n = -\log_{10}|\epsilon_n|",
        title="Global Convergence Analysis",
        legend=:topleft,
        markershape=:circle,
        markerstrokewidth=0
    )
    
    # Calculate the y-values for the fitted line
    y_fit = @. intercept + slope * n_iter
    
    plot!(
        global_plot,
        n_iter,
        y_fit,
        label=L"Linear Fit (Slope $\approx$ $(round(slope, digits=3)))"
    )
    
    save_graph(global_plot, "global-convergence-an", save_dir)
    display(global_plot)
    readline()
    
    return global_plot
end

end
