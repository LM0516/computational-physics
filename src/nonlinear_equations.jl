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
    plot_init()
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

    save_plot(p, label, save_dir)
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
The convergence rate is a measure of how much an error shrinks with each step.

Convergence rate formula:
q = ln(d_{k+1} / d_k) / (ln(d_k / d_{k-1}))
Asymptotic error constant formula:
C = (|ϵ_{k+1}|) / (|ϵ_k|^q = (|x_{k+1} - r|) / (|x_k - r|^q)
"""
function conv_rate_and_asymp_const(x_history::Array{Float64}, root::Float64)
    # === Convergence rates (Algorithm 2: using successive differences) ===
    d = abs.(diff(x_history))

    # Truncate array instead of filtering to maintain sequence integrity.
    # This can avoid cancellation errors when calculating `q_values`
    valid_idx = length(d)
    for i in 1:(length(d)-1)
        if d[i] < 1e-7 || d[i+1] >= d[i]
            valid_idx = i
            break
        end
    end

    # Keep only the valid, monotonically decreasing part of the sequence
    d_valid = d[1:valid_idx]
    # @show d_valid
    n = length(d_valid)

    if n < 3
        @warn "Not enough iterations to estimate convergence rate."
        return [NaN], [NaN]
    end

    # Pre-allocate memory for q_values
    num_q = n - 2
    q_values = Vector{Float64}(undef, num_q)

    for k in 1:num_q
        idx = k + 1
        q_values[k] = log(d_valid[idx+1] / d_valid[idx]) / log(d_valid[idx] / d_valid[idx-1])
    end

    # === Asymptotic constant ===
    errors = abs.(x_history .- root)

    # Keep only errors above the numerical noise threshold
    valid_errors = filter(e -> e > 1e-14, errors)
    m = length(valid_errors)

    if m < 2
        return q_values, [NaN]
    end

    # Pre-allocate memory for C_values
    num_c = m - 1
    C_values = Vector{Float64}(undef, num_c)

    last_q = q_values[end]

    for k in 1:num_c
        C_values[k] = valid_errors[k+1] / (valid_errors[k]^last_q)
    end

    return q_values, C_values
end

"""
Generates a plot to visually analyze the convergence of an iterative method.
It plots `log10(|ε_{n+1}|)` against `log10(|ε_n|)` and fits a line to determine
the order of convergence and the asymptotic error constant.
"""
function plot_convergence_analysis(x, r; skip_initial::Int=0, save_dir::String="output/figures")
    plot_init()
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

    conv_plot = scatter_generic(
        log_eps_n,
        log_eps_n1,
        label=L"Data ($\log_{10}|\epsilon_{n+1}|) vs (\log_{10}|\epsilon_n|))$",
        xlabel=L"\log_{10}|\epsilon_n|",
        ylabel=L"\log_{10}|\epsilon_{n+1}|",
        # title="Convergence Analysis (Skipping first $skip_initial iterations)",
        legend=:topleft,
        markershape=:circle,
        markerstrokewidth=0
    )

    # Calculate the y-values for the fitted line
    y_fit = @. log10_C + q * log_eps_n

    plot_add!(
        conv_plot,
        log_eps_n,
        y_fit,
        label="Least-Squares Fit (slope q ≈ $(round(q, digits=3)))",
    )

    save_plot(conv_plot, "conv-analysis-with-$skip_initial-skipped", save_dir)
    display(conv_plot)
    readline()
    return conv_plot
end

"""
    newton_method(f, df, a, b; kwargs...) -> (root, trajectory)

Finds a root of the function `f` within the interval `[a, b]` using the Newton-Raphson 
algorithm, with an optional fallback to the bisection method.

The Newton iteration follows the formula:
\$\$x_{n+1} = x_n - \\frac{f(x_n)}{f'(x_n)}\$\$

# Arguments
- `f::Function`: The objective function to find the root for.
- `df::Function`: The first derivative of the function `f`.
- `a::Real`, `b::Real`: The lower and upper bounds of the search interval.

# Keywords
- `x_init::Real=b`: The starting point for the Newton iteration.
- `xtol::Float64=1e-13`: The tolerance for the change in x (|x_curr - x_prev|).
- `ftol::Float64=1e-13`: The tolerance for the function value (|f(x)|).
- `max_iter::Int=100`: The maximum number of Newton iterations allowed.
- `bracketing::Bool=true`: If true, triggers a `bisection` fallback if Newton's 
  method iterates outside the bounds `[a, b]`.

# Returns
- `root::Float64`: The estimated root of the function.
- `trajectory::Vector{Float64}`: The sequence of x-values generated during the iteration.

# Note
If Newton's method fails to converge within `max_iter`, a warning is issued. If 
`bracketing` is enabled and the method diverges outside the interval `[a, b]`, 
the function automatically calls `bisection(f, a, b)` to guarantee a solution.
"""
function newton_method(f::Function, df::Function, a::Real, b::Real; x_init::Real=b, xtol::Float64=1e-13, ftol::Float64=1e-13, max_iter::Int=100, bracketing::Bool=true)

    x_prev = float(x_init)
    x_curr = x_prev - f(x_prev) / df(x_prev)

    trajectory = Vector{Float64}(undef, max_iter + 1)
    trajectory[1] = x_prev
    trajectory[2] = x_curr
    iter = 1

    while iter < max_iter && (abs(x_curr - x_prev) >= xtol || abs(f(x_curr)) >= ftol)
        x_prev = x_curr
        x_curr = x_prev - f(x_prev) / df(x_prev)
        iter += 1
        trajectory[iter+1] = x_curr
    end

    if iter == max_iter
        @warn "Maximum number of iterations reached at x=$x_curr"
        @warn "Residual f(x) = $(f(x_curr))"
    end

    if bracketing && (x_curr < a || x_curr > b)
        println("Newton diverged outside [$a, $b], using bisection fallback...")
        x_curr = bisection(f, a, b)
    end

    return x_curr, trajectory[1:iter+1]
end


"""
Finds a root of the function `f` using the secant method.
"""
function secant_method(f::Function, x_0::Real, x_1::Real; xtol::Float64=1e-14, ftol::Float64=1e-14, max_iter::Int=100, bracketing::Bool=true)
    x_hist = Vector{Float64}(undef, max_iter + 2)
    x_hist[1] = float(x_0)
    x_hist[2] = float(x_1)

    f_prev = f(x_hist[1])
    f_curr = f(x_hist[2])

    # Check if we are already at a root before looping
    if abs(f_curr) < ftol
        return x_hist[2], x_hist[1:2]
    end

    curr_idx = 2
    for iter in 1:max_iter
        diff_f = f_curr - f_prev
        # Avoid division by zero
        if abs(diff_f) < 1e-18
            break
        end

        x_next = x_hist[curr_idx] - f_curr * (x_hist[curr_idx] - x_hist[curr_idx-1]) / diff_f
        curr_idx += 1
        x_hist[curr_idx] = x_next

        # Check convergence criteria
        if abs(x_next - x_hist[curr_idx-1]) < xtol || abs(f(x_next)) < ftol
            break
        end

        # Update values for next iteration
        f_prev = f_curr
        f_curr = f(x_next)
    end

    # Trim the history to the actual number of iterations performed
    final_history = x_hist[1:curr_idx]
    x_final = final_history[end]

    # Bracketing check (Bisection fallback)
    if bracketing && (f(x_0) * f(x_1) < 0)
        low, high = min(x_0, x_1), max(x_0, x_1)
        if x_final < low || x_final > high
            # bisection(f, a, b) must be defined elsewhere
            x_final = bisection(f, x_0, x_1)
        end
    end

    return x_final, final_history
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
function plot_global_convergence(x::Vector{Float64}, r::Float64; save_dir::String="test", plot_name::String="global-convergence-an")
    plot_init()

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

    global_plot = scatter_generic(
        n_iter,
        a_n,
        label=L"Data $(n, a_n)$",
        xlabel=L"Iteration $n$",
        ylabel=L"a_n = -\log_{10}|\epsilon_n|",
        # title="Global Convergence Analysis",
        legend=:topleft,
        markershape=:circle,
        markerstrokewidth=0
    )

    # Calculate the y-values for the fitted line
    y_fit = @. intercept + slope * n_iter

    plot_add!(
        global_plot,
        n_iter,
        y_fit,
        label=LaTeXString("Linear Fit (Slope \$\\approx $(round(slope, digits=3))\$)")
    )

    save_plot(global_plot, plot_name, save_dir)
    display(global_plot)
    readline()

    return global_plot
end
