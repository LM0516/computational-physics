include("../../modules/interpolations.jl")
using .Interpolations
using Plots
using LaTeXStrings

function analyze_smoothness()
    """
    Part (a): Analyze how many continuous derivatives f_m(x) = |x|^m has
    """
    m_values = [1, 3, 5, 7, 9, 11]
    
    println("\n" * "="^70)
    println("PART (a): Smoothness Analysis of f_m(x) = |x|^m")
    println("="^70)
    
    println("\nFor f_m(x) = |x|^m, the function has ν = floor(m) continuous derivatives.")
    println("The νth derivative is discontinuous at x=0.\n")
    
    for m in m_values
        nu = floor(Int, m)
        println("m = $m:  f_m has $nu continuous derivatives ($(nu)th deriv discontinuous at x=0)")
    end
    println()
end

function estimate_convergence_rate(n_values, errors)
    """
    Estimate the convergence rate by fitting log(error) vs log(n)
    Returns the slope (which represents the power law exponent)
    """
    log_n = log.(n_values)
    log_errors = log.(errors)
    
    # Linear fit: log(error) = a + slope * log(n)
    # slope ≈ -ν for algebraic convergence O(n^{-ν})
    slope = (log_errors[end] - log_errors[1]) / (log_n[end] - log_n[1])
    
    return slope
end

function main()
    # Part (a): Analyze smoothness
    analyze_smoothness()
    
    # Part (b): Compute convergence and plot
    println("="^70)
    println("PART (b): Computing convergence rates...")
    println("="^70)
    
    a = -1
    b = 1
    n_values = collect(10:10:100)
    m_values = [1, 3, 5, 7, 9, 11]
    
    # Create single log-log plot for all m values
    p_errors = plot(
        xlabel="n", 
        ylabel=L"||f_m - p||_{\infty}",
        title="Convergence for |x|^m with Chebyshev Nodes (Algebraic Convergence)",
        legend=:bottomleft,
        size=(1000, 700),
        xaxis=:log,      # Log scale on x-axis
        yaxis=:log       # Log scale on y-axis
    )
    
    println("\nConvergence rates (slopes on log-log plot):\n")
    println("m  | ν (continuous derivs) | Observed slope | Predicted slope")
    println("-"^60)
    
    for m in m_values
        f_m(x) = abs(x)^m
        nu = floor(Int, m)
        
        errors = Float64[]
        
        # Compute error for each n
        for n in n_values
            # Generate Chebyshev nodes of second kind on [-1, 1]
            # Chebyshev 2:
            x_nodes = [- cos((i - 1)/(n - 1)*π) for i in 1:n]
            # Chebyshev 1:
            #=x_nodes = [- cos((2*i - 1) * π / (2*n)) for i in 1:n]=#
            
            # Transform from [-1,1] to [a,b]
            x_nodes = @. (b - a) / 2 * x_nodes + (a + b) / 2
            
            # Evaluate interpolation at many points
            x_eval_range = range(a, b, length=4000)
            p_interp = [barycentric_lagrange(x, x_nodes, f_m, type="chebyshev") 
                       for x in x_eval_range]
            
            # Compute infinity norm
            max_error = maximum(abs.(f_m.(x_eval_range) - p_interp))
            push!(errors, max_error)
        end
        
        # Estimate convergence rate
        observed_slope = estimate_convergence_rate(n_values, errors)
        predicted_slope = -nu  # Theory predicts O(n^{-ν})
        
        println("$m  | $nu                      | $(round(observed_slope, digits=3))           | $predicted_slope")
        
        # Plot this curve
        plot!(p_errors, n_values, errors,
              label="m = $m (predicted slope ≈ $predicted_slope)",
              marker=:circle,
              linewidth=2.5,
              markersize=3)
    end
    
    println()
    
    display(p_errors)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
