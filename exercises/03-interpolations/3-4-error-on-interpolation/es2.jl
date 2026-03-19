using ComputationalPhysics
using Plots
using LaTeXStrings

function estimate_K(n_values, errors)
    # Take natural log of errors
    log_errors = log.(errors)
    
    # Linear fit: use first and last points for slope
    n1, n_last = n_values[1], n_values[end]
    err1, err_last = log_errors[1], log_errors[end]
    
    slope = (err_last - err1) / (n_last - n1)
    
    # slope = -ln(K), so K = exp(-slope)
    K = exp(-slope)
    
    return K, slope
end

function plot_interpolation(n_values::Array{Int}, a::Int, b::Int, f::Function, function_name; interp_type="chebyshev", fig_name="plot")
    # Create evaluation points
    x_eval_range = range(a, b, length=4000)
    inf_norm = Vector{Float64}(undef, length(n_values))

    # Set plot title based on interpolation type
    title_str = interp_type == "chebyshev" ? "Barycentric Chebyshev Interpolation" : "Barycentric Lagrange Interpolation"
    
    # Plot for f
    # BUG: There is a problem while plotting the log fuctins, if there is a non-positive value 
    # the plot is not dispalyed correctly. Try to came up with a different approach in plotting
    # the log-scale graphs. This is a recurring problem.
    p1 = plot_generic(x_eval_range, make_log_safe(f.(x_eval_range)),
        label="y = $function_name",
        xlabel="x", ylabel=L"\log(y)", title=title_str,
        yaxis=:log)

    for (idx, n) in enumerate(n_values)
        # Create nodes based on interpolation type
        if interp_type == "chebyshev"
            # Chebyshev nodes in [-1, 1]
            x_nodes = [cos((2*i - 1) * π / (2*n)) for i in 1:n]
            # Transform from [-1,1] to [a,b]
            x_nodes = @. (b - a) / 2 * x_nodes + (a + b) / 2
        else  # lagrange (equally spaced)
            x_nodes = collect(range(a, b, length=n))
        end

        # Evaluate interpolation at many points
        p_interp = [barycentric_lagrange(x, x_nodes, f, type=interp_type) for x in x_eval_range]

        # Plot interpolation
        plot_add!(p1, x_eval_range, make_log_safe(p_interp),
            label="n = $n nodes", linestyle=:dot, legend=false,
            seriestype=:steppre, linealpha=:0.5)

        # Compute infinity norm
        norm = @. abs(f(x_eval_range) - p_interp)
        inf_norm[idx] = maximum(norm)
    end

    K, slope = estimate_K(n_values, inf_norm)
    println("Function: $function_name | Type: $interp_type | K ≈ $K (slope ≈ $slope)")

    # Plot infinity norm with log scale
    p2 = plot_generic(n_values, inf_norm,
              linewidth=1, xlabel="n", ylabel=L"||f - p||_{\infty}",
        yaxis=:log, marker=:circle,
        title="Interpolation Error", legend=false)

    # Plot both plots
    p = multi_plot(p1, p2, layout=(1, 2), size=(1200, 600))
    save_plot(p, "$fig_name-$interp_type", "3-4")
    display(p)
    readline()
end

function main()
    a = -1
    b = 1
    n = collect(4:4:60)

    f_a(x) = 1 / (25 * x^2 + 1.0)
    f_b(x) = tanh(5 * x + 2)
    f_c(x) = cosh(sin(x))
    f_d(x) = sin(cosh(x))

    println("=== CHEBYSHEV NODES ===")
    plot_interpolation(n, a, b, f_a, L"\frac{1}{(25x^2 + 1)}", interp_type="chebyshev", fig_name="fa")
    plot_interpolation(n, a, b, f_b, L"\tanh(5x + 2)", interp_type="chebyshev", fig_name="fb")
    plot_interpolation(n, a, b, f_c, L"\cosh(\sin(x))", interp_type="chebyshev", fig_name="fc")
    plot_interpolation(n, a, b, f_d, L"\sin(\cosh(x))", interp_type="chebyshev", fig_name="fd")
    
    println("\n=== EQUIDISTANT NODES ===")
    plot_interpolation(n, a, b, f_a, L"\frac{1}{(25x^2 + 1)}", interp_type="lagrange", fig_name="fa")
    plot_interpolation(n, a, b, f_b, L"\tanh(5x + 2)", interp_type="lagrange", fig_name="fb")
    plot_interpolation(n, a, b, f_c, L"\cosh(\sin(x))", interp_type="lagrange", fig_name="fc")
    plot_interpolation(n, a, b, f_d, L"\sin(\cosh(x))", interp_type="lagrange", fig_name="fd")
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
