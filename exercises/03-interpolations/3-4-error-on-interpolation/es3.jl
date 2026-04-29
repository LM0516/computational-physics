using ComputationalPhysics
using LinearAlgebra
using Plots
using LaTeXStrings

global cheb = ComputationalPhysics.Chebyshev()
global lag = ComputationalPhysics.Lagrange()

function plot_interpolation(n_values::Int, a::Float64, b::Float64, f::Function, foo; interp_type=cheb)
    # Create evaluation points
    x_eval_range = range(a, b, length=1000)

    # Set plot title based on interpolation type
    title_str = interp_type == cheb ? "Barycentric Chebyshev Interpolation" : "Barycentric Lagrange Interpolation"

    # Plot for f in log-scale
    p1 = plot_generic(x_eval_range, make_log_safe(f.(x_eval_range)),
        label="y = $foo", linewidth=1,
        xlabel="x", ylabel=L"\log(y)", title=title_str,
        yaxis=:log, legend=true)

    # Create nodes based on interpolation type
    if interp_type == cheb
        # Chebyshev nodes in [-1, 1]
        x_nodes = [cos((2 * i - 1) * π / (2 * n_values)) for i in 1:n_values]
        # Transform from [-1,1] to [a,b]
        x_nodes = @. (b - a) / 2 * x_nodes + (a + b) / 2
    else  # lagrange (equally spaced)
        x_nodes = collect(range(a, b, length=n_values))
    end

    # Evaluate interpolation at many points
    p_interp = [barycentric_lagrange(x, x_nodes, f, method=interp_type) for x in x_eval_range]

    # Plot interpolation
    plot_add!(p1, x_eval_range, make_log_safe(p_interp),
        label="n = $n_values nodes", linewidth=1, linestyle=:dash)

    # Plot nodes
    scatter_add!(p1, x_nodes, f.(x_nodes),
        label="", markersize=4, color=:black, alpha=0.5)

    # Maximum error
    err = @. abs(f.(x_eval_range) - p_interp)
    max_error = maximum(err)

    return p1, max_error
end

function main()
    a = 0.0
    b = 2 * π
    n = 40

    f(x) = cosh(sin(x))
    p, max_error = plot_interpolation(n, a, b, f, L"\cosh(\sin(x))", interp_type=cheb)

    println("Maximum error: $max_error")

    save_plot(p, "cosh(sinh(x)) graph", "3-4")
    display(p)
    readline()
end


if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
