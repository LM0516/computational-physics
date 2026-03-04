include("../../modules/interpolations.jl")
using .Interpolations
using Plots
using LaTeXStrings

function plot_interpolation(n_values::Array{Int}, a::Int, b::Int, f::Function, foo)
    # Create evaluation points
    x_eval_range = range(a, b, length=200)

    # Plot for f
    p = plot(x_eval_range, f.(x_eval_range),
        label="y = $foo", linewidth=1,
        xlabel="x", ylabel="y", title="Barycentric Lagrange Interpolation",
        legend=:bottomright)

    for n in n_values
        # Create equally spaced nodes
        x_nodes = collect(range(a, b, length=n))

        # Evaluate interpolation at many points
        p_interp = [barycentric_lagrange(x, x_nodes, f, type="lagrange") for x in x_eval_range]

        # Plot interpolation
        plot!(p, x_eval_range, p_interp,
            label="n = $n nodes", linewidth=1, linestyle=:dash)

        # Plot nodes
        scatter!(p, x_nodes, f.(x_nodes),
            label="", markersize=4, color=:black, alpha=0.5)
    end

    display(p)

    readline()
end

function main()
    f_a(x) = log(x)
    f_b(x) = tanh(x)
    f_c(x) = cosh(x)
    f_d(x) = abs(x)

    plot_interpolation([2, 3, 4], 1, 10, f_a, "log(x)")
    plot_interpolation([2, 3, 4], -3, 2, f_b, "tanh(x)")
    plot_interpolation([2, 3, 4], -1, 3, f_c, "cosh(x)")
    plot_interpolation([3, 5, 7], -2, 1, f_d, "|x|")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
