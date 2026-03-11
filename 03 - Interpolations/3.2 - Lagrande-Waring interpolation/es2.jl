include("../../modules/interpolations.jl")
include("../../modules/ConsistentPlots.jl")

using .Interpolations
using .ConsistentPlots
using Plots
using LaTeXStrings

function plot_interpolation(n_values::Array{Int}, a::Int, b::Int, f::Function, function_name)
    println("=== Evaluating $function_name ===")
    # Create evaluation points
    x_eval_range = range(a, b, length=200)
    chi2 = Vector{Float64}(undef, length(n_values))


    # Plot for f
    p1 = plot(x_eval_range, f.(x_eval_range), label=L"y = %$function_name", xlabel=L"x", ylabel=L"y")
    p2 = plot(linewidth=1, yscale=:log10, xlabel=L"x", ylabel=L"Absolute Error $|P(x) - f(x)|$")

    for (i, n) in enumerate(n_values)
        # Create equally spaced nodes
        x_nodes = collect(range(a, b, length=n))

        # Evaluate interpolation at many points
        p_interp = [barycentric_lagrange(x, x_nodes, f, type="lagrange") for x in x_eval_range]
        p_real = f.(x_eval_range)
        chi2[i] = chi_square(p_interp, p_real)

        # Plot interpolation
        plot_add!(p1, x_eval_range, p_interp, label=L"$n = %$n$ nodes", linewidth=1, linestyle=:dash)

        # Plot nodes
        scatter!(p1, x_nodes, f.(x_nodes),
            label="", markersize=4, color=:black, alpha=0.5)

        # Calculate absolute residuals
        # We add eps(Float64) (machine epsilon) to avoid taking log10(0) at the nodes
        residuals = abs.(p_interp .- p_real) .+ eps(Float64)
        
        plot_add!(p2, x_eval_range, residuals, label=L"$n = %$n$ error")
    end

    p = plot(p1, p2, layout = (1, 2), size = (900, 400))
    println("Chi squares: ", chi2)
    function_name == "|x|" ? save_graph(p, "abs_x-plot", "3-2") : save_graph(p, "$function_name-plot", "3-2")
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
