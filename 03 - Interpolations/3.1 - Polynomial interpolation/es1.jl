include("../../modules/linear_systemsV2.jl")
include("../../modules/interpolations.jl")
include("../../modules/ConsistentPlots.jl")

using .LinearSystemsV2
using .Interpolations
using .ConsistentPlots
using LinearAlgebra
using Plots
using LaTeXStrings

function animations()
    for k in 10:10:100
        anim = @animate for n in 4:4:40
            h = 2 / n
            i = collect(1:n+1)
            x = @. -1 + (i - 1) * h

            f = @. cos(k * x)

            # Vandermonde matrix 
            V = [x[i]^j for i = 1:n+1, j = 0:n]

            coeffs = solve_least_squares(Matrix{Float64}(V), Array{Float64}(f))

            x_fit = range(-1, 1, length=200)
            y_fit = [sum(coeffs[j+1] * xi^j for j = 0:n) for xi in x_fit]
            y_true = @. cos(k * x_fit)

            # Plot
            scatter(x, f, label="Data points", markersize=4)
            plot!(x_fit, y_fit, label="Polynomial fit (degree $n)",
                xlabel="x", ylabel="f(x)", legend=:topright)
        end
        gif(anim, "gif/test_$k.gif", fps=2)
    end
end

function plot_condition_number(ns, condition_numbers)
    p = plot_generic(
        ns, condition_numbers,
        marker = :circle,
        linewidth = 2,
        markersize = 6,
        yscale = :log10,
        xlabel = "n (polynomial degree)",
        ylabel = "Condition Number (log scale)",
        legend = false,
    )
    save_graph(p, "condition-numbers", "3-1")
    display(p)
    readline()
end

function plot_errors(ns, ks, errors, x0_coefs)
    p1 = plot(
        xlabel = "n (polynomial degree)",
        ylabel = "Errors",
    )

    heatmap!(p1, ns, ks, log10.(errors),
        xlabel="n",
        ylabel="k",
        colorbar_title="log10 error"
    )

    p2 = plot(
        xlabel = "n (polynomial degree)",
        ylabel = "x = 0 coefficients",
    )

    for (ki, k) in enumerate(ks)
        plot!(p2, ns, x0_coefs[ki, :],
            marker = :circle,
            linewidth = 1.5,
            markersize = 4,
            label = "k = $k"
        )
    end

    p = plot(p1, p2, layout = (1, 2), size = (900, 400))
    save_graph(p, "errors-plots", "3-1")
    display(p)
    readline()
end

function plot_chi2_heatmap(ns, ks, x0_coefs)
    p = plot(
        xlabel = "n (polynomial degree)",
        ylabel = L"\chi^2",
    )

    heatmap!(p, ns, ks, x0_coefs,
        xlabel="n",
        ylabel="k",
        title="Value of interpolating polynomial at x = 0",
        colorbar_title=L"p_n(0)"
    )

    save_graph(p, "chi-square", "3-1")
    display(p)
    readline()
end

"""
Compute the Vandermonde matrix for a given n
"""
function vandermonde(n)
    # Parameters
    h = 2 / n
    i = collect(1:n+1)
    x = @. -1 + (i - 1) * h

    # Vandermonde Matrix
    V = [x[i]^j for i = 1:n+1, j = 0:n]
    c_number = cond(V)
    return V, c_number
end

function main()
    ns = collect(4:4:40)
    ks = collect(10:10:100)
    # Create an empty Vector of matrices for the Vandermonde matrices
    vandermonde_matrices = Vector{Matrix{Float64}}(undef, length(ns))
    # Create an empty Vector for the condition numbers
    condition_numbers = zeros(length(ns))
    # Creata an empty Matrix for the x = 0 coeffs
    x0_coefs = zeros(length(ks), length(ns))
    # Create an empty Matrix for the chi2
    chi2_values = zeros(length(ks), length(ns))

    # Compose a Vecror of Matrices for every 'n' Vandermonde matrx
    for (j, n) in enumerate(ns)
        vandermonde_matrices[j], condition_numbers[j] = vandermonde(n)
    end

    for (ki, k) in enumerate(ks)
        for (j, n) in enumerate(ns)
            i = collect(1:n+1)
            h = 2 / n
            x = @. -1 + (i - 1) * h
            f = @. cos(k * x)
            coeffs = solve_least_squares(Matrix{Float64}(vandermonde_matrices[j]), Array{Float64}(f))
            x0_coefs[ki, j] = coeffs[1]

            x_fit = range(-1, 1, length=200)
            y_fit = [sum(coeffs[j+1] * xi^j for j = 0:n) for xi in x_fit]
            y_true = @. cos(k * x_fit)

            chi2_values[ki, j] = chi_square(y_fit, y_true)
        end
    end

    # Calculate the errors
    errors = abs.(x0_coefs - ones(length(ks), length(ns)))

    # Plotting the results
    plot_condition_number(ns, condition_numbers)
    plot_errors(ns, ks, errors, x0_coefs)
    plot_chi2_heatmap(ns, ks, x0_coefs)
end

if abspath(PROGRAM_FILE) == @__FILE__
    initialize_style()
    main()
end
