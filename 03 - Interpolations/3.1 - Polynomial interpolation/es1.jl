include("../../modules/linear_systemsV2.jl")
include("../../modules/interpolations.jl")

using .LinearSystemsV2
using .Interpolations
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

    # Debug
    #=display(chi2_values)=#
    #=readline()=#
    #=display(x0_coefs)=#
    #=readline()=#

    p1 = plot(
        ns, condition_numbers,
        marker = :circle,
        linewidth = 2,
        markersize = 6,
        color = :royalblue,
        yscale = :log10,
        xlabel = "n (polynomial degree)",
        ylabel = "Condition Number (log scale)",
        title = "Vandermonde Matrix Condition Number vs n",
        legend = false,
        grid = true,
        framestyle = :box
    )

    p2 = plot(
        xlabel = "n (polynomial degree)",
        ylabel = "χ² (log scale)",
        title = "Least-Squares χ² vs n for cos(kx)",
        #=yscale = :log10,=#
        legend = :topright,
        grid = true,
        framestyle = :box,
        color_palette = :tab10
    )

    for (ki, k) in enumerate(ks)
        plot!(p2, ns, chi2_values[ki, :],
            marker = :circle,
            linewidth = 1.5,
            markersize = 4,
            label = "k = $k"
        )
    end

    p3 = plot(
        xlabel = "n (polynomial degree)",
        ylabel = "x = 0 coefficients",
        title = "Polynomial coefficients for x = 0",
        legend = :topright,
        grid = true,
        framestyle = :box,
        color_palette = :tab10
    )

    for (ki, k) in enumerate(ks)
        plot!(p3, ns, x0_coefs[ki, :],
            marker = :circle,
            linewidth = 1.5,
            markersize = 4,
            label = "k = $k"
        )
    end

    combined = plot(p1, p3, layout = (1, 2), size = (1000, 700), dpi = 150)
    display(combined)
    readline()
    display(p2)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
