include("../../modules/linear_systemsV2.jl")
include("../../modules/interpolations.jl")

using .LinearSystemsV2
using .Interpolations
using LinearAlgebra
using Plots
using LaTeXStrings

function main()
    a = 1
    chi2 = zeros(10, 10)
    p_value = zeros(10, 10)
    for k in 10:10:100
        b = 1
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

            dof = length(y_fit) - length(coeffs)
            chi2[a][b] = chi_square(y_true, y_fit)
            p_value[a][b] = p_value(chi2[a][b], dof)
            b += 1
        end
        gif(anim, "gif/test_$k.gif", fps=2)
        a += 1
    end
    println("Chi-square: ")
    display(chi2)
    println("p-value: ")
    display(p_value)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
