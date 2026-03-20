using ComputationalPhysics
using LinearAlgebra
using Plots
using LaTeXStrings

function main()
    g(t) = exp.(sin.(t .- 1))
    t = [(2π * i) / 60 for i in 1:60]

    b = g(t)

    # === Polynomial Fit ===
    # y(t) = c1 + c2*t + ... + c7*t^6
    P = hcat([t .^ i for i in 0:6]...)  # columns: 1, t, t^2, ..., t^6
    coeffs = solve_least_squares(Matrix{Float64}(P), Vector{Float64}(b))

    println("\nPolynomial fit coefficients: ", coeffs)

    # === FT Fit ===
    # y(t) = d1 + d2 cos(t) + d3 sin(t) + d4 cos(2t) + d5 sin(2t)
    # build F with columns for each basis function
    F = hcat(ones(length(t)), cos.(t), sin.(t), cos.(2t), sin.(2t))
    coeffs_ft = solve_least_squares(Matrix{Float64}(F), Vector{Float64}(b))

    println("\nFourier Transform coefficients: ", coeffs_ft)

    # === Plot ===
    t_plot = range(minimum(t), stop=maximum(t), length=400)
    y_pol = [sum(coeffs[j+1] * x^j for j in 0:6) for x in t_plot]
    F_plot = hcat(ones(length(t_plot)), cos.(t_plot), sin.(t_plot), cos.(2 .* t_plot), sin.(2 .* t_plot))
    y_ft = F_plot * coeffs_ft
    y_true = g(t_plot)

    p1 = plot_generic(t_plot, y_true, label=L"g(t)=e^{\sin(t-1)}", lw=2, xlabel="t", ylabel="y", title="Least squares fits")
    scatter_add!(p1, t, b, label="data", ms=4)
    plot_add!(p1, t_plot, y_pol, label="polynomial fit (deg 6)", lw=2, ls=:dash)
    plot_add!(p1, t_plot, y_ft, label="FT fit (2 harmonics)", lw=2, ls=:dash)

    chi2_poly = chi_square(y_pol, y_true)
    chi2_ft = chi_square(y_ft, y_true)
    dof_poly = length(y_pol) - length(coeffs)
    dof_ft = length(y_ft) - length(coeffs_ft)
    chi2r_poly = chi2_poly / dof_poly
    chi2r_ft = chi2_ft / dof_ft
    p_poly = p_value(chi2_poly, dof_poly)
    p_ft = p_value(chi2_ft, dof_ft)

    println()
    println("=== Polynomial Fit Analysis ===")
    println()
    println("Chi-square: ", round(chi2_poly, digits=4))
    println("Degrees of freedom: ", dof_poly)
    println("Reduced chi-square: ", round(chi2r_poly, digits=4))
    println("P-value: ", round(p_poly, digits=4))

    println("=== Fourier Transform Fit Analysis ===")
    println()
    println("Chi-square: ", round(chi2_ft, digits=4))
    println("Degrees of freedom: ", dof_ft )
    println("Reduced chi-square: ", round(chi2r_ft, digits=4))
    println("P-value: ", round(p_ft, digits=4))

    save_plot(p1, "multi-fit", "2-6")
    display(p1)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
