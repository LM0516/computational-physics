using ComputationalPhysics
using LinearAlgebra
using Plots
using LaTeXStrings

function main()
    R = [57.59, 108.11, 149.57, 227.84, 778.14, 1427.0, 2870.3, 4499.9] # Distance from Sun
    τ = [87.99, 224.7, 365.26, 686.98, 4332.24, 10759.0, 30684.0, 60188.0] # Orbital period in days
    # Linear system: τ = c * R^α
    x = log.(R)
    y = log.(τ)
    X = [ones(length(R)) log.(R)]
    # display(X)

    a, α = solve_least_squares(Matrix{Float64}(X), Vector{Float64}(y))
    c = exp(a)

    println("a: ", a)
    println("α: ", α)
    println("c: ", c)
    
    # === Chi-square test and p-value for fit ===
    τ_predicted = @. c * R^α
    chi2 = chi_square(τ, τ_predicted)

    dof = length(τ) - 2
    chi2_reduced = chi2 / dof
    p = p_value(chi2, dof)

    println("Chi-square: ", round(chi2, digits=4))
    println("Degrees of freedom: ", dof)
    println("Reduced chi-square: ", round(chi2_reduced, digits=4))
    println("P-value: ", round(p, digits=4))

    # === Plot in log-log scale ===
    p1 = scatter_generic(x, y, label=L"Data (\log)", xlabel=L"\ln R", ylabel=L"\ln \tau", legend=:topleft)
    plot_add!(p1, x, X * [a, α], label="Linear fit", lw=2)

    # === Plot in original scale ===
    R_fit = range(minimum(R), maximum(R), length=200)
    τ_fit = @. c * R_fit^α

    p2 = plot_generic(R, τ, label="Data", seriestype=:scatter, xlabel=L"R \, (\mathrm{Mkm})", ylabel=L"\tau \, (\mathrm{days})")
    plot_add!(p2, R_fit, τ_fit, label=latexstring("Fit \$\\tau = $(round(c, sigdigits=3)) R^{$(round(α, digits=3))}\$"), lw=2)

    plt = multi_plot(
        p1, p2,
        layout=(1, 2),
        size=(1400, 1000),
    )

    save_plot(plt, "kepler-fit", "2-6")
    display(plt)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
