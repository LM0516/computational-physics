include("../../modules/linear_systemsV2.jl")
using .LinearSystemsV2
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
    display(X)

    a, α = solve_least_squares(Matrix{Float64}(X), Vector{Float64}(y))
    c = exp(a)

    # --- Plot in log-log scale ---
    p1 = scatter(x, y, label=L"Data (\log)", xlabel=L"\ln R", ylabel=L"\ln \tau", legend=:topleft)
    plot!(x, X * [a, α], label="Linear fit", lw=2)

    # --- Plot in original scale ---
    R_fit = range(minimum(R), maximum(R), length=200)
    τ_fit = @. c * R_fit^α

    p2 = plot(R, τ, seriestype=:scatter, label="Data", xlabel=L"R \, (\mathrm{Mkm})", ylabel=L"\tau \, (\mathrm{days})")
    plot!(p2, R_fit, τ_fit, label=latexstring("Fit \$\\tau = $(round(c, sigdigits=3)) R^{$(round(α, digits=3))}\$"), lw=2)

    plt = plot(
        p1, p2,
        layout=(1, 2),
        size=(1400, 1000),
    )
    display(plt)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
