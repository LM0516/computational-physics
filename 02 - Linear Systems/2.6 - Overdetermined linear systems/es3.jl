include("../../modules/linear_systemsV2.jl")
using .LinearSystemsV2
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

    # === FFT Fit ===
    # y(t) = d1 + d2 cos(t) + d3 sin(t) + d4 cos(2t) + d5 sin(2t)
    # build F with columns for each basis function
    F = hcat(ones(length(t)), cos.(t), sin.(t), cos.(2t), sin.(2t))
    coeffs_fft = solve_least_squares(Matrix{Float64}(F), Vector{Float64}(b))

    # === Plot ===
    t_plot = range(minimum(t), stop=maximum(t), length=400)
    y_pol = [sum(coeffs[j+1] * x^j for j in 0:6) for x in t_plot]
    # evaluate FFT fit on the dense plotting grid
    F_plot = hcat(ones(length(t_plot)), cos.(t_plot), sin.(t_plot), cos.(2 .* t_plot), sin.(2 .* t_plot))
    y_fft = F_plot * coeffs_fft
    y_true = g(t_plot)

    p1 = plot(t_plot, y_true, label=L"g(t)=e^{\sin(t-1)}", lw=2, xlabel="t", ylabel="y", title="Least squares fits")
    scatter!(p1, t, b, label="data", ms=4)
    plot!(p1, t_plot, y_pol, label="polynomial fit (deg 6)", lw=2, ls=:dash)
    plot!(p1, t_plot, y_fft, label="FFT fit (2 harmonics)", lw=2, ls=:dash)

    # Chi-Square
    χ2_poly = sum(((y_pol .- y_true) .^ 2) ./ y_true)
    χ2_fft = sum(((y_fft .- y_true) .^ 2) ./ y_true)
    @show χ2_poly χ2_fft

    display(p1)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
