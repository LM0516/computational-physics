include("../../modules/linear_systemsV2.jl")
include("../../modules/interpolations.jl")
include("../../modules/ConsistentPlots.jl")

using .LinearSystemsV2
using .Interpolations
using .ConsistentPlots
using LinearAlgebra
using Plots
using LaTeXStrings

function LU_decomposition(P::Matrix{Float64}, b::Vector{Float64})
    # === LU Factorization via Normal Equations ===
    println("Normal Equations with LU")
    println("-"^80)
    coeffs_lu = solve_least_squares(Matrix{Float64}(P), Vector{Float64}(b))
    println("Coefficient c₁₅ = : ", coeffs_lu[15])
    println("Error: ", abs(coeffs_lu[15] - 1.0))
    println()
    return coeffs_lu
end

function QR_decomposition(P::Matrix{Float64}, b::Vector{Float64})
    # === QR Decomposition ===
    println("QR Decomposition")
    println("-"^80)
    Q, R = qr_mgs(Matrix{Float64}(P))
    z = Q' * b
    coeffs_qr = R \ z
    println("Coefficient c₁₅ = ", coeffs_qr[15])
    println("Error: ", abs(coeffs_qr[15] - 1.0))
    println()
    return coeffs_qr
end

function Q_less(P::Matrix{Float64}, b::Vector{Float64})
    # === Q-less QR Decomposition ===
    println("Q-less QR Decomposition")
    println("-"^80)
    R_qless, z_qless = qr_mgs_augmented(Matrix{Float64}(P), Vector{Float64}(b))
    coeffs_qless = R_qless \ z_qless
    println("Coefficient c₁₅ = ", coeffs_qless[15])
    println("Error: ", abs(coeffs_qless[15] - 1.0))
    println()
    return coeffs_qless
end

function main()
    println("="^80)
    println("Least Squares Polynomial Fitting")
    println("="^80)
    println()

    h(t) = exp(sin(4*t)) / 2006.787453080206 
    t = [i/99 for i in 0:99]

    b = h.(t)

    # === Construct Vandermonde matrix for polynomial fit ===
    # y(t) = c1 + c2*t + ... + c15*t^14 (15 coefficients, degree 14)
    P = hcat([t .^ i for i in 0:14]...)

    # Fits
    coeffs_lu = LU_decomposition(P, b)
    coeffs_qr = QR_decomposition(P, b)
    coeffs_qless = Q_less(P, b)

    # === Plot ===
    t_plot = range(0, stop=1, length=400)
    y_true = h.(t_plot)
    
    # Evaluate polynomial fits
    y_lu = [sum(coeffs_lu[j+1] * x^j for j in 0:14) for x in t_plot]
    y_qr = [sum(coeffs_qr[j+1] * x^j for j in 0:14) for x in t_plot]
    y_qless = [sum(coeffs_qless[j+1] * x^j for j in 0:14) for x in t_plot]

    # Fit goodness
    fit_goodness(y_lu, y_true, coeffs_lu)
    fit_goodness(y_qr, y_true, coeffs_qr)
    fit_goodness(y_qless, y_true, coeffs_qless)
    
    p1 = plot_generic(t_plot, y_true, label=L"h(t)=e^{\sin(4t)}/2006.78", lw=2, 
              xlabel="t", ylabel="y", title="Least Squares Polynomial Fits (degree 14)")
    scatter!(p1, t, b, label="data points", ms=3, alpha=0.6)
    
    plot_add!(p1, t_plot, y_lu, label="(a) Normal Equations", lw=1.5, ls=:dash)
    plot_add!(p1, t_plot, y_qr, label="(b) QR Decomposition", lw=1.5, ls=:dot)
    plot_add!(p1, t_plot, y_qless, label="(c) Q-less QR", lw=1.5, ls=:dashdot)
    
    save_graph(p1, "multiple-least-square", "2-7")
    display(p1)

    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    initialize_style()
    main()
end
