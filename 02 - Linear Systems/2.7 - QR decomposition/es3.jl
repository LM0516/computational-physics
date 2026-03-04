include("../../modules/linear_systemsV2.jl")
using .LinearSystemsV2
using LinearAlgebra
using Plots
using LaTeXStrings
using Printf

function main()
    println("="^80)
    println("EXERCISE 3: Least Squares Polynomial Fitting")
    println("="^80)
    println()
    
    # Define function and generate data
    h(t) = exp(sin(4*t)) / 2006.787453080206 
    t = [i/99 for i in 0:99]  # 100 data points from 0 to 1
    
    b = h.(t)  # Apply function element-wise
    
    # === Construct Vandermonde matrix for polynomial fit ===
    # y(t) = c1 + c2*t + ... + c15*t^14 (15 coefficients, degree 14)
    P = hcat([t .^ i for i in 0:14]...)  # columns: 1, t, t^2, ..., t^14
    
    println("Number of data points: $(length(t))")
    println("Polynomial degree: 14 (15 coefficients: c₁, c₂, ..., c₁₅)")
    println()
    
    # === Method (a): LU Factorization via Normal Equations ===
    println("METHOD (a): Normal Equations with LU/Cholesky")
    println("-"^80)
    coeffs_lu = solve_least_squares(Matrix{Float64}(P), Vector{Float64}(b))
    @printf("Coefficient c₁₅ = %.15f\n", coeffs_lu[15])
    @printf("Error: %.2e\n", abs(coeffs_lu[15] - 1.0))
    println()
    
    # === Method (b): QR Decomposition ===
    println("METHOD (b): QR Decomposition")
    println("-"^80)
    Q, R = qr_mgs(Matrix{Float64}(P))
    z = Q' * b
    coeffs_qr = R \ z
    @printf("Coefficient c₁₅ = %.15f\n", coeffs_qr[15])
    @printf("Error: %.2e\n", abs(coeffs_qr[15] - 1.0))
    println()
    
    # === Method (c): Q-less QR Decomposition ===
    println("METHOD (c): Q-less QR Decomposition (Most Stable)")
    println("-"^80)
    R_qless, z_qless = qr_mgs_augmented(Matrix{Float64}(P), Vector{Float64}(b))
    coeffs_qless = R_qless \ z_qless
    @printf("Coefficient c₁₅ = %.15f\n", coeffs_qless[15])
    @printf("Error: %.2e\n", abs(coeffs_qless[15] - 1.0))
    println()
    
    # === Summary ===
    println("="^80)
    println("COMPARISON SUMMARY")
    println("="^80)
    println("Method                          c₁₅ value           Absolute Error")
    println("-"^80)
    @printf("(a) Normal Equations            %.15f    %.2e\n", coeffs_lu[15], abs(coeffs_lu[15] - 1.0))
    @printf("(b) Thin QR                     %.15f    %.2e\n", coeffs_qr[15], abs(coeffs_qr[15] - 1.0))
    @printf("(c) Q-less QR (Recommended)     %.15f    %.2e\n", coeffs_qless[15], abs(coeffs_qless[15] - 1.0))
    println("="^80)
    println()
    
    # === Plot ===
    t_plot = range(0, stop=1, length=400)
    y_true = h.(t_plot)
    
    # Evaluate polynomial fits
    y_lu = [sum(coeffs_lu[j+1] * x^j for j in 0:14) for x in t_plot]
    y_qr = [sum(coeffs_qr[j+1] * x^j for j in 0:14) for x in t_plot]
    y_qless = [sum(coeffs_qless[j+1] * x^j for j in 0:14) for x in t_plot]
    
    p1 = plot(t_plot, y_true, label=L"h(t)=e^{\sin(4t)}/2006.78", lw=2, 
              xlabel="t", ylabel="y", title="Least Squares Polynomial Fits (degree 14)",
              legend=:topleft)
    scatter!(p1, t, b, label="data points", ms=3, alpha=0.6)
    
    plot!(p1, t_plot, y_lu, label="(a) Normal Equations", lw=1.5, ls=:dash)
    plot!(p1, t_plot, y_qr, label="(b) QR Decomposition", lw=1.5, ls=:dot)
    plot!(p1, t_plot, y_qless, label="(c) Q-less QR", lw=1.5, ls=:dashdot)
    
    display(p1)

    readline()
 
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
