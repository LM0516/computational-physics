using ComputationalPhysics
using Plots
using LaTeXStrings

function main()
    f = (t, u) -> -2t * u
    g = t -> 2 * exp(-t^2)
    u0 = 2
    a = 0
    b = 2

    ns = 30:30:300
    errors_rk4 = Vector{Float64}(undef, length(ns))
    errors_ie2 = Vector{Float64}(undef, length(ns))
    evals_rk4 = Vector{Float64}(undef, length(ns))
    evals_ie2 = Vector{Float64}(undef, length(ns))

    for (i, n) in enumerate(ns)
        t_rk4, u_rk4 = rk4(f, a, b, n, u0)
        t_ie2, u_ie2 = ie2(f, a, b, n, u0)

        errors_rk4[i] = maximum(abs.(u_rk4 .- g.(t_rk4)))
        errors_ie2[i] = maximum(abs.(u_ie2 .- g.(t_ie2)))

        # Number of function evaluations
        evals_rk4[i] = 4 * n
        evals_ie2[i] = 2 * n
    end

    p = plot_generic(evals_rk4, errors_rk4, label="RK4", xscale=:log10, yscale=:log10, marker=:circle)
    plot_add!(p, evals_ie2, errors_ie2, label="IE2", xscale=:log10, yscale=:log10, marker=:square)

    # Reference lines
    # For IE2 (2nd order), error ~ C * h^2 ~ C * (1/evals)^2. So log(error) ~ -2 log(evals)
    # For RK4 (4th order), error ~ C * h^4 ~ C * (1/evals)^4. So log(error) ~ -4 log(evals)

    ref_ie2 = errors_ie2[1] * (evals_ie2[1] ./ evals_ie2) .^ 2
    ref_rk4 = errors_rk4[1] * (evals_rk4[1] ./ evals_rk4) .^ 4

    plot_add!(p, evals_ie2, ref_ie2, label=L"O(h^2)", linestyle=:dash, color=:black)
    plot_add!(p, evals_rk4, ref_rk4, label=L"O(h^4)", linestyle=:dot, color=:black)

    xlabel!("Number of function evaluations")
    ylabel!("Max Error")
    title!("Convergence of IE2 and RK4")

    save_plot(p, "IE2-vs-RK4", "6-3")
    display(p)
    readline()
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
