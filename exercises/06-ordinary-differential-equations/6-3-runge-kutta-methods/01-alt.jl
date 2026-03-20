include("../../modules/ode.jl")
using .ODE
using Plots

function main()
    f = (t, u) -> -2t * u
    g = t -> 2 * exp(-t^2)
    u0 = 2
    a = 0
    b = 2
    errors_rk4 = zeros(10)
    errors_ie2 = zeros(10)

    for (i, n) in enumerate(30:30:300)
        t_rk4, u_rk4 = rk4(f, a, b, n, u0)
        t_ie2, u_ie2 = ie2(f, a, b, n, u0)

        errors_rk4[i] = maximum(abs.(u_rk4 .- g.(t_rk4)))
        errors_ie2[i] = maximum(abs.(u_ie2 .- g.(t_ie2)))
    end

    p = plot(errors_rk4, label="RK4", xscale=:log10, yscale=:log10, marker=:circle)
    plot!(p, errors_ie2, label="IE2", xscale=:log10, yscale=:log10, marker=:square)

    display(p)
    println("Press Enter to exit...")
    readline()
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
