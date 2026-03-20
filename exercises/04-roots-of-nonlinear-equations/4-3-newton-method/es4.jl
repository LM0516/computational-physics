using ComputationalPhysics
using Plots
using LaTeXStrings

# TODO: Check if this exercise is correct, the plot is strange
function main()
    τ = 1.7610
    ϵ = 0.2230
    a = 0
    b = 2π
    θ = zeros(100)

    for (i, t) in enumerate(range(0, τ, length=100))
        f = ψ -> ψ - ϵ * sin(ψ) - 2π * t / τ
        df = ψ -> 1 - ϵ * cos(ψ)
        sol = newton_method(f, df, a, b)[1]
        θ[i] = mod(2 * atan(sqrt((1 + ϵ) / (1 - ϵ)) * tan(sol / 2)), 2π)
        println("t = $t, sol = $sol")
    end

    p = plot_generic(range(0, τ, length=100), θ, label=L"\theta(t)")
    xlabel!(p, L"t")
    ylabel!(p, L"\theta(t)")
    save_plot(p, "full-orbit", "4-3")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
