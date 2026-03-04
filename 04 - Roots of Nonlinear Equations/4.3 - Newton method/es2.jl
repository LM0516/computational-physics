include("../../modules/nonlinear_equatioins.jl")
using .NonlinearEquations
using Plots

function solutions(f::Function, df::Function, a, b, f_eqation)
    println("Findig solutions...")
    solution, x = newton_method(f, df, a, b)
    q, C = convergence(x, solution, skip=1)
    @show solution

    println("Solution found! Plotting the data...")

    x_vals = range(a, b, length=200)
    p = plot(x_vals, f.(x_vals), label=f_eqation)
    scatter!([solution], [0.0], label=round(solution; digits=4))
    hline!([0.0], label="y=0", linestyle=:dash)

    return p
end

function main()
    # === x^2 = e^{-x} ===
    f1 = x -> (x^2 - exp(-x))
    df1 = x -> (2x + exp(-x))
    a1 = -2
    b1 = 2

    # === 2x = tan(x) ===
    f2 = x -> (2x - tan(x))
    df2 = x -> (2 - sec(x)^2)
    a2 = 0.2
    b2 = 1.4

    # === e^{x + 1} = 2 + x ===
    f3 = x -> (exp(x + 1) - 2 - x)
    df3 = x -> (exp(x + 1) - 1)
    a3 = -2
    b3 = 2

    # Better Plots
    p1 = solutions(f1, df1, a1, b1, L"y = x^2 - e^{-x}")
    p2 = solutions(f2, df2, a2, b2, L"2x - \tan(x)")
    p3 = solutions(f3, df3, a3, b3, L"e^{x + 1} - 2 - x")

    p = plot(p1, p2, p3, layout=(1, 3), size=(1200, 600))
    display(p)
    readline()

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
