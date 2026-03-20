using ComputationalPhysics
using Plots
using LaTeXStrings

function solutions(f::Function, a, b, label)
    println("Findig solutions...")
    sol, x = secant_method(f, a, b)
    @show sol

    println("Studying the convergence...")
    q, C = convergence(x, sol, skip=1)

    println("Plotting the data...")
    x_vals = range(a, b, length=200)
    p = plot_generic(x_vals, f.(x_vals), label=label)
    scatter_add!(p, [sol], [0.0], label="Solution")
    hline!([0.0], label="y=0", linestyle=:dash)

    return p
end

function main()
    f1 = x -> x^2 - exp(-x)
    a1 = -2
    b1 = 2

    f2 = x -> 2x - tan(x)
    a2 = 0.2
    b2 = 1.4

    f3 = x -> exp(x + 1) - 2 - x
    a3 = -2
    b3 = 2

    p1 = solutions(f1, a1, b1, L"x^2 - e^{-x}")
    p2 = solutions(f2, a2, b2, L"2x - \tan(x)")
    p3 = solutions(f3, a3, b3, L"e^{x + 1} - 2 - x")

    p = multi_plot(p1, p2, p3, layout=(1, 3), size=(1200, 600))
    save_plot(p, "secant-meghod-plot", "4-4")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
