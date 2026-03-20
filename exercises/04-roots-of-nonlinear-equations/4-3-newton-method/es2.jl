using ComputationalPhysics
using Plots
using LaTeXStrings

function solutions(f::Function, df::Function, a, b, f_eqation, save_dir, plot_name)
    println("Findig solutions...")
    r, x = newton_method(f, df, a, b)
    q, C = convergence(x, r, skip=1)
    @show r

    println("Findig global convergence...")

    # p_global_convergence = plot_global_convergence(x, r, save_dir=save_dir, plot_name=plot_name)

    println("Solution found! Plotting the data...")

    x_vals = range(a, b, length=200)
    p = plot_generic(x_vals, f.(x_vals), xlabel=L"x", ylabel=L"y", label=f_eqation)
    scatter_add!(p, [r], [0.0], label=round(r; digits=4))
    hline!([0.0], label=L"y=0", linestyle=:dash)

    return p
end

function main()
    dir_name = "4-3"
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
    p1 = solutions(f1, df1, a1, b1, L"y = x^2 - e^{-x}", dir_name, "global-conv-1")
    p2 = solutions(f2, df2, a2, b2, L"2x - \tan(x)", dir_name, "global-conv-2")
    p3 = solutions(f3, df3, a3, b3, L"e^{x + 1} - 2 - x", dir_name, "global-conv-3")

    p = multi_plot(p1, p2, p3, layout=(1, 3), size=(1200, 600))

    save_plot(p, "plot-solution-newton-method", dir_name)
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
