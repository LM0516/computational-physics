using ComputationalPhysics
using Plots
using LaTeXStrings

function main()
    f = x -> x^(-2) - sin(x)
    df = x -> -2x^(-3) - cos(x)
    a = 0.5
    b = 10
    x_init = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    solutions = []

    for i in x_init
        sol, x = newton_method(f, df, a, b; x_init=i)

        println("Solution for x_init = $i: $sol")
        q_values, C_values = conv_rate_and_asymp_const(x, sol)
        @show q_values[end], C_values[end]

        push!(solutions, sol)
    end

    x_vals = range(a, b, length=200)
    p = plot_generic(x_vals, f.(x_vals), xlabel=L"x", ylabel=L"y", label=L"f(x) = x^{-2} - \sin(x)")
    scatter_add!(p, solutions, [0.0 for _ in solutions], label="Solutions")
    hline!([0.0], label="y=0", linestyle=:dash)
    save_plot(p, "multiple-newton-method-roots", "4-3")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
