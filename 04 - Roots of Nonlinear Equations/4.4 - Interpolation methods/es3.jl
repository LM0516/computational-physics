include("../../modules/nonlinear_equatioins.jl")
using .NonlinearEquations
using Plots
using LaTeXStrings

function solutions(f, x0, x1, x2, a, b, label)
    println("\n--- Solving $label ---")
    sol, history = inverse_quadratic_interpolation(f, x0, x1, x2)
    println("Solution: $sol")

    try
        q, C = convergence(history, sol; skip=2)
        println("Order of convergence q ≈ $q")
    catch e
        println("Could not calculate convergence: $e")
    end

    # Saving the plot
    println("Plotting the data...")
    x_vals = range(a, b, length=200)
    p = plot(x_vals, f.(x_vals), label=label)
    scatter!([sol], [0.0], label="Solution")
    hline!([0.0], label="y=0", linestyle=:dash)

    return p
end

function main()
    # 1. x^2 = e^{-x} => f(x) = x^2 - e^{-x}
    f1 = x -> x^2 - exp(-x)
    a1 = -2
    b1 = 2
    p1 = solutions(f1, 0.0, 0.5, 1.0, a1, b1, L"x^2 - e^{-x}")

    # 2. 2x = tan(x) => f(x) = 2x - tan(x)
    f2 = x -> 2x - tan(x)
    a2 = -0.2
    b2 = 1.4
    p2 = solutions(f2, 1.0, 1.2, 1.4, a2, b2, L"2x - tan(x)")

    # 3. e^{x+1} = 2 + x => f(x) = e^{x+1} - 2 - x
    f3 = x -> exp(x + 1) - 2 - x
    a3 = -2
    b3 = 2
    p3 = solutions(f3, -1.2, -0.5, 0.5, a3, b3, L"e^{x+1} - 2 - x")

    println("Plotting the results...")
    p = plot(p1, p2, p3, layout=(1, 3), size=(1200, 600))
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
