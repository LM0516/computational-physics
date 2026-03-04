include("../../modules/nonlinear_equatioins.jl")
using .NonlinearEquations
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
        q, C = convergence(x, sol, skip=1)

        push!(solutions, sol)
    end

    x_vals = range(a, b, length=200)
    p = plot(x_vals, f.(x_vals), label=L"f(x) = x^{-2} - \sin(x)")
    scatter!(solutions, [0.0 for _ in solutions], label="Solutions")
    hline!([0.0], label="y=0", linestyle=:dash)
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
