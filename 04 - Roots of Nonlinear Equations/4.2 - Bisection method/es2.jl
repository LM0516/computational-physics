include("../../modules/nonlinear_equatioins.jl")
using .NonlinearEquations

function main()
    # Main function
    f = (x -> x^3 - 7*x^2 + 14*x - 6)
    # Domain of the solution
    a = 0.0
    b = 1.0
    label = "x^3 - 7x^2 + 14x - 6"
    # How many initial values it has to skip before calculating the convergence rate
    n = collect(1:2:20)
    # Plot the solution of the equation  calculated with the bisection method.
    x, r = bisection_plot(f, a, b, label, save_dir="4-2") 
    plot_global_convergence(x, r, save_dir="4-2")
    for i in n
        # Plot of the convergence; it skips `n` parameters.
        plot_convergence_analysis(x, r, skip_initial=i, save_dir="4-2")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
