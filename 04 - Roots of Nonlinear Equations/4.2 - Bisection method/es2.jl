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
    # Automatically check if there are any errors and prints them.
    try 
        # Plot the solution of the equation  calculated with the bisection method.
        x, r = bisection_plot(f, a, b, label) 
        for i in n
            # Plot of the convergence; it skips `n` parameters.
            plot_convergence_analysis(x, r, skip_initial=i)
        end
    catch e
        println("Error occurred: $e")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
