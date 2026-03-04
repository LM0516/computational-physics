include("../../modules/nonlinear_equatioins.jl")
using .NonlinearEquations

function main()
    f = (x -> x^3 - 7*x^2 + 14*x - 6)
    a = 0.0
    b = 1.0
    label = "x^3 - 7x^2 + 14x - 6"
    n = collect(1:2:20)
    try 
        x, r = bisection_plot(f, a, b, label) 
        for i in n
            plot_convergence_analysis(x, r, skip_initial=i)
        end
    catch e
        println("Error occurred: $e")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
