using ComputationalPhysics
using Plots
using SpecialFunctions
using LaTeXStrings

function solutions(f::Function, a::Real, b::Real, exact::Real, function_name; fig_name::String="test")
    fajer_solutions = Vector{Float64}(undef, 19)
    clenshaw_curtis_solutions = Vector{Float64}(undef, 19)
    glq_solutions = Vector{Float64}(undef, 19)

    println("Calculating the solutions...")
    for (i, n) in enumerate(4:2:40)
        fajer_solutions[i] = fejer_rule(f, n, a, b) # Fayer-rule integral solutions
        clenshaw_curtis_solutions[i] = clenshaw_curtis_rule(f, n) # Clenshaw-Curtis integral solutions
        glq_solutions[i] = glq_integral(f, a, b, n) # Gauss-Legendre quadrature integral
    end
    @show fajer_solutions[end]
    @show clenshaw_curtis_solutions[end]
    @show glq_solutions[end]

    println("Calculating the errors...")
    fajer_err = @. abs(fajer_solutions - exact)
    cc_err = @. abs(clenshaw_curtis_solutions - exact)
    glq_err = @. abs(glq_solutions - exact)

    println("Plotting the errors...")
    p = scatter_generic(1:length(fajer_err), make_log_safe(fajer_err), label="Fayer rule", yscale=:log10)
    # scatter_add!(p, 1:length(cc_err), make_log_safe(cc_err), label="Clenshaw-Curtis")
    scatter_add!(p, 1:length(glq_err), make_log_safe(glq_err), label="Gauss-Legendre")
    xlabel!(p, "Number of nodes")
    ylabel!(p, "Error")
    # title!(function_name)
    save_plot(p, "clenshaw-curtis-plot-$fig_name", "5-3")
    display(p)
    readline()
end

function main()
    # exp(-4x)
    f1 = x -> exp(-4x)
    a1 = -1
    b1 = 1
    exact_sol1 = sinh(4) / 2
    @show exact_sol1
    solutions(f1, a1, b1, exact_sol1, L"e^{-4x}", fig_name="1")

    # exp(-9x^2)
    f2 = x -> exp(-9x^2)
    a2 = -1
    b2 = 1
    exact_sol2 = sqrt(pi) * erf(3) / 3
    @show exact_sol2
    solutions(f2, a2, b2, exact_sol2, L"\exp(-9x^2)", fig_name="2")

    # sech(x)
    f3 = x -> sech(x)
    a3 = -1
    b3 = 1
    exact_sol3 = 2 * atan(sinh(1))
    @show exact_sol3
    solutions(f3, a3, b3, exact_sol3, "sech(x)", fig_name="3")

    # 1 / (1 + 9x^2)
    f4 = x -> 1 / (1 + 9x^2)
    a4 = -1
    b4 = 1
    exact_sol4 = 2 / 3 * atan(3)
    @show exact_sol4
    solutions(f4, a4, b4, exact_sol4, L"\frac{1}{1 + 9x^2}", fig_name="4")

    # x^2 sin(8x)
    f5 = x -> x^2 * sin(8x)
    a5 = π / 2
    b5 = π
    exact_sol5 = -(3 * π^2) / 32
    @show exact_sol5
    solutions(f5, a5, b5, exact_sol5, L"x^2 \sin(8x)", fig_name="5")
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
