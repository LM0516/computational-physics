using ComputationalPhysics
using Plots
using SpecialFunctions
using LaTeXStrings

function solutions(f::Function, exact::Real, function_name; fig_name::String="test")
    errors = Float64[]
    num_nodes_values = Vector{Float64}(undef, 29)
    errors = zeros(Float64, 29)
    i = 1

    for n in 4:2:60
        N = Int(n / 2)
        println("Calculating the integral for N = $N")
        int = double_exponential_quadrature(f, N)
        current_err = abs(int - exact)

        # Handle zero errors for logarithmic scale: replace with machine epsilon
        if current_err == 0.0
            errors[i] = eps(Float64)
        else
            errors[i] = current_err
        end
        num_nodes_values[i] = N

        println("Calculated value: $int")
        println("Exact value: $exact")
        println("="^40)
        i += 1
    end

    p = scatter_generic(num_nodes_values, errors, label=function_name, yaxis=:log, marker=:circle, legend=:bottomleft)
    xlabel!(p, "Number of nodes (N)")
    ylabel!(p, "Absolute Error")
    title!(p, "Error vs. Number of Nodes")
    save_plot(p, "double-exponential-quadrature-$fig_name", "5-3")
    display(p)
    readline()
end

function main()
    # Transformations and their derivatives for double exponential quadrature
    # Transformation for (0, Inf)
    x_0_inf(t) = exp(π / 2 * sinh(t))
    dx_0_inf_dt(t) = x_0_inf(t) * (π / 2 * cosh(t))

    # Alternative transformation for (0, Inf)
    x_0_inf_alt(t) = exp(t - exp(-t))
    dx_0_inf_alt_dt(t) = x_0_inf_alt(t) * (1 + exp(-t))

    # Transformation for (-Inf, Inf)
    x_inf_inf(t) = sinh(π / 2 * sinh(t))
    dx_inf_inf_dt(t) = cosh(π / 2 * sinh(t)) * (π / 2 * cosh(t))

    # Original integrand functions f(x)
    f1(x) = 1 / (1 + x^2 + x^4)
    f2(x) = exp(-x^2) * cos(x)
    f3(x) = (1 + x^2)^(-2 / 3)
    f4(x) = 1 / (1 + x^2)
    f5(x) = exp(-x) / sqrt(x)

    # Transformed integrand for double exponential quadrature
    g1(t) = f1(x_inf_inf(t)) * dx_inf_inf_dt(t)
    g2(t) = f2(x_inf_inf(t)) * dx_inf_inf_dt(t)
    g3(t) = f3(x_inf_inf(t)) * dx_inf_inf_dt(t)
    g4(t) = f4(x_0_inf(t)) * dx_0_inf_dt(t)
    g5(t) = f5(x_0_inf_alt(t)) * dx_0_inf_alt_dt(t)

    # Exact integrals
    exact_int1 = π / sqrt(3)
    exact_int2 = exp(-1 / 4) * sqrt(π)
    exact_int3 = sqrt(π) * gamma(1 / 6) / gamma(2 / 3)
    exact_int4 = π / 2
    exact_int5 = sqrt(π)

    # Solutions
    solutions(g1, exact_int1, L"\frac{1}{1 + x^2 + x^4}", fig_name="1")
    solutions(g2, exact_int2, L"\exp(-x^2) \cos(x)", fig_name="2")
    solutions(g3, exact_int3, L"(1 + x^2)^{-\frac{2}{3}}", fig_name="3")
    solutions(g4, exact_int4, L"\frac{1}{1 + x^2}", fig_name="4")
    solutions(g5, exact_int5, L"\frac{\exp(-x)}{\sqrt{x}}", fig_name="5")

end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
