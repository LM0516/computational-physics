include("../../modules/numerical_integration.jl")
using .NumericalIntegration
using Plots
using LaTeXStrings

function solutions(a::Real, b::Real, f::Function, nodes, label, exact; rule::String="trapezoidal")
    s = Vector{Float64}(undef, length(nodes))

    # Compute solutions 
    if rule == "trapezoidal"
        for j in 1:length(nodes)
            s[j] = composite_trapezoidal(a, b, f, Int(nodes[j]))
        end
    elseif rule == "simpson"
        for j in 1:length(nodes)
            s[j] = composite_simpson(a, b, f, Int(nodes[j]))
        end
    else
        error("Invalid rule: $rule")
    end

    # Plot solutions
    p = scatter(nodes, s, label=label)
    hline!([exact], label="Exact solution: $(round(exact, digits=4))", linestyle=:dash)
    xlabel!("Nodes")
    ylabel!("Solution")

    return s, p
end

function errors(s::Vector{Float64}, exact::Float64)
    println("-"^40)

    relative_errors = @. log10(abs(s - exact) / abs(exact))
    absolute_errors = @. log10(abs(s - exact))

    # Plot relative errors
    relative_plots = plot(relative_errors, label="Relative errors")
    xlabel!("Nodes")
    ylabel!("Relative error")

    # Plot absolute errors
    absolute_plots = plot(absolute_errors, label="Absolute errors")
    xlabel!("Nodes")
    ylabel!("Absolute error")

    return relative_plots, absolute_plots
end

function main()
    k = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    nodes = @. 10 * 2^k
    rule = "trapezoidal"

    println("=== $rule rule for integration ===")
    println("-"^40)

    println("=== First integral: x log(1 + x) ===")
    f1 = x -> x * log(1 + x)
    a1 = 0.0
    b1 = 1.0
    exact1 = 1 / 4
    s1, p1 = solutions(a1, b1, f1, nodes, L"x log(1 + x)", exact1, rule=rule)
    println("Solution: $(s1[end])")
    rp1, ap1 = errors(s1, exact1)

    println("=== Second integral: x^2 tan^-1(x) ===")
    f2 = x -> x^2 * atan(x)
    a2 = 0.0
    b2 = 1.0
    exact2 = (π - 2 + 2log(2)) / 12
    s2, p2 = solutions(a2, b2, f2, nodes, L"x^2 tan^{-1}(x)", exact2, rule=rule)
    println("Solution: $(s2[end])")
    rp2, ap2 = errors(s2, exact2)

    println("=== Third integral: e^x cos(x) ===")
    f3 = x -> exp(x) * cos(x)
    a3 = 0.0
    b3 = π / 2
    exact3 = (exp(π / 2) - 1) / 2
    s3, p3 = solutions(a3, b3, f3, nodes, L"e^x cos(x)", exact3, rule=rule)
    println("Solution: $(s3[end])")
    rp3, ap3 = errors(s3, exact3)

    println("=== Fourth integral: (tan^-1 (√2 + x^2))/(1 + x^2)(√2 + x^2) ===")
    f4 = x -> atan(sqrt(2 + x^2)) / ((1 + x^2) * sqrt(2 + x^2))
    a4 = 0.0
    b4 = 1.0
    exact4 = 5π^2 / 96
    s4, p4 = solutions(a4, b4, f4, nodes, L"\frac{\tan^{-1} (\sqrt{2 + x^2})}{(1 + x^2)\sqrt{2 + x^2}}", exact4, rule=rule)
    println("Solution: $(s4[end])")
    rp4, ap4 = errors(s4, exact4)

    println("=== Fifth integral: √x log(x) ===")
    f5 = x -> sqrt(x) * log(x)
    a5 = 0.0
    b5 = 1.0
    exact5 = -4 / 9
    s5, p5 = solutions(a5, b5, f5, nodes, L"\sqrt{x} log(x)", exact5, rule=rule)
    println("Solution: $(s5[end])")
    rp5, ap5 = errors(s5, exact5)

    println("=== Sixth integral: √(1 - x^2) ===")
    f6 = x -> sqrt(1 - x^2)
    a6 = 0.0
    b6 = 1.0
    exact6 = π / 4
    s6, p6 = solutions(a6, b6, f6, nodes, L"\sqrt{1 - x^2}", exact6, rule=rule)
    println("Solution: $(s6[end])")
    rp6, ap6 = errors(s6, exact6)

    # Show the results
    p_results = plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), size=(1200, 800))
    display(p_results)
    readline()

    # Show the relative errors
    p_relative = plot(rp1, rp2, rp3, rp4, rp5, rp6, layout=(2, 3), size=(1200, 800))
    display(p_relative)
    readline()

    # Show the absolute errors
    p_absolute = plot(ap1, ap2, ap3, ap4, ap5, ap6, layout=(2, 3), size=(1200, 800))
    display(p_absolute)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

