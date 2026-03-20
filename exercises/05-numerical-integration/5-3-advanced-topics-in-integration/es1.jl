include("../../modules/numerical_integration.jl")
using .NumericalIntegration
using Plots
using SpecialFunctions

function solutions(f::Function, a::Real, b::Real, exact::Real, function_name)
    sol = zeros(19)
    alt = zeros(19)
    i = 1

    println("Calculating the solutions...")
    for n in 4:2:40
        sol[i] = fajer_rule(f, n, a, b)
        alt[i] = clenshaw_curtis_rule(f, n)
        i += 1
    end

    println("Calculating the errors...")
    err = @. log(abs(sol - exact))
    alt_err = @. log(abs(alt - exact))

    println("Plotting the errors...")
    p = scatter(err, label=function_name)
    scatter!(alt_err, label="Clenshaw-Curtis")
    xlabel!(p, "Number of nodes")
    ylabel!(p, "Error")
    display(p)
    readline()
end

function main()
    # exp(-4x)
    f1 = x -> exp(-4x)
    a1 = -1
    b1 = 1
    exact_sol1 = sinh(4) / 2
    solutions(f1, a1, b1, exact_sol1, L"\frac{\sinh(4)}{2}")

    # exp(-9x^2)
    f2 = x -> exp(-9x^2)
    a2 = -1
    b2 = 1
    exact_sol2 = sqrt(pi) * erf(3) / 3
    solutions(f2, a2, b2, exact_sol2, L"\exp(-9x^2)")

    # sech(x)
    f3 = x -> sech(x)
    a3 = -1
    b3 = 1
    exact_sol3 = 2 * atan(sinh(1))
    solutions(f3, a3, b3, exact_sol3, L"sech(x)")

    # 1 / (1 + 9x^2)
    f4 = x -> 1 / (1 + 9x^2)
    a4 = -1
    b4 = 1
    exact_sol4 = 2 / 3 * atan(3)
    solutions(f4, a4, b4, exact_sol4, L"\frac{1}{1 + 9x^2}")

    # x^2 sin(8x)
    f5_original = x -> x^2 * sin(8x)
    a5_original = π / 2
    b5_original = π
    # Coordinate transformation: y = (b_orig - a_orig)/2 * x + (b_orig + a_orig)/2
    # Here, x is the variable in [-1, 1]
    # y = (π - π/2)/2 * x + (π + π/2)/2 = (π/4) * x + (3π/4)
    # The integral ∫_a^b f(y) dy = ∫_-1^1 f(y(x)) * (b_orig - a_orig)/2 dx
    # So, the new function to be integrated over [-1, 1] is f(y(x)) * (b_orig - a_orig)/2
    f5 = x -> f5_original((b5_original - a5_original)/2 * x + (b5_original + a5_original)/2) * (b5_original - a5_original)/2
    a5 = -1
    b5 = 1
    exact_sol5 = -(3 * π^2) / 32
    solutions(f5, a5, b5, exact_sol5, L"x^2 \sin(8x)")
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
