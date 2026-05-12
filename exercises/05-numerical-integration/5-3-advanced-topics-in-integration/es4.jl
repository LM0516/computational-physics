using ComputationalPhysics
using Plots
using SpecialFunctions
using LaTeXStrings

function deq_solutions(f::Function, a::Real, b::Real, N::Int)
    # Transformation for finite interval [a, b] to real line for DE rule
    # We use the standard DE transformation x(t) = tanh(π/2 * sinh(t)) which maps (-inf, inf) to (-1, 1)
    # Then we map (-1, 1) to (a, b) linearly: y(x) = (b-a)/2 * x + (b+a)/2

    # Combined transformation y(t):
    y(t) = (b - a) / 2 * tanh(π / 2 * sinh(t)) + (b + a) / 2

    # Derivative dy/dt:
    # dy/dx = (b-a)/2
    # dx/dt = (π/2 * cosh(t)) / cosh(π/2 * sinh(t))^2
    dy_dt(t) = (b - a) / 2 * (π / 2 * cosh(t)) / (cosh(π / 2 * sinh(t)))^2

    g(t) = f(y(t)) * dy_dt(t)

    I = double_exponential_quadrature(g, N)
    return I
end

function plot_errors(glq_integrals::Vector{Float64}, deq_integrals::Vector{Float64}, exact, function_eq; fig_name="test")
    println("Calculating the errors...")
    glq_err = @. abs(glq_integrals - exact)
    deq_err = @. abs(deq_integrals - exact)

    println("Plotting the errors...")
    p = scatter_generic(1:length(glq_err), make_log_safe(glq_err), label="Gauss-Legendre integration", yscale=:log10)
    scatter_add!(p, 1:length(deq_err), make_log_safe(deq_err), label="Double exponential quadrature integration")
    xlabel!(p, "Number of nodes")
    ylabel!(p, "Errors")
    save_plot(p, "comparison-gl-deq-$fig_name", "5-3")
    display(p)
    readline()
end

function main()
    functions = [
        (
            f=x -> 1 / (sqrt(1 - x^2)),
            exact=pi / 2,
            label=L"\frac{1}{\sqrt{1 - x^2}}",
            a=0,
            b=1
        ),
        (
            f=x -> sqrt(x) * log(x),
            exact=-4 / 9,
            label=L"\sqrt{x} \log(x)",
            a=0,
            b=1
        ),
        (
            f=x -> sqrt(1 - x^2),
            exact=pi / 4,
            label=L"\sqrt{1 - x^2}",
            a=0,
            b=1
        ),
        (
            f=x -> (log(x))^2,
            exact=2,
            label=L"(\log(x))^2",
            a=0,
            b=1
        ),
        (
            f=x -> log(cos(x)),
            exact=-π / 2 * log(2),
            label=L"\log(\cos(x))",
            a=0,
            b=π / 2
        ),
        (
            f=x -> sqrt(tan(x)),
            exact=π / sqrt(2),
            label=L"\sqrt{\tan(x)}",
            a=0,
            b=π / 2
        )
    ]

    for (count, fun) in enumerate(functions)
        println("Function: $(fun.label)")
        i = 1
        glq_sol = Vector{Float64}(undef, 29)
        deq_sol = Vector{Float64}(undef, 29)
        for (i, n) in enumerate(4:2:60)
            N = n / 2
            glq_sol[i] = glq_integral(fun.f, fun.a, fun.b, Int(n)) # compute the solution using the library
            deq_sol[i] = deq_solutions(fun.f, fun.a, fun.b, Int(N)) # compute the solution with the correct transformation
        end
        # Comparing the integrals solutions
        @show glq_sol[end]
        @show deq_sol[end]
        @show fun.exact
        plot_errors(glq_sol, deq_sol, fun.exact, fun.label; fig_name=count)
        println("="^60)
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
