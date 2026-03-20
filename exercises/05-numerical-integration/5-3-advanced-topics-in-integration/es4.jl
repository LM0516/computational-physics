include("../../modules/numerical_integration.jl")
using .NumericalIntegration
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

function glq_solutions(f::Function, a::Real, b::Real, n::Int)
    # Get nodes and weights for [-1, 1]
    nodes, weights = gauss_legendre_quadrature(n)

    # Map nodes from [-1, 1] to [a, b]
    # x = (b-a)/2 * ξ + (b+a)/2
    mapped_nodes = @. (b - a) / 2 * nodes + (b + a) / 2

    # Scale weights
    # dx = (b-a)/2 * dξ
    mapped_weights = weights .* (b - a) / 2

    return sum(mapped_weights .* f.(mapped_nodes))
end

function plot_errors(glq_integrals::Vector{Float64}, deq_integrals::Vector{Float64}, exact::Float64, function_eq)
    println("Calculating the errors...")
    glq_err = @. log(abs(glq_integrals - exact))
    deq_err = @. log(abs(deq_integrals - exact))

    println("Plotting the errors...")
    p = scatter(glq_err, label="Gauss-Legendre integration", title=function_eq)
    scatter!(deq_err, label="Double exponential quadrature integration")
    xlabel!(p, "Number of nodes")
    ylabel!(p, "Errors")
    display(p)
    readline()
end

function main()
    functions = [
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
            exact=eps(),
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
    # NOTE: Ho assegnato come limite inveriore eps() perchè altrimenti la funzione decadrebbe troppo lentamente,
    # sarebbe interessante valutare se la versione con eps() è meglio di quella con la safeguard oppure no.
    # La versione con gli estermi giusti da warning perchè decade troppo lentamente e quindi gli viene assegnato
    # un valore di t_m prestabilito, tm = 10, scelto in modo arbitrario.

    for f in functions
        println("Function: $(f.label)")
        i = 1
        glq_sol = zeros(Float64, 29)
        deq_sol = zeros(Float64, 29)
        for n in 4:2:60
            N = n / 2
            glq_sol[i] = glq_solutions(f.f, f.a, f.b, Int(n))
            deq_sol[i] = deq_solutions(f.f, f.a, f.b, Int(N))
            i += 1
        end
        plot_errors(glq_sol, deq_sol, f.exact, f.label)
        println("="^60)
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
