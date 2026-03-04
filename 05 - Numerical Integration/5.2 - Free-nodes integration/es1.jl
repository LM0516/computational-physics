include("../../modules/numerical_integration.jl")
using .NumericalIntegration
using GaussQuadrature

function check_results(r::Vector{Float64}, w::Vector{Float64}, n::Int)
    println("Checking results for n = $n...")
    real_nodes, real_weights = legendre(n)
    sort!(r)
    sort!(real_nodes)
    sort!(w)
    sort!(real_weights)

    for k in 1:n
        if abs(r[k] - real_nodes[k]) > 1e-14 || abs(w[k] - real_weights[k]) > 1e-14
            println("Error at index $k")
        end
    end
    println("Results are correct.")

    println("Nodes: ", r)
    println("Weights: ", w)
end

function main()
    n = 100
    r, w = gauss_legendre_quadrature(n)
    check_results(r, w, n)
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
