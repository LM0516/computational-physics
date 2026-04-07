using ComputationalPhysics
using GaussQuadrature

function check_results(r::Vector{Float64}, w::Vector{Float64}, n::Int)
    println("\nChecking results for n = $n with GaussQuadrature package...")
    real_nodes, real_weights = legendre(n)

    # Sort calculated nodes and weights
    p = sortperm(r)
    correct_nodes = r[p]
    correct_weights = w[p]
    # Sort real nodes and weights
    p_real = sortperm(real_nodes)
    sorted_real_nodes = real_nodes[p_real]
    sorted_real_weights = real_weights[p_real]

    correct = true

    for k in 1:n
        if abs(correct_nodes[k] - sorted_real_nodes[k]) > 1e-14 || abs(correct_weights[k] - sorted_real_weights[k]) > 1e-14
            @error "Error at index $k"
            correct = false
        end
    end

    correct ? println("Check succesfull") : nothing
end

function main()
    n = [2, 3, 4, 5]
    for i in n
        r, w = gauss_legendre_quadrature(i)
        check_results(r, w, i)

        println("Nodes: ", r)
        println("Weights: ", w)
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
