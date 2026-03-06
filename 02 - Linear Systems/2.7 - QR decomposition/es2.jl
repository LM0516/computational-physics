include("../../modules/linear_systemsV2.jl")
using .LinearSystemsV2
using LinearAlgebra

function solver(A::Matrix)
    println("QR-decomposition for: ", A)
    Q, R = qr_mgs(Matrix{Float64}(A))

    println()
    println("Q matrix: ", Q)
    println("R matrix: ", R)
    println()
    # `isapprox` (or its operator `≈`) handles floating-point errors
    A ≈ Q * R ? println("QR-decomposition correct") : println("Something's wrong")

    # Check the actual residual norm directly
    println("Absolute residual: ", norm(A - Q * R))
    println("Relative residual: ", norm(A - Q * R) / norm(A))
end

function main()
    A1 = [2 -1
        0 1
        -2 2]

    A2 = [2 3 4
        4 5 10
        4 8 2]

    A3 = [2 -1 4
        0 3 5
        7 2 -6
        1 -4 0]

    println("Modified Gram-Schmidt QR Decomposition Implementation")
    solver(A1)
    solver(A2)
    solver(A3)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
