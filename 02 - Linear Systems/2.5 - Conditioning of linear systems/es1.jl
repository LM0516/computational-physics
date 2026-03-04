include("../../modules/linear_systems.jl")
using .LinearSystems
using LinearAlgebra

function checkALU(A, L, U)
    println("A - LU = ")
    display(A - L * U)
    if norm(A - L * U) < 1e-10
        println("LU factorization is correct within tolerance.")
    else
        println("LU factorization is NOT correct.")
    end
end

function conditioning(A::Matrix, ε::Float64)
    x = [1, 1]

    b = A * x

    println("\n=== LU Decomposition without row pivoting for ε = ", ε, "and M = ", A, " ===")
    println("A = ")
    display(A)
    println("b = ")
    display(b)

    println("\n=== Precision and residual error ===")
    println("Precision of x: ", norm(x - A \ b))
    println("Residual error: ", norm(b - A * (A \ b)))

    println("\n=== LU Factorization ===")
    L, U = LUdec(A)
    checkALU(A, L, U)

    println("\n=== Condition number for ε = ", ε, " ===")
    println("Condition number of A:", condition_number(A))
    println("Condition number of L: ", condition_number(L))
    println("Condition number of U: ", condition_number(U))

    println("\n=== Condition number check with LinearAlgebra ===")
    println("Condition number of A (LinearAlgebra): ", abs(cond(A, 1) - condition_number(A)) < 1e-10 ? "Match" : "Mismatch")
    println("Condition number of L (LinearAlgebra): ", abs(cond(L, 1) - condition_number(L)) < 1e-10 ? "Match" : "Mismatch")
    println("Condition number of U (LinearAlgebra): ", abs(cond(U, 1) - condition_number(U)) < 1e-10 ? "Match" : "Mismatch")

    println("\n=== Checking if L and U are ill-conditioned ===")
    println("L is ill-conditioned: ", condition_number(L) > 1e10 ? "Yes" : "No")
    println("U is ill-conditioned: ", condition_number(U) > 1e10 ? "Yes" : "No")

    println("")
    println("-"^60)
end

function main()
    ϵ = [1e-12, 1e-20]

    for ε in ϵ
        A = [
            -ε 1
            1 -1
        ]
        PA = [
            1 -1
            -ε 1
        ]
        conditioning(A, ε)
        println("="^60)
        conditioning(PA, ε)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
