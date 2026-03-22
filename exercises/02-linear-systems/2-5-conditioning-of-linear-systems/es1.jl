using ComputationalPhysics
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

    println("\n=== LU Decomposition without row pivoting for ε = ", ε, " and M = ", A, " ===")

    L, U = lu_decomposition(A)

    println("A = ")
    display(A)
    println("b = ")
    display(b)

    # Solve Ly = b (forward substitution)
    y = L \ b
    # Solve Ux = y (backward substitution)
    x_approx = U \ y

    println("\n=== Precision and residual error ===")

    # NOTE: The notation `A \ b` uses partioal pivoting to see a true comparison I have to 
    # compute `U \ (L \ b)` for a true comparison.

    println("Precision of x: ", norm(x - A \ b))
    println("Residual error: ", norm(b - A * (A \ b)))
    println("Precision of x (no pivoting): ", norm(x - x_approx))
    println("Residual error (no pivoting): ", norm(b - A * x_approx))

    println("\n=== LU Factorization ===")
    println("L:")
    display(L)
    println("U:")
    display(U)
    checkALU(A, L, U)


    println("\n=== Condition number for ε = ", ε, " ===")
    println("Condition number of A:", matrix_condition_number(A))
    println("Condition number of L: ", matrix_condition_number(L))
    println("Condition number of U: ", matrix_condition_number(U))

    println("\n=== Condition number check with LinearAlgebra ===")
    println("Condition number of A (LinearAlgebra): ", abs(cond(A, 1) - matrix_condition_number(A)) < 1e-10 ? "Match" : "Mismatch")
    println("Condition number of L (LinearAlgebra): ", abs(cond(L, 1) - matrix_condition_number(L)) < 1e-10 ? "Match" : "Mismatch")
    println("Condition number of U (LinearAlgebra): ", abs(cond(U, 1) - matrix_condition_number(U)) < 1e-10 ? "Match" : "Mismatch")

    println("\n=== Checking if L and U are ill-conditioned ===")
    println("L is ill-conditioned: ", matrix_condition_number(L) > 1e10 ? "Yes" : "No")
    println("U is ill-conditioned: ", matrix_condition_number(U) > 1e10 ? "Yes" : "No")
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
        println("")
        println("="^60)
        conditioning(PA, ε)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
