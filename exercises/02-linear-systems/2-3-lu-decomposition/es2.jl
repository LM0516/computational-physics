using ComputationalPhysics
using LinearAlgebra

function main()
    A1 = [
        2 3 4
        4 5 10
        4 8 2
    ]
    A2 = [
        6 -2 -4 4
        3 -3 -6 1
        -12 8 21 -8
        -6 0 -10 7
    ]
    A3 = [
        1 4 5 -5
        -1 0 -1 -5
        1 3 -1 2
        1 -1 5 -1
    ]
    L1, U1 = lu_decomposition(Matrix{Float64}(A1))
    L2, U2 = lu_decomposition(Matrix{Float64}(A2))
    L3, U3 = lu_decomposition(Matrix{Float64}(A3))

    println("="^60)
    println("LU-Factorization")
    println("="^60)

    println("L1 and U1")
    display(L1)
    display(U1)

    println("L2 and U2")
    display(L2)
    display(U2)

    println("L3 and U3")
    display(L3)
    display(U3)

    println("-"^60)

    println("Check A1:")
    display(A1 - L1 * U1)
    println("Check A2:")
    display(A2 - L2 * U2)
    println("Check A3:")
    display(A3 - L3 * U3)

    println("="^60)
    println("Determinant with LU-Decomposition")
    println("="^60)

    detA1 = det(L1) * det(U1)
    println("Check on A1: ", detA1 - det(A1))

    detA2 = det(L2) * det(U2)
    println("Check on A2: ", detA2 - det(A2))

    detA3 = det(L3) * det(U3)
    println("Check on A3: ", detA3 - det(A3))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
