include("../../modules/linear_systems.jl")
using .LinearSystems
using LinearAlgebra

function PLU_solve(A, label::String)
    println("-"^60)
    println("=== ", label, " ===")
    println("-"^60)
    P, L, U = LUdecrp(Matrix{Float64}(A))

    println("L:")
    display(L)
    println("U:")
    display(U)
    println("Check:")
    display(P*A - L * U)

    println("Condizion number of $label: ", condition_number(Matrix{Float64}(A)))

    detA = det(P) * det(L) * det(U)
    println("Check on $label: ", detA - det(A))
end

function main()
    println("\n" * "="^60)
    println("=== Exercise 2 and 4 ===")
    println("="^60)
    # Es1
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

    PLU_solve(A1, "A1")
    PLU_solve(A2, "A2")
    PLU_solve(A3, "A3")

    #=# LU-Factorization-Row Pivoting=#
    #=P1, L1, U1 = LUdecrp(Matrix{Float64}(A1))=#
    #=P2, L2, U2 = LUdecrp(Matrix{Float64}(A2))=#
    #=P3, L3, U3 = LUdecrp(Matrix{Float64}(A3))=#
    #==#
    #=println("A1:")=#
    #=println("\nL, U:")=#
    #=display(L1)=#
    #=display(U1)=#
    #=println("Check:")=#
    #=display(P1*A1 - L1 * U1)=#
    #==#
    #=println("-"^60)=#
    #=println("A2:")=#
    #=println("\nL, U:")=#
    #=display(L2)=#
    #=display(U2)=#
    #=println("Check:")=#
    #=display(P2*A2 - L2 * U2)=#
    #==#
    #=println("-"^60)=#
    #=println("A3:")=#
    #=println("\nL, U:")=#
    #=display(L3)=#
    #=display(U3)=#
    #=println("Check:")=#
    #=display(P3*A3 - L3 * U3)=#
    #==#
    #=println("-"^60)=#
    #=println("=== Condition Numbers ===")=#
    #=println("-"^60)=#
    #=println("Condizion number A1: ", condition_number(Matrix{Float64}(A1)))=#
    #=println("Condizion number A2: ", condition_number(Matrix{Float64}(A2)))=#
    #=println("Condizion number A3: ", condition_number(Matrix{Float64}(A3)))=#

    println("\n" * "="^60)
    println("=== Exercise 3 ===")
    println("="^60)
    T(x, y) = [
        1 0 x
        0 1 y
        0 0 1
    ]

    R(θ) = [
        cos(θ) sin(θ) 0
        -sin(θ) cos(θ) 0
        0 0 1
    ]
    
    # Compute A and z
    A = T(3, -1) * R(π / 5) * T(-3, 1)
    z = [2, 2, 1]

    println("-"^60)
    println("=== Find b = A*z ===")
    println("-"^60)

    b = A * z
    println("\nb = A * z = ")
    display(b)
    println()

    println("\n" * "-"^60)
    println("=== Find LU factorization of A ===")
    println("-"^60)

    # LU-Factorization of A
    P, L, U = LUdecrp(Matrix{Float64}(A))

    println("\nL = ")
    display(L)
    println("\nU = ")
    display(U)

    println("\nVerification: P*A - L*U = ")
    display(P*A - L * U)

    println("\n" * "-"^60)
    println("=== Solve Ax = b using LU factorization ===")
    println("-"^60)

    # Solve L*y = b using forward substitution
    y = forwardsub(L, b)
    println("\nSolve L*y = b")
    println("y = ")
    display(y)

    # Solve U*x = y using backward substitution
    x = backsub(U, y)
    println("\nSolve U*x = y")
    println("x = ")
    display(x)

    # Compute x - z
    diff = x - z
    println("\nx - z = ")
    display(diff)
    println()

    # Verification
    println("\nVerification: A*x - b = ")
    display(A * x - b)

    #=println("\n" * "="^60)=#
    #=println("Exercise 4")=#
    #=println("="^60)=#
    #==#
    #=detA1 = det(P1) * det(L1) * det(U1)=#
    #=println("Check on A1: ", detA1 - det(A1))=#
    #==#
    #=detA2 = det(P2) * det(L2) * det(U2)=#
    #=println("Check on A2: ", detA2 - det(A2))=#
    #==#
    #==#
    #=detA3 = det(P3) * det(L3) * det(U3)=#
    #=println("Check on A3: ", detA3 - det(A3))=#
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
