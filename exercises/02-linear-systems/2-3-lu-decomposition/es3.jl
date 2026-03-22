using ComputationalPhysics

function main()
    # Define transformation matrices
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
    z = [2.0, 2.0, 1.0]

    println("="^60)
    println("Find b = A*z")
    println("="^60)

    b = A * z
    println("\nb = A * z = ")
    display(b)
    println()

    println("\n" * "="^60)
    println("Find LU factorization of A")
    println("="^60)

    # LU-Factorization of A
    L, U = lu_decomposition(Matrix{Float64}(A))

    println("\nL = ")
    display(L)
    println("\nU = ")
    display(U)

    println("\nVerification: A - L*U = ")
    display(A - L * U)

    println("\n" * "="^60)
    println("Solve Ax = b using LU factorization")
    println("="^60)

    # Solve L*y = b using forward substitution
    y = forward_substitution(L, b)
    println("\nSolve L*y = b")
    println("y = ")
    display(y)

    # Solve U*x = y using backward substitution
    x = backward_substitution(U, y)
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

    # Condition number
    cond_number = matrix_condition_number(A)
    println("\nCondition number: $cond_number")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
