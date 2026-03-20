using ComputationalPhysics
using LinearAlgebra

function check(M, x, b)
    #== Checks the accuracy of the solution by computing M*x - b ==#
    c = M * x .- b
    if all(abs.(c) .< 1e-10)
        println("Correct! The solution is accurate with at least 1e-10 precision.")
    else
        max_error = maximum(abs.(c))
        println("Try again please. Maximum error: $max_error")
    end
end

function main()
    # === Forward substitution ===
    L1 = [-2 0 0; 1 -1 0; 3 2 1]
    L2 = [4 0 0 0; 1 -2 0 0; -1 4 4 0; 2 -5 5 1]
    b1 = [-4, 2, 1]
    b2 = [-4, 1, -3, 5]

    x1 = forward_substitution(L1, b1)
    x2 = forward_substitution(L2, b2)
    println("Solution for system 1: ", x1)
    println("Solution for system 2: ", x2)

    # === Backward substitution ===
    U1 = [3 1 0; 0 -1 -2; 0 0 3]
    U2 = [3 1 0 6; 0 -1 -2 7; 0 0 3 4; 0 0 0 5]
    c1 = [1, 1, 6]
    c2 = [4, 1, 1, 5]

    y1 = backward_substitution(U1, c1)
    y2 = backward_substitution(U2, c2)
    println("Solution for system 3: ", y1)
    println("Solution for system 4: ", y2)

    # === Random matrix ===
    n = 4
    L_random = LowerTriangular(rand(n, n))
    B_random = rand(n, 3)

    X = forward_substitution_multiple(Matrix(L_random), B_random)
    println("\nRandom L matrix:")
    display(L_random)
    println("\nRandom B matrix:")
    display(B_random)
    println("\nSolution X:")
    display(X)
    println()
    check(Matrix(L_random), X, B_random)

    # === Checks ===
    println("\n=== Verification ===")
    check(L1, x1, b1)
    check(L2, x2, b2)
    check(U1, y1, c1)
    check(U2, y2, c2)
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
