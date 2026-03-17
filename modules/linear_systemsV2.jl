module LinearSystemsV2

using LinearAlgebra

export forward_substitution, forward_substitution_multiple,
    backward_substitution, backward_substitution_multiple,
    lu_decomposition, lu_decomposition_pivoting,
    plu_factorization, determinant_plu,
    matrix_condition_number, solve_least_squares,
    qr_mgs, test_qr, qr_mgs_augmented,
    pure_qr_algorithm, matrix_infinity_norm

"""
Solve the lower triangular system Lx = b using forward substitution.
"""
function forward_substitution(L::Matrix, b::Vector)
    n = length(b)
    x = zeros(Float64, n)
    for i in 1:n
        x[i] = (b[i] - dot(L[i, 1:i-1], x[1:i-1])) / L[i, i]
    end
    return x
end

"""
Solve the lower triangular system LX = B for multiple right-hand sides using forward substitution.
"""
function forward_substitution_multiple(L::Matrix, B::Matrix)
    n_row, n_col = size(B)
    X = zeros(Float64, n_row, n_col)
    for j in 1:n_col
        X[:, j] = forward_substitution(L, B[:, j])
    end
    return X
end

"""
Solve the upper triangular system Ux = b using backward substitution.
"""
function backward_substitution(U::Matrix, b::Vector)
    n = length(b)
    x = zeros(Float64, n)
    for i in n:-1:1
        if i < n
            x[i] = (b[i] - dot(U[i, i+1:n], x[i+1:n])) / U[i, i]
        else
            x[i] = b[i] / U[i, i]
        end
    end
    return x
end

"""
Solve the upper triangular system UX = B for multiple right-hand sides using backward substitution.
"""
function backward_substitution_multiple(U::Matrix, B::Matrix)
    n_row, n_col = size(B)
    X = zeros(Float64, n_row, n_col)
    for j in 1:n_col
        X[:, j] = backward_substitution(U, B[:, j])
    end
    return X
end

"""
Compute an LU decomposition of a square matrix A without pivoting: A = LU.
Uses Gaussian elimination without row pivoting to decompose A into the product of 
a lower triangular matrix L and an upper triangular matrix U.
"""
function lu_decomposition(A::Matrix{Float64})
    n = size(A, 1)
    L = diagm(ones(n))
    U = zeros(n, n)
    A_work = copy(A)

    for i in 1:n-1
        U[i, :] = A_work[i, :]
        L[:, i] = A_work[:, i] / U[i, i]
        A_work -= L[:, i] * U[i, :]'
    end
    U[n, n] = A_work[n, n]

    return L, U
end

"""
Compute an LU decomposition of a square matrix A using partial row pivoting: PA = LU.
"""
function lu_decomposition_pivoting(A::Matrix{Float64})
    n = size(A, 1)
    L_prime = zeros(n, n)
    U = zeros(n, n)
    pivot_rows = zeros(Int, n)
    A_work = copy(A)

    # Perform decomposition with pivoting
    for i in 1:n
        pivot_rows[i] = argmax(abs.(A_work[:, i]))
        U[i, :] = A_work[pivot_rows[i], :]
        L_prime[:, i] = A_work[:, i] / U[i, i]
        A_work -= L_prime[:, i] * U[i, :]'
    end

    # Build permutation matrix P
    P = Matrix{Float64}(I, n, n)
    P = P[pivot_rows, :]

    # Compute final L with unit diagonal
    L = P * L_prime

    return P, L, U
end

"""
Compute a PLU factorization of a square matrix A using partial row pivoting.
"""
function plu_factorization(A::Matrix{Float64})
    n = size(A, 1)
    L_prime = zeros(n, n)
    U = zeros(n, n)
    p = collect(1:n)
    A_work = copy(A)

    for i in 1:n
        # Find pivot in column i from row i onward
        pivot_index = argmax(abs.(A_work[i:n, i])) + i - 1

        # Swap in permutation vector
        p[i], p[pivot_index] = p[pivot_index], p[i]

        # Swap rows in working matrix
        A_work[i, :], A_work[pivot_index, :] = A_work[pivot_index, :], A_work[i, :]

        # Store in U and L_prime
        U[i, :] = A_work[i, :]
        L_prime[:, i] = A_work[:, i] / U[i, i]

        # Update working matrix
        A_work -= L_prime[:, i] * U[i, :]'
    end

    return L_prime, U, p
end

"""
Compute the determinant of a square matrix using PLU factorization.
"""
function determinant_plu(A::Matrix{Float64})
    _, U, p = plu_factorization(A)

    # Determinant is product of diagonal of U
    det = prod(diag(U))

    # Count number of swaps to adjust sign
    num_swaps = 0
    for i in 1:length(p)
        while p[i] != i
            # Swap to move element to correct position
            p[p[i]], p[i] = p[i], p[p[i]]
            num_swaps += 1
        end
    end

    # Adjust sign based on number of swaps
    det *= (-1)^num_swaps

    return det
end

"""
Compute the matrix 1-norm (maximum absolute column sum).
"""
function matrix_one_norm(A::Matrix{Float64})
    return maximum(sum(abs, A, dims=1))
end

"""
Compute the matrix ∞-norm (maximum absolute row sum).
"""
function matrix_infinity_norm(A::Matrix{Float64})
    return maximum(sum(abs, A, dims=2))
end

"""
Compute the condition number of a square invertible matrix A.
"""
function matrix_condition_number(A::Matrix{Float64}, norm_type::Symbol=:one)
    if norm_type == :one
        norm_func = matrix_one_norm
    elseif norm_type == :infinity
        norm_func = matrix_infinity_norm
    else
        error("Unsupported norm type. Use :one or :infinity.")
    end

    return norm_func(A) * norm_func(inv(A))
end

"""
Solve the linear least squares problem min ||Ax - b||₂ using normal equations.
"""
function solve_least_squares(A::Matrix, b::Vector)
    # Form normal equations: A'Ax = A'b
    N = A' * A  # Gram matrix (n×n)
    z = A' * b  # Projected right-hand side (length n)

    # Solve Nx = z using LU decomposition
    L, U = lu_decomposition(Matrix{Float64}(N))
    c = forward_substitution(L, z)
    x = backward_substitution(U, c)

    return x
end

"""
Compute the thin QR decomposition of matrix A using the Modified Gram-Schmidt algorithm.
Returns Q̂ (m×n ONC matrix) and R̂ (n×n upper triangular matrix) such that A = Q̂R̂.
"""
function qr_mgs(A)
    m, n = size(A)

    # Check dimensions
    if m < n
        error("Matrix must have m ≥ n for thin QR decomposition")
    end

    # Initialize Q and R
    Q = zeros(m, n)
    R = zeros(n, n)
    V = copy(A)

    # Modified Gram-Schmidt procedure
    for j in 1:n
        # Start with v_j^(1) = a_j (already in V[:, j])

        # Apply projectors P_⊥q_i for i = 1, ..., j-1
        for i in 1:j-1
            # Compute r_ij = q_i^T * v_j^(i)
            R[i, j] = dot(Q[:, i], V[:, j])

            # Update v_j^(i+1) = v_j^(i) - q_i * (q_i^T * v_j^(i))
            V[:, j] -= R[i, j] * Q[:, i]
        end

        # Compute r_jj = ||v_j||
        R[j, j] = norm(V[:, j])

        # Normalize to get q_j = v_j / ||v_j||
        Q[:, j] = V[:, j] / R[j, j]
    end

    return Q, R
end


# Q-less QR for augmented matrix
function qr_mgs_augmented(A::Matrix{T}, b::Vector{T}) where T<:Real
    m, n = size(A)
    A_aug = hcat(A, b)
    Q_aug, R_aug = qr_mgs(A_aug)

    R = R_aug[1:n, 1:n]
    z = R_aug[1:n, n+1]

    return R, z
end

# Pure QR Algorithm
function pure_qr_algorithm(A; max_iter=1000, tol=1e-10)
    n = size(A, 1)
    A_k = copy(A)
    Q_total = Matrix{Float64}(I, n, n)

    for k in 1:max_iter
        # QR decomposition of A^(k-1)
        Q_k, R_k = qr_mgs(A_k)

        # Recombine in reverse order: A^(k) = R^(k) Q^(k)
        A_k_new = R_k * Q_k

        # Accumulate eigenvector matrix
        Q_total = Q_total * Q_k

        # Check convergence (off-diagonal elements should go to zero)
        off_diag_norm = 0.0
        for i in 1:n
            for j in 1:n
                if i != j
                    off_diag_norm += abs(A_k_new[i, j])
                end
            end
        end

        if off_diag_norm < tol
            println("Converged after $k iterations")
            return diag(A_k_new), Q_total, k
        end

        A_k = A_k_new
    end

    println("Warning: Did not converge within $max_iter iterations")
    return diag(A_k), Q_total, max_iter
end

end # module LinearSystems
