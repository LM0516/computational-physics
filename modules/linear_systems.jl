module LinearSystems

using LinearAlgebra

export forwardsub, forwardcompl, backsub, backcompl, LUdec, LUdecrp, plufact, detplu, condition_number, least_squares, ∞norm

"""
Solves the equation Lx = b for x using forward substitution.
"""
function forwardsub(L::Matrix, b)
    len = length(b)
    x = zeros(Float64, len)
    for i in 1:len
        x[i] = (b[i] - dot(L[i, 1:i-1], x[1:i-1])) / L[i, i]
    end
    return x
end

"""
Solves the equation LX = B for X using forward substitution for multiple right-hand sides.
"""
function forwardcompl(L::Matrix, B::Matrix)
    n_row, n_col = size(B)
    X = zeros(Float64, n_row, n_col)
    for j in 1:n_col
        X[:, j] = forwardsub(L, B[:, j])
    end
    return X
end

"""
Solves the equation Ux = b for x using backward substitution.
"""
function backsub(U::Matrix, b)
    len = length(b)
    x = zeros(Float64, len)
    for i in len:-1:1
        if i < len
            x[i] = (b[i] - dot(U[i, i+1:len], x[i+1:len])) / U[i, i]
        else
            x[i] = b[i] / U[i, i]
        end
    end
    return x
end

"""
Solves the equation UX = B for X using backward substitution for multiple right-hand sides.
"""
function backcompl(U::Matrix, B::Matrix)
    n_row, n_col = size(B)
    X = zeros(Float64, n_row, n_col)
    for j in 1:n_col
        X[:, j] = backward_sub(U, B[:, j])
    end
    return X
end

"""
Compute an LU decomposition of a square Float64 matrix A without pivoting.
"""
function LUdec(A::Matrix{Float64})::Tuple{Matrix{Float64}, Matrix{Float64}}
    dim = size(A, 1)
    L = Matrix{Float64}(I, dim, dim)
    U = zeros(dim, dim)

    A_work = copy(A)

    for i in 1:dim-1
        U[i, :] = A_work[i, :]
        L[:, i] = A_work[:, i] / U[i, i]
        A_work -= L[:, i] * U[i, :]'
    end
    U[dim, dim] = A_work[dim, dim]
    return L, U
end

"""
Compute an LU decomposition of a square Float64 matrix A using partial (row) pivoting.
"""
function LUdecrp(A::Matrix{Float64})
    #= LU decomposition with row pivoting =#
    n = size(A, 1)
    L_prime, U = zeros(n, n), zeros(n, n)
    pivot_rows = zeros(Int, n)
    A_work = copy(A)
    # Perform the decomposition
    for i in 1:n
        # Select pivot row
        pivot_rows[i] = argmax(abs.(A_work[:, i]))
        # Swap rows in working matrix
        U[i, :] = A_work[pivot_rows[i], :]
        # Build L'
        L_prime[:, i] = A_work[:, i] / U[i, i]
        # Update working matrix
        A_work -= L_prime[:, i] * U[i, :]'
    end
    # Build permutation matrix P and final L
    P = Matrix(I, n, n)
    # Reorder rows of P according to pivot_rows
    P = P[pivot_rows, :]
    L = P * L_prime
    return P, L, U
end

# BUG: Sistemare algorimo, questo è sbagliato
"""
Compute an LU decomposition of a square Float64 matrix A using partial (row) pivoting.
"""
function plufact(A::Matrix{Float64})
    n = size(A, 1)
    L_prime, U = zeros(n, n), zeros(n, n)
    p = collect(1:n)
    A_work = copy(A)
    for i in 1:n
        pivot_index = argmax(abs.(A_work[i:n, i])) + i - 1
        p[i], p[pivot_index] = p[pivot_index], p[i]
        A_work[i, :], A_work[pivot_index, :] = A_work[pivot_index, :], A_work[i, :]
        U[i, :] = A_work[i, :]
        L_prime[:, i] = A_work[:, i] / U[i, i]
        A_work -= L_prime[:, i] * U[i, :]'
    end
    return L_prime, U, p
end

"""
Compute an LU decomposition of a square Float64 matrix A using partial (row) pivoting.
"""
function detplu(A)
    _, U, p = plufact(A)
    d = prod(diag(U))
    S = 0
    for i in 1:length(p)
        while p[i] != i
            p[p[i]], p[i] = p[i], p[p[i]]
            S += 1
        end
    end
    d *= (-1)^S
    return d
end

"""
Compute the matrix one norm (induced 1-norm), defined as the maximum absolute column sum of A.
"""
function onenorm(A::Matrix{Float64})
    return maximum(sum(abs, A, dims=1))
end

"""
Compute the matrix infinity norm (induced ∞-norm), defined as the maximum absolute row sum of A.
"""
function ∞norm(A::Matrix{Float64})
    return maximum(sum(abs, A, dims=2))
end

"""
Calculate the condition number using `onenorm` or `infinity norm`.
The smallest value a condition number can have is 1. In that case, the relative perturbation of 
the solution has the same size as that of the data.
"""
function condition_number(A::Matrix{Float64}, norm_type::Symbol=:one)
    if norm_type == :one
        nor_function = onenorm
    elseif norm_type == :infinity
        nor_function = ∞norm
    else
        error("Unsupported norm type. Use :one or :infinity.")
    end
    return nor_function(A) * nor_function(inv(A))
end

function least_squares(A::Matrix, b)
    N = A' * A
    z = A' * b

    L, U = LUdec(Matrix{Float64}(N))
    c = forwardsub(L, z)
    x = backsub(U, c)

    return x
end

end # module LinearSystems
