using ComputationalPhysics

"""
Test the QR decomposition by checking:
1. A ≈ QR
2. Q^T Q ≈ I (orthonormality)

# Arguments
- `A`: Input matrix
- `tol`: Tolerance for verification (default: 1e-10)
"""
function test_qr(A::Matrix; tol=1e-10)
    println("Testing QR decomposition for matrix of size $(size(A))")
    println("=" ^ 60)
    
    # Compute QR decomposition
    Q, R = qr_mgs(A)
    
    # Test 1: A ≈ QR
    reconstruction_error = norm(A - Q * R)
    println("1. Reconstruction test: ||A - QR|| = $reconstruction_error")
    @assert reconstruction_error < tol "Reconstruction failed!"
    
    # Test 2: Q^T Q ≈ I
    QtQ = Q' * Q
    orthonormality_error = norm(QtQ - I)
    println("2. Orthonormality test: ||Q^T Q - I|| = $orthonormality_error")
    @assert orthonormality_error < tol "Orthonormality failed!"
    
    # Test 3: R is upper triangular
    is_upper = all(R[i, j] ≈ 0 for i in 1:size(R, 1) for j in 1:size(R, 2) if i > j)
    println("3. Upper triangular test: R is upper triangular = $is_upper")
    @assert is_upper "R is not upper triangular!"
    
    println("\nAll tests passed! ✓")
    println("=" ^ 60)
    
    return Q, R
end
