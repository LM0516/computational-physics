using ComputationalPhysics
using Printf

function main()
    # Define the matrix
    A_eigen = [5.0  1.0  0.0  0.0
               1.0  4.0  1.0  0.0
               0.0  1.0  3.0  1.0
               0.0  0.0  1.0  2.0]
    
    println("Original matrix A:")
    display(A_eigen)
    
    # Run the Pure QR algorithm to obtain eigenvalues and eigenvectors 
    eigenvalues_qr, eigenvectors_qr, n_iter = pure_qr_algorithm(A_eigen)
    
    println()
    println("Results from Pure QR Algorithm:")
    println("-"^80)
    println("Eigenvalues (sorted in descending order):")
    sorted_idx = sortperm(eigenvalues_qr, rev=true)
    for (i, idx) in enumerate(sorted_idx)
        @printf("  λ%d = %.15f\n", i, eigenvalues_qr[idx])
    end
    println()
    
    # Verification using Julia's built-in eigen function
    println("=== Verification ===")
    println("-"^80)
    
    eigen_result = eigen(A_eigen)
    eigenvalues_true = sort(eigen_result.values, rev=true)
    eigenvectors_true = eigen_result.vectors
    
    println("Built-in results")
    for (i, λ) in enumerate(eigenvalues_true)
        @printf("  λ%d = %.15f\n", i, λ)
    end
    println()
    
    println("Error in eigenvalues:")
    eigenvalues_qr_sorted = sort(eigenvalues_qr, rev=true)
    for i in 1:length(eigenvalues_qr_sorted)
        error = abs(eigenvalues_qr_sorted[i] - eigenvalues_true[i])
        @printf("  |λ%d_QR - λ%d_true| = %.2e\n", i, i, error)
    end
    println()
    
    # Check that we found true eigenvectors: A*v = λ*v
    println("VERIFICATION: A*v = λ*v for each eigenpair")
    println("-"^80)
    
    for i in 1:size(eigenvectors_qr, 2)
        idx = sorted_idx[i]
        v = eigenvectors_qr[:, idx]
        λ = eigenvalues_qr[idx]
        
        Av = A_eigen * v
        λv = λ * v
        
        error = norm(Av - λv)
        @printf("  Eigenvector %d: ||A*v - λ*v|| = %.2e: ", i, error)
        if error < 1e-8
            println("correct")
        else
            println("wrong")
        end
    end
    
    println()
    println("="^80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
