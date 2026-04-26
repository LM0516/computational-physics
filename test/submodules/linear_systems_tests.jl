using Test
using ComputationalPhysics
using LinearAlgebra

@testset "Linear systems" begin
    L = [2.0 0.0 0.0; -1.0 3.0 0.0; 4.0 2.0 1.0]
    b = [2.0, 5.0, 6.0]
    x_true = L \ b
    @test forward_substitution(L, b) ≈ x_true atol=1e-12

    B = hcat(b, 2b)
    @test forward_substitution_multiple(L, B) ≈ L \ B atol=1e-12

    U = [2.0 -1.0 4.0; 0.0 3.0 2.0; 0.0 0.0 1.0]
    y = [5.0, 4.0, 3.0]
    @test backward_substitution(U, y) ≈ U \ y atol=1e-12
    Y = hcat(y, 2y)
    @test backward_substitution_multiple(U, Y) ≈ U \ Y atol=1e-12

    A = [4.0 3.0; 6.0 3.0]
    L_lu, U_lu = lu_decomposition(A)
    @test L_lu * U_lu ≈ A atol=1e-10

    P, Lp, Up = lu_decomposition_pivoting(A)
    @test P * A ≈ Lp * Up atol=1e-10

    Lprime, Uprime, p = plu_factorization(A)
    @test A[p, :] ≈ Lprime * Uprime atol=1e-10
    @test determinant_plu(A) ≈ det(A) atol=1e-10

    M = [1.0 -2.0; 3.0 4.0]
    @test ComputationalPhysics.matrix_one_norm(M) == 6.0
    @test matrix_infinity_norm(M) == 7.0
    @test matrix_condition_number(M, :one) ≈ opnorm(M, 1) * opnorm(inv(M), 1) atol=1e-10
    @test_throws ErrorException matrix_condition_number(M, :two)

    A_ls = [1.0 1.0; 1.0 2.0; 1.0 3.0]
    b_ls = [1.0, 2.0, 2.0]
    @test solve_least_squares(A_ls, b_ls) ≈ A_ls \ b_ls atol=1e-10

    Q, R = qr_mgs(A_ls)
    @test Q' * Q ≈ Matrix(I, 2, 2) atol=1e-12
    @test Q * R ≈ A_ls atol=1e-12

    R_aug, z_aug = qr_mgs_augmented(A_ls, b_ls)
    @test R_aug ≈ R atol=1e-12
    @test backward_substitution(R_aug, z_aug) ≈ A_ls \ b_ls atol=1e-10

    eigvals_qr, eigvecs_qr, iterations = redirect_stdout(devnull) do
        pure_qr_algorithm([2.0 1.0; 1.0 2.0])
    end
    @test sort(eigvals_qr) ≈ [1.0, 3.0] atol=1e-8
    @test eigvecs_qr' * eigvecs_qr ≈ Matrix(I, 2, 2) atol=1e-8
    @test iterations < 1000
end
