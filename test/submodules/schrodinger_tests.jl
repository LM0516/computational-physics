using Test
using ComputationalPhysics
using LinearAlgebra
using SparseArrays

@testset "Schrodinger equation" begin
    x_grid = collect(range(-1.0, 1.0; length=11))
    dx = x_grid[2] - x_grid[1]
    H = t_independent_hamiltonian(x_grid, x -> 0.5 * x^2)
    @test issparse(H)
    @test Matrix(H) ≈ Matrix(H)' atol=1e-12
    @test size(H) == (11, 11)

    t_grid = collect(0.0:0.01:0.05)
    psi0 = ComplexF64.(exp.(-x_grid.^2))
    psi0 ./= sqrt(sum(abs2.(psi0)) * dx)
    snap = schrodinger_solver1D(x_grid, t_grid, x -> 0.0, psi0)
    norms = [sum(abs2.(psi)) * dx for psi in snap]
    @test maximum(abs.(norms .- 1.0)) < 2e-6

    psi_test = ComplexF64[1, 1, 0, 0, 2]
    snap_test = [psi_test, psi_test]
    x_test = [-2.0, -1.0, 0.0, 1.0, 2.0]
    R, T = compute_RT(1.0, 0.0, 1.0, snap_test, x_test, -0.5, 0.5, x -> 0.0)
    @test R ≈ 2.0 atol=1e-12
    @test T ≈ 4.0 atol=1e-12
end
