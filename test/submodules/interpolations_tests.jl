using Test
using ComputationalPhysics

@testset "Interpolations" begin
    f(x) = x^3 - 2x + 1
    x_nodes = collect(range(-1.0, 1.0; length=5))
    @test barycentric_lagrange(0.3, x_nodes, f) ≈ f(0.3) atol=1e-12
    @test barycentric_lagrange(x_nodes[3], x_nodes, f) == f(x_nodes[3])

    cheb_nodes = cos.(range(0.0, π; length=5))
    @test barycentric_lagrange(0.1, collect(cheb_nodes), f; method=ComputationalPhysics.Chebyshev()) ≈ f(0.1) atol=1e-12

    observed = [1.0, 2.0, 3.0]
    expected = [1.0, 2.0, 3.0]
    @test chi_square(observed, expected) == 0.0
    @test chi_square([1.0, 1.0], [0.0, 2.0]) == 0.5
    @test p_value(0.0, 3) == 1.0

    chi2, dof, chi2r, p = redirect_stdout(devnull) do
        fit_goodness(observed, expected, [1.0, 2.0])
    end
    @test (chi2, dof, chi2r, p) == (0.0, 1, 0.0, 1.0)
end
