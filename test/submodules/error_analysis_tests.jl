using Test
using ComputationalPhysics
using Statistics

@testset "Error analysis" begin
    @test mclaurin_series(1.0, 14) ≈ exp(1) - 1 atol=1e-10

    x = Float64[1, 2, 3, 4, 5]
    @test var_double_pass(x) ≈ var(x; corrected=true)
    @test var_single_pass(x) ≈ var(x; corrected=true)

    x_shifted = 1.0e8 .+ Float64[1, 2, 3, 4, 5]
    stable_var = var(x_shifted; corrected=true)
    @test var_double_pass(x_shifted) ≈ stable_var rtol=1e-8
    @test abs(var_single_pass(x_shifted) - stable_var) > 1e-6

    @test kappa_f(2.0, x -> x^3) ≈ 3.0 atol=1e-12
end
