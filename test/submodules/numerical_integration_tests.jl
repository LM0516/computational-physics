using Test
using ComputationalPhysics

@testset "Numerical integration" begin
    exact = 2.0
    f(x) = sin(x)
    ns = [20, 40, 80, 160]

    trap_errors = [abs(composite_trapezoidal(0.0, π, f, n) - exact) for n in ns]
    simp_errors = [abs(composite_simpson(0.0, π, f, n) - exact) for n in ns]
    trap_order = empirical_order(ns, trap_errors)
    simp_order = empirical_order(ns, simp_errors)
    @test trap_order ≈ 2.0 atol=0.2
    @test simp_order ≈ 4.0 atol=0.2

    trap_m = 12
    counted_trap, trap_calls = counted_scalar_function(f)
    composite_trapezoidal(0.0, π, counted_trap, trap_m)
    @test trap_calls[] == 2 * (trap_m + 1)

    simpson_m = 12
    counted_simp, simp_calls = counted_scalar_function(f)
    composite_simpson(0.0, π, counted_simp, simpson_m)
    @test simp_calls[] == 2 * simpson_m + 1

    glq_poly_n = 3
    nodes, weights = glq(glq_poly_n)
    @test sort(nodes) ≈ [-sqrt(3 / 5), 0.0, sqrt(3 / 5)] atol=1e-10
    @test sum(weights) ≈ 2.0 atol=1e-12
    glq_exact = 0.0
    glq_poly_val = glq_integral(x -> x^5, -1.0, 1.0, glq_poly_n)
    glq_poly_error = abs(glq_poly_val - glq_exact)
    @test glq_poly_val ≈ glq_exact atol=1e-12

    glq_cost_n = 4
    counted_glq, glq_calls = counted_scalar_function(x -> x^2)
    glq_integral(counted_glq, -1.0, 1.0, glq_cost_n)
    @test glq_calls[] == glq_cost_n

    fejer_n = 8
    counted_fejer, fejer_calls = counted_scalar_function(x -> x^2)
    fejer_exact = 2 / 3
    fejer_val = fejer_rule(counted_fejer, fejer_n, -1.0, 1.0)
    fejer_error = abs(fejer_val - fejer_exact)
    @test fejer_val ≈ fejer_exact atol=1e-12
    @test fejer_calls[] == fejer_n

    cc_n = 8
    counted_cc, cc_calls = counted_scalar_function(x -> x^2)
    cc_exact = 2 / 3
    cc_val = clenshaw_curtis_rule(counted_cc, cc_n, -1.0, 1.0)
    cc_error = abs(cc_val - cc_exact)
    @test cc_val ≈ cc_exact atol=1e-12
    @test cc_calls[] == cc_n + 1

    de_n = 100
    counted_de, de_calls = counted_scalar_function(t -> exp(-t^2))
    de_val = double_exponential_quadrature(counted_de, de_n)
    de_exact = sqrt(π)
    de_error = abs(de_val - de_exact)
    @test de_val ≈ de_exact atol=1e-6
    expected_search_calls = let
        calls = 0
        for k in 1:ceil(Int, 10.0 / 0.2)
            t = k * 0.2
            calls += 1
            if abs(exp(-t^2)) < 1e-15
                calls += 1
                break
            end
        end
        calls
    end
    @test de_calls[] == 1 + 2 * de_n + expected_search_calls

    print_method_table("Quadrature Metrics", [
        ("Trapezoidal", trap_order, "$(trap_calls[]) evals (m=$trap_m)"),
        ("Simpson", simp_order, "$(simp_calls[]) evals (m=$simpson_m)"),
        ("GLQ", glq_poly_error, "$(glq_calls[]) evals (n=$glq_cost_n)"),
        ("Fejer", fejer_error, "$(fejer_calls[]) evals (n=$fejer_n)"),
        ("Clenshaw-Curtis", cc_error, "$(cc_calls[]) evals (n=$cc_n)"),
        ("Double exponential", de_error, "$(de_calls[]) evals (N=$de_n)"),
    ]; headers=("Function", "Convergence / accuracy", "Computational cost"))
end
