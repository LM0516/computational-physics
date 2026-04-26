using Test
using ComputationalPhysics

@testset "ODE solvers" begin
    f(t, u) = u
    exact = exp(1)
    ns = [10, 20, 40, 80]

    euler_errors = Float64[]
    ie2_errors = Float64[]
    rk4_errors = Float64[]

    for n in ns
        _, u_euler = euler_method(f, 0.0, 1.0, n, 1.0)
        _, u_ie2 = ie2(f, 0.0, 1.0, n, 1.0)
        _, u_rk4 = rk4(f, 0.0, 1.0, n, 1.0)
        push!(euler_errors, abs(u_euler[end] - exact))
        push!(ie2_errors, abs(u_ie2[end] - exact))
        push!(rk4_errors, abs(u_rk4[end] - exact))
    end

    euler_order = empirical_order(ns, euler_errors)
    ie2_order = empirical_order(ns, ie2_errors)
    rk4_order = empirical_order(ns, rk4_errors)
    @test euler_order ≈ 1.0 atol=0.1
    @test ie2_order ≈ 2.0 atol=0.1
    @test rk4_order ≈ 4.0 atol=0.2

    cost_n = 25

    counted_euler, euler_calls = counted_scalar_function(f)
    euler_method(counted_euler, 0.0, 1.0, cost_n, 1.0)
    @test euler_calls[] == cost_n

    counted_ie2, ie2_calls = counted_scalar_function(f)
    ie2(counted_ie2, 0.0, 1.0, cost_n, 1.0)
    @test ie2_calls[] == 2 * cost_n

    counted_rk4, rk4_calls = counted_scalar_function(f)
    rk4(counted_rk4, 0.0, 1.0, cost_n, 1.0)
    @test rk4_calls[] == 4 * cost_n

    g(t, u) = [u[2], -u[1]]
    _, u_vec = rk4(g, 0.0, π / 2, 100, [1.0, 0.0])
    @test u_vec[end][1] ≈ 0.0 atol=1e-6
    @test u_vec[end][2] ≈ -1.0 atol=1e-6

    print_method_table("ODE Metrics", [
        ("Euler", euler_order, "$(euler_calls[]) evals (n=$cost_n)"),
        ("IE2", ie2_order, "$(ie2_calls[]) evals (n=$cost_n)"),
        ("RK4", rk4_order, "$(rk4_calls[]) evals (n=$cost_n)"),
    ]; headers=("Function", "Convergence rate", "Computational cost"))
end
