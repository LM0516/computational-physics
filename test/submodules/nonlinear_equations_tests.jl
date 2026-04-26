using Test
using ComputationalPhysics

@testset "Nonlinear equations" begin
    history_orders(history, root) = begin
        errors = abs.(history .- root)
        [
            log(errors[k + 1] / errors[k]) / log(errors[k] / errors[k - 1])
            for k in 2:length(errors)-1
            if errors[k - 1] > 0 && errors[k] > 0 && errors[k + 1] > 0
        ]
    end

    f(x) = x^2 - 2
    df(x) = 2x
    root_true = sqrt(2)

    root_bis = bisection(f, 0.0, 2.0; tol=1e-12)
    @test root_bis ≈ root_true atol=1e-10

    root_newton, hist_newton = newton_method(f, df, 0.0, 2.0; x_init=1.5, xtol=1e-14, ftol=1e-14)
    @test root_newton ≈ root_true atol=1e-12
    q_newton, c_newton = conv_rate_and_asymp_const(hist_newton, root_true)
    newton_order = last(sample_tail(history_orders(hist_newton, root_true)))
    @test last(sample_tail(q_newton)) ≈ 2.0 atol=0.15
    @test all(isfinite, sample_tail(c_newton))
    @test newton_order ≈ 2.0 atol=0.15

    root_secant, hist_secant = secant_method(f, 1.0, 2.0; xtol=1e-14, ftol=1e-14)
    @test root_secant ≈ root_true atol=1e-12
    secant_order = last(sample_tail(history_orders(hist_secant, root_true)))
    @test secant_order ≈ 1.618 atol=0.2

    root_iqi, hist_iqi = inverse_quadratic_interpolation(f, 0.0, 1.0, 2.0; xtol=1e-14, ftol=1e-14)
    @test root_iqi ≈ root_true atol=1e-12
    @test length(hist_iqi) <= length(hist_secant) + 2
    @test abs(hist_iqi[end] - root_true) < abs(hist_iqi[end - 1] - root_true)

    synthetic_root = 1.0
    q_true = 2.0
    c_true = 0.5
    eps = Float64[1e-1]
    for _ in 1:7
        push!(eps, c_true * eps[end]^q_true)
    end
    synthetic_hist = synthetic_root .+ eps
    q_est, c_est = redirect_stdout(devnull) do
        convergence(synthetic_hist, synthetic_root; skip=1)
    end
    @test q_est ≈ q_true atol=0.05
    @test c_est ≈ c_true atol=0.05

    theoretical_bisection_steps = ceil(Int, log2((2.0 - 0.0) / 1e-14))
    @test length(hist_newton) < length(hist_secant) < theoretical_bisection_steps

    print_method_table("Nonlinear Metrics", [
        ("Newton", newton_order, length(hist_newton)),
        ("Secant", secant_order, length(hist_secant)),
        ("IQI", missing, length(hist_iqi)),
        ("Bisection", missing, theoretical_bisection_steps),
    ]; headers=("Function", "Convergence rate", "Computational cost"))

    print_metric_table("Nonlinear Details", [
        ("Newton C tail", sample_tail(c_newton)),
        ("Synthetic q", q_est),
        ("Synthetic C", c_est),
    ])
end
