using ComputationalPhysics
using Plots
using LaTeXStrings

function main()
    functions = [
        (
            f=(t, u) -> -2t * u,
            a=0,
            b=2,
            u0=2,
            g=t -> 2 * exp(-t^2),
            title=L"u' = -2t u"
        ),
        (
            f=(t, u) -> u + t,
            a=0,
            b=1,
            u0=2,
            g=t -> -1 - t + 3 * exp(t),
            title=L"u' = u + t"
        ),
        (
            f=(t, u) -> t^2 / ((1 + t^3) * u),
            a=0,
            b=3,
            u0=1,
            g=t -> sqrt(1 + (2 / 3) * log(1 + t^3)),
            title=L"u' = \frac{t^2}{(1 + t^3) u}"
        )
    ]
    for (fig_name, f) in enumerate(functions)
        n = 320
        t, u = euler_method(f.f, f.a, f.b, n, f.u0)
        p1 = plot_generic(t, u, label="Numerical (Euler)", xlabel="Time (t)", ylabel="u(t)")
        plot_add!(p1, t, f.g.(t), label="Analytical (Exact)", linestyle=:dash)
        save_plot(p1, "euler-method-$fig_name", "6-2")
        display(p1)
        readline()

        ks = 2:10
        ns = [10 * 2^k for k in ks]
        errors = zeros(length(ns))
        final_errors = zeros(length(ns))
        for (i, n) in enumerate(ns)
            t, u = euler_method(f.f, f.a, f.b, n, f.u0)
            errors[i] = maximum(abs, u .- f.g.(t))
            final_errors[i] = abs(u[end] - f.g.(t[end]))
        end

        # First-order reference line: C / n, scaled to pass through the first data point
        ref_line = errors[1] .* (ns[1] ./ ns)

        # Plot on log-log scale to see convergence
        p2 = scatter_generic(ns, errors, xaxis=:log10, yaxis=:log10,
            xlabel="Number of steps (n)", ylabel="Error",
            label="$(f.title) Max Error")
        scatter_add!(p2, ns, final_errors, label="Final errors")
        plot_add!(p2, ns, ref_line, label=L"\mathcal{O}(n^{-1})", linestyle=:dash)
        save_plot(p2, "euler-method-error-$fig_name", "6-2")
        display(p2)
        readline()
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
