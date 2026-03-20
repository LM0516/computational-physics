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
            title=L"u' = \sqrt{\frac{t^2}{(1 + t^3) \cdot u}}"
        )
    ]
    fig_name = 1
    for f in functions
        n = 320
        t, u = euler_method(f.f, f.a, f.b, n, f.u0)
        p = plot_generic(t, u, label="Numerical (Euler)", xlabel="Time (t)", ylabel="u(t)", title=f.title)
        plot_add!(p, t, f.g.(t), label=L"Analytical (Exact)", linestyle=:dash)
        save_plot(p, "euler-method-$fig_name", "6-2")
        display(p)
        readline()
        fig_name+=1
    end

    fig_name = 1
    for f in functions
        ks = 2:10
        ns = [10 * 2^k for k in ks]
        errors = zeros(length(ns))
        for (i, n) in enumerate(ns)
            t, u = euler_method(f.f, f.a, f.b, n, f.u0)
            errors[i] = maximum(abs, u .- f.g.(t))
        end
        
        # Plot on log-log scale to see convergence
        p = plot_generic(ns, errors, xaxis=:log10, yaxis=:log10,
                 xlabel="Number of steps (n)", ylabel="Error", 
                 title=f.title, label="Error", marker=:circle)
        save_plot(p, "euler-method-error-$fig_name", "6-2")
        display(p)
        readline()
        fig_name+=1
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
