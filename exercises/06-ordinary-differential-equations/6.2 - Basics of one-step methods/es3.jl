include("../../modules/ode.jl")
using .ODE
using Plots
using LaTeXStrings

function main()
    functions = [
        (
            f=(t, u) -> [u[2], sin(2t) - 9u[1]],
            a=0.0,
            b=2π,
            u0=[2.0, 1.0],
            g=t -> [
                (1 / 5) * sin(3t) + 2 * cos(3t) + (1 / 5) * sin(2t), # y(t)
                (3 / 5) * cos(3t) - 6 * sin(3t) + (2 / 5) * cos(2t)  # y'(t)
            ],
            title=L"y'' + 9y = \sin(2t)"
        ),
        (
            f=(t, u) -> [u[2], 4t + 4u[1]],
            a=0.0,
            b=1.5,
            u0=[2.0, -1.0],
            g=t -> [
                exp(2t) + exp(-2t) - t,       # y(t)
                2 * exp(2t) - 2 * exp(-2t) - 1    # y'(t)
            ],
            title=L"y'' - 4y = 4t"
        ),
        (
            f=(t, u) -> [u[2], t - 4u[2] - 4u[1]],
            a=0.0,
            b=4.0,
            u0=[1.0, 0.75],
            g=t -> [
                (3t + 5 / 4) * exp(-2t) + (t - 1) / 4,      # y(t)
                (3 - 2 * (3t + 5 / 4)) * exp(-2t) + 1 / 4     # y'(t)
            ],
            title=L"y'' + 4y' + 4y = t"
        )
    ]

    n = 1000

    for (i, func) in enumerate(functions)
        t, u = euler_method(func.f, func.a, func.b, n, func.u0)

        # Extract components
        y_num = [val[1] for val in u]
        yp_num = [val[2] for val in u]

        # Analytical solution
        exact = func.g.(t)
        y_exact = [val[1] for val in exact]
        yp_exact = [val[2] for val in exact]

        # Plot solution and derivative
        p1 = plot(t, y_num, label="Numerical y(t)", xlabel="Time (t)", ylabel="Value", title=func.title, lw=2)
        plot!(p1, t, yp_num, label="Numerical y'(t)", lw=2, linestyle=:dash)
        plot!(p1, t, y_exact, label="Exact y(t)", color=:black, linestyle=:dot)
        plot!(p1, t, yp_exact, label="Exact y'(t)", color=:grey, linestyle=:dot)

        display(p1)
        println("Press Enter to see error plot for problem $(i)...")
        readline()

        # Plot errors
        error_y = abs.(y_num .- y_exact)
        error_yp = abs.(yp_num .- yp_exact)

        p2 = plot(t, error_y, label="Error y(t)", xlabel="Time (t)", ylabel="Absolute Error", title="Error for $(func.title)", lw=2)
        plot!(p2, t, error_yp, label="Error y'(t)", lw=2, linestyle=:dash)

        display(p2)
        println("Press Enter to continue...")
        readline()
    end
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
