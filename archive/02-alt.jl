include("../../modules/ode.jl")
using .ODE
using Plots

function main()
    f_a = (t, u) -> begin
        x, y = u
        dx = -4y + x * (1 - x^2 - y^2)
        dy = 4x + y * (1 - x^2 - y^2)
        return [dx, dy]
    end

    f_b = (t, u) -> begin
        x, y = u
        dx = -4y - 1 / 4 * x * (1 - x^2 - y^2) * (4 - x^2 - y^2)
        dy = 4x - 1 / 4 * y * (1 - x^2 - y^2) * (4 - x^2 - y^2)
        return [dx, dy]
    end

    a = 0.0
    b = 10.0
    n = 1000

    ics_a = [[0.1, 0.0], [0.0, 1.9]]
    p_a = plot(title="Phase Plane (a)", xlabel="x", ylabel="y", aspect_ratio=:equal)
    for u0 in ics_a
        t, u = rk4(f_a, a, b, n, u0)
        x = [ui[1] for ui in u]
        y = [ui[2] for ui in u]
        plot!(p_a, x, y, label="IC: $u0")
    end

    ics_b = [[0.95, 0.0], [0.0, 1.05], [-2.5, 0.0]]
    p_b = plot(title="Phase Plane (b)", xlabel="x", ylabel="y", aspect_ratio=:equal)
    for u0 in ics_b
        t, u = rk4(f_b, a, b, n, u0)
        x = [ui[1] for ui in u]
        y = [ui[2] for ui in u]
        plot!(p_b, x, y, label="IC: $u0")
    end

    display(plot(p_a, p_b, layout=(1, 2), size=(800, 400)))
    readline()
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
