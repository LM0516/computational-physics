using ComputationalPhysics
using Plots

function main()
    f = (t, u) -> begin
        alpha = 0.1
        beta = 0.25
        y, z = u
        dy = y * (1 - alpha * y) - (y * z) / (1 + beta * y)
        dz = -z + (y * z) / (1 + beta * y)
        return [dy, dz]
    end

    a = 0.0
    b = 60.0
    n = 1000
    ics = [1, 0.01]

    t, u = rk4(f, a, b, n, ics)
    y = [ui[1] for ui in u]
    z = [ui[2] for ui in u]

    p_time = plot_generic(t, y, label="Prey (y)", xlabel="Time (t)", ylabel="Population", linewidth=1.5)
    plot_add!(p_time, t, z, label="Predator (z)", linewidth=1.5)
    title!(p_time, "Population Dynamics")

    p_phase = plot(title="Phase Plane ", xlabel="y", ylabel="z", aspect_ratio=:equal)
    plot_add!(p_phase, y, z, label="IC: $ics")

    p = multi_plot(p_time, p_phase)
    save_plot(p, "predator-prey-model", "6-3")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end

