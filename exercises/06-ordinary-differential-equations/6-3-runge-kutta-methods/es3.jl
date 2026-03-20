include("../../modules/ode.jl")
using .ODE
using Plots

function main()
    f = (t, u) -> begin
        alpha = 0.1
        beta = 0.25
        y, z = u
        dy = y * (1 - alpha * y) - (y * z) / (1 + beta * y)
        dz = - z + (y * z) / (1 + beta * y)
        return [dy, dz]
    end

    a = 0.0
    b = 60.0
    n = 1000
    ics = [1, 0.01]
    p = plot(title="Phase Plane ", xlabel="y", ylabel="z", aspect_ratio=:equal)
    t, u = rk4(f, a, b, n, ics)
    z = [ui[1] for ui in u]
    y = [ui[2] for ui in u]
    plot!(p, y, z, label="IC: $ics")
 
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

