using ComputationalPhysics

function gaussian_packet(x, x_0, k_0, σ)
    norm = (2 * π * σ^2)^(-0.25)
    space = exp(-(x - x_0)^2 / (4 * σ^2))
    phase = exp(im * k_0 * x)
    return norm * space * phase
end

function V_func(x)
    V = x > -1.0 && x < 1.0 ? 2.0 : 0.0
    return V
end

function main()
    # Space grid
    x_min, x_max = -100.0, 100.0
    dx = 0.1
    x_grid = collect(x_min:dx:x_max)

    # Time grid
    t_min, t_max = 0.0, 40.0
    dt = 0.01
    t_grid = collect(t_min:dt:t_max)

    # Initial packet parameters
    x_0 = -20.0
    k_0 = 2.0
    σ = 2.0

    # Create initial wavefunction
    psi_init = gaussian_packet.(x_grid, x_0, k_0, σ)
    snap = schrodinger_solver1D(x_grid, t_grid, V_func, psi_init)

    # === Plotting ===
    target_times = [0.0, 10.0, 30.0]
    p = plot_snapshots(target_times, snap, x_grid, t_min, dt, dx)

    save_plot(p, "schrodinger-equaiton-3", "6-5")
    display(p)
    readline()

    # === R and T coefficients ===
    target_time = 30.0
    R, T = compute_RT(target_time, t_min, dt, snap, x_grid, -1.0, 1.0, V_func; graph=true)

    println("Reflection Coefficient (R):  ", round(R, digits=4))
    println("Transmission Coefficient (T):", round(T, digits=4))
    println("Total Probability (R + T):   ", round(R + T, digits=4))

end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
