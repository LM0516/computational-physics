using ComputationalPhysics

function gaussian_packet(x, x_0, k_0, σ)
    norm = (2 * π * σ^2)^(-0.25)
    space = exp(-(x - x_0)^2 / (4 * σ^2))
    phase = exp(im * k_0 * x)
    return norm * space * phase
end

function V_func(x)
    return 0.0
end

function main()
    # Space grid
    x_min, x_max = -10.0, 10.0
    dx = 0.1
    x_grid = collect(x_min:dx:x_max)

    # Time grid
    t_min, t_max = 0.0, 20.0
    dt = 0.01
    t_grid = collect(t_min:dt:t_max)

    # Initial packet parameters
    x_0 = 0.0
    k_0 = 0.0
    σ = 1.0

    # Create initial wavefunction
    psi_init = gaussian_packet.(x_grid, x_0, k_0, σ)
    snap = schrodinger_solver1D(x_grid, t_grid, V_func, psi_init)

    # === Plotting ===
    target_times = [0.0, 2.0, 3.0, 4.0, 5.0]
    plot_snapshots(target_times, snap, x_grid, t_min, dt, dx)

end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    main()
end
