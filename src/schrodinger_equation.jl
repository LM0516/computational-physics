const ħ = 1.0
const m = 1.0

function t_independent_hamiltonian(x_grid::Vector{Float64}, potential_func::Function)
    N = length(x_grid)
    dx = x_grid[2] - x_grid[1]

    # Kinetic energy coefficients
    h_coeff = 5 * ħ^2 / (4m * dx^2)
    b = -ħ^2 / (24m * dx^2)

    # Potential vector
    V_vec = potential_func.(x_grid)

    # Construct Sparse Matrix H
    d0 = h_coeff .+ V_vec
    d1 = fill(16 * b, N - 1)
    d2 = fill(-b, N - 2)

    H = spdiagm(
        0 => d0,
        1 => d1,
        -1 => d1,
        2 => d2,
        -2 => d2
    )

    return H
end

function schrodinger_solver1D(x_grid::Vector{Float64}, t_grid::Vector{Float64}, potential_func::Function, psi_init::Vector{ComplexF64})
    met = first(methods(potential_func))
    if met.nargs - 1 == 1
        H = t_independent_hamiltonian(x_grid, potential_func)
        M = (-im / ħ) * H
        dt = t_grid[2] - t_grid[1]

        snap = Vector{Vector{ComplexF64}}(undef, length(t_grid) + 1)
        snap[1] = psi_init

        for i in 1:length(t_grid)
            k1 = M * snap[i]
            k2 = M * (snap[i] .+ (0.5 * dt) .* k1)
            k3 = M * (snap[i] .+ (0.5 * dt) .* k2)
            k4 = M * (snap[i] .+ dt .* k3)

            snap[i+1] = snap[i] .+ (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        end

        return snap
    elseif met.nargs - 1 == 2
        @warn "The potential function is time dependent, using the time dependent solver"
        exit()
    else
        @error "The potential function has an invalid number of arguments ($met.nargs)"
        exit()
    end
end

function compute_RT(t_target::Float64, t_min::Float64, dt::Float64, snap::Vector{Vector{ComplexF64}}, x_grid::Vector{Float64}, barrier_left::Float64, barrier_right::Float64, V_barrier::Function; graph::Bool=false)
    idx_time = round(Int, (t_target - t_min) / dt)
    psi = snap[idx_time]
    dx = x_grid[2] - x_grid[1]
    prob_density = abs2.(psi)

    idx_R = x_grid .< barrier_left
    idx_T = x_grid .> barrier_right

    R = sum(prob_density[idx_R]) * dx
    T = sum(prob_density[idx_T]) * dx

    if graph
        p_barrier = 1.0 - (R + T)
        println("Probability inside barrier:  ", round(p_barrier, digits=4))
        # E. Optional: Plotting
        p = plot(x_grid, abs2.(psi),
            label="t=30",
            title="Scattering at t=30",
            xlabel="x", ylabel=L"|\psi(x,t)|^{2}",
            xlims=(-60, 60), # Zoom in to see the packets
            fill=(0, 0.2, :blue)
        )

        # Draw the barrier area for visualization
        vplot = [V_barrier(x) for x in x_grid]
        plot!(p, x_grid, vplot .* 0.1, label="Scaled V(x)", color=:red, alpha=0.5)

        display(p)
        readline()
    end

    return R, T
end

function plot_snapshots(target_times::Vector{Float64}, snap::Vector{Vector{ComplexF64}}, x_grid::Vector{Float64}, t_min::Float64, dt::Float64, dx::Float64)
    p = plot(xlabel="Position (x)",
        ylabel=L"|\psi(x,t)|^2",
        legend=:topright)

    for t_target in target_times
        idx = round(Int, (t_target - t_min) / dt) + 1

        if 1 <= idx <= length(snap)
            # Extract wavefunction at that index
            psi_t = snap[idx]

            # Calculate Probability Density: |ψ|^2
            prob_density = abs2.(psi_t)

            # Add to the plot (plot! adds to the existing figure)
            plot!(p, x_grid, prob_density, label="t = $t_target", lw=2)

            # Exercise 1(b): Verify Norm
            current_norm = sum(prob_density) * dx
            println("  t=$t_target: Norm ≈ $current_norm")
        else
            println("  Warning: Time $t_target is outside the simulation range.")
        end
    end

    display(p)
    readline()
end
