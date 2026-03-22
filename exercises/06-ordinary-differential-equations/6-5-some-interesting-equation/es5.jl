using ComputationalPhysics
using Plots
using LaTeXStrings

function gaussian_packet(x, x_0, k_0, σ)
    norm = (2 * π * σ^2)^(-0.25)
    space = exp(-(x - x_0)^2 / (4 * σ^2))
    phase = exp(im * k_0 * x)
    return norm * space * phase
end

function V_func(x)
    V = x > -1.0 && x < 1.0 ? -10.0 : 0.0
    return V
end

function main()
    # Space grid
    x_min, x_max = -5.0, 5.0
    dx = 0.05
    x_grid = collect(x_min:dx:x_max)

    # Energies E_n 
    H = t_independent_hamiltonian(x_grid, V_func)
    eigenvalues, eigenvectors, _ = pure_qr_algorithm(H, max_iter=600)

    # Sort and Filter Bound States (E < 0)
    bound_indices = findall(e -> e < 0, eigenvalues)
    E_bound = eigenvalues[bound_indices]
    φ_bound = eigenvectors[:, bound_indices]

    # Sort by energy level
    p = sortperm(E_bound)
    E_sorted = E_bound[p]
    φ_sorted = φ_bound[:, p]

    p = plot_generic(x_grid, V_func.(x_grid), label="Potential V(x)", color=:black, lw=2)

    for n in 1:length(E_sorted)
        # Normalize: wave_function / sqrt(dx)
        φ_n = φ_sorted[:, n] ./ sqrt(dx)

        plot!(p, x_grid, abs2.(φ_n) .+ E_sorted[n], label = L"\mathrm{State}\ %$n: |\varphi_{%$n}|^2 + E_{%$n}")
        #=plot_add!(p, x_grid, abs2.(φ_n) .+ E_sorted[n], label=latexstring("\\text{State } $n: |\\varphi_{$n}|^2 + E_{$n}"))=#
        #=plot_add!(p, x_grid, abs2.(φ_n) .+ E_sorted[n], label=L"\text{State } %$n: |\varphi_{%$n}|^2 + E_{%$n}")=#
    end

    save_plot(p, "schrodinger-equation-5", "6-5")
    display(p)
    readline()

end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    plot_init()
    main()
end
