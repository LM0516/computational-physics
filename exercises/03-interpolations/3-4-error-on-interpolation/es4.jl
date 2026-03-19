using ComputationalPhysics
using Plots
using LaTeXStrings

function phi_forward(x)
    """Map from x in (-1,1) to z on real line"""
    return 2*x / (1 - x^2)
end

function phi_inverse(z)
    """Map from z on real line back to x in (-1,1)"""
    # From solving z = 2x/(1-x^2):
    # z(1-x^2) = 2x
    # z - zx^2 = 2x
    # zx^2 + 2x - z = 0
    # Using quadratic formula: x = (-2 ± sqrt(4 + 4z^2))/(2z)
    # For the branch in (-1,1): x = (-1 + sqrt(1 + z^2))/z
    
    if z ≈ 0
        return 0.0
    else
        return (-1 + sqrt(1 + z^2)) / z
    end
end

function main()
    n = 30
    
    # Generate Chebyshev nodes on (-1, 1)
    # Using first kind: x_i = cos((2i-1)π/(2n))
    x_nodes_cheb = [cos((2*i - 1) * π / (2*n)) for i in 1:n]
    
    # Transform to real line
    z_nodes = phi_forward.(x_nodes_cheb)
    
    # Function to interpolate
    f(z) = (z^2 - 2*z + 2)^(-1)
    
    # Get function values at transformed nodes
    y_nodes = f.(z_nodes)
    
    # Evaluation points on real line
    z_eval = range(-6.0, 6.0, length=1000)
    
    # Transform to x-space for interpolation
    x_eval = phi_inverse.(z_eval)
    
    # Evaluate interpolant
    p_interp = [barycentric_lagrange(x, x_nodes_cheb, (xi -> f(phi_forward(xi))), type="chebyshev") 
                for x in x_eval]
    
    # Plot
    p = plot_generic(z_eval, f.(z_eval),
        label=L"f(z) = (z^2 - 2z + 2)^{-1}", linewidth=1,
        xlabel="z", ylabel="f(z)",
        title="Real Line Interpolation via Chebyshev Nodes")
    
    plot_add!(p, z_eval, p_interp,
        label="Interpolant q(z)", linewidth=2, linestyle=:dash)
    
    # Mark nodes (only those in [-6,6])
    z_nodes_in_range = z_nodes[findall(z -> -6 ≤ z ≤ 6, z_nodes)]
    scatter_add!(p, z_nodes_in_range, f.(z_nodes_in_range),
        label="Nodes in [-6,6]")
    
    # Compute error
    error = abs.(f.(z_eval) - p_interp)
    max_error = maximum(error)
    println("Maximum interpolation error: $max_error")
    
    save_plot(p, "cheb_interpolation_plot", "3-4")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init() 
    main()
end
