include("../../modules/interpolations.jl")
using .Interpolations
using Plots
using LaTeXStrings

"""
    plot_interpolation(n_values, a, b, f, f_latex; interp_type="chebyshev")

Create comprehensive interpolation visualization with error analysis.

# Arguments
- `n_values::Vector{Int}`: Array of node counts to test
- `a::Real`: Left interval boundary
- `b::Real`: Right interval boundary
- `f::Function`: Function to interpolate
- `f_latex::LaTeXString`: LaTeX representation for legend
- `interp_type::String`: "chebyshev" or "lagrange"

# Returns
- Tuple of (interpolation_plot, error_plot, error_heatmap)
"""
function plot_interpolation(n_values::Vector{Int}, a::Real, b::Real, 
                            f::Function, f_latex::LaTeXString; 
                            interp_type::String="chebyshev")
    
    # Validation
    @assert interp_type in ["chebyshev", "lagrange"] "interp_type must be 'chebyshev' or 'lagrange'"
    
    # Create high-resolution evaluation points
    x_eval = range(a, b, length=40)
    f_true = f.(x_eval)
    
    # Storage for error analysis
    inf_norms = zeros(Float64, length(n_values))
    l2_norms = zeros(Float64, length(n_values))
    error_matrix = zeros(Float64, length(x_eval), length(n_values))
    
    # Configuration
    title_str = interp_type == "chebyshev" ? 
                "Chebyshev Barycentric Interpolation" : 
                "Lagrange Barycentric Interpolation"
    
    # Color scheme for different n values
    colors = palette(:viridis, length(n_values))
    
    # === Plot 1: Interpolation Comparison ===
    p1 = plot(x_eval, f_true,
        label="f(x) = $f_latex", 
        linewidth=2.5,
        color=:black,
        xlabel="x", 
        ylabel="f(x)", 
        title=title_str,
        legend=:best,
        dpi=300,
        size=(800, 500),
        margin=5Plots.mm)
    
    # Interpolate for each n and plot
    for (idx, n) in enumerate(n_values)
        # Generate nodes
        x_nodes = generate_nodes(n, a, b, interp_type)
        y_nodes = f.(x_nodes)
        
        # Compute interpolation
        p_interp = [barycentric_lagrange(x, x_nodes, f, type=interp_type) for x in x_eval]
        
        # Plot interpolation (only show a few to avoid clutter)
        if idx in [1, length(n_values)÷2, length(n_values)]
            plot!(p1, x_eval, p_interp,
                label="n = $n", 
                linewidth=1.5, 
                linestyle=:dash,
                color=colors[idx],
                alpha=0.8)
            
            # Show nodes for selected n
            scatter!(p1, x_nodes, y_nodes,
                label="",
                markersize=3,
                color=colors[idx],
                alpha=0.6,
                markerstrokewidth=0)
        else
            # Calculate without plotting
            nothing
        end
        
        # Compute errors
        errors = abs.(f_true .- p_interp)
        inf_norms[idx] = maximum(errors)
        l2_norms[idx] = sqrt(sum(errors.^2) / length(errors))
        error_matrix[:, idx] = log10.(errors .+ 1e-16)  # Add small constant for log stability
    end
    
    # === Plot 2: Error Convergence ===
    p2 = plot(xlabel="Number of nodes (n)", 
              ylabel="Error",
              title="Convergence Analysis",
              yscale=:log10,
              legend=:best,
              dpi=300,
              size=(800, 500),
              margin=5Plots.mm,
              grid=true,
              minorgrid=true)
    
    # L∞ norm
    plot!(p2, n_values, inf_norms,
          linewidth=2.5,
          marker=:circle,
          markersize=5,
          label=L"||f - p_n||_{\infty}",
          color=:red,
          markerstrokewidth=0)
    
    # L² norm
    plot!(p2, n_values, l2_norms,
          linewidth=2.5,
          marker=:square,
          markersize=5,
          label=L"||f - p_n||_{2}",
          color=:blue,
          markerstrokewidth=0)
    
    # Add reference lines for convergence rates
    if length(n_values) > 3
        # Exponential convergence reference (for Chebyshev)
        if interp_type == "chebyshev" && inf_norms[end] < inf_norms[1]
            n_ref = n_values[end÷2:end]
            exp_ref = inf_norms[end÷2] .* exp.(-0.5 .* (n_ref .- n_values[end÷2]))
            plot!(p2, n_ref, exp_ref,
                  linewidth=1.5,
                  linestyle=:dot,
                  label="Exponential ref.",
                  color=:gray,
                  alpha=0.7)
        end
    end
    
    # === Plot 3: Error Heatmap ===
    p3 = heatmap(n_values, x_eval, error_matrix,
                 xlabel="Number of nodes (n)",
                 ylabel="x",
                 title="Spatial Error Distribution (log₁₀ scale)",
                 colorbar_title=L"\log_{10}|error|",
                 color=:green,
                 dpi=300,
                 size=(800, 500),
                 margin=5Plots.mm)
    
    return p1, p2, p3
end

"""
    generate_nodes(n, a, b, interp_type)

Generate interpolation nodes based on type.
"""
function generate_nodes(n::Int, a::Real, b::Real, interp_type::String)
    if interp_type == "chebyshev"
        # Chebyshev nodes (Chebyshev-Gauss-Lobatto points)
        nodes = [cos((2*i - 1) * π / (2*n)) for i in 1:n]
        # Map from [-1,1] to [a,b]
        return @. (b - a) / 2 * nodes + (a + b) / 2
    else
        # Equally spaced nodes
        return collect(range(a, b, length=n))
    end
end

"""
    create_comparison_grid(n_values, a, b, functions; save_path=nothing)

Create a comprehensive comparison grid for multiple functions.
"""
function create_comparison_grid(n_values::Vector{Int}, a::Real, b::Real, 
                                functions::Vector{Tuple{Function, LaTeXString}};
                                save_path::Union{String, Nothing}=nothing)
    
    n_funcs = length(functions)
    plots_array = []
    
    println("Generating plots for $n_funcs functions...")
    
    for (idx, (f, f_latex)) in enumerate(functions)
        println("  Processing function $idx: $f_latex")
        
        # Chebyshev interpolation
        p_cheb, err_cheb, heat_cheb = plot_interpolation(
            n_values, a, b, f, f_latex, interp_type="chebyshev"
        )
        
        # Lagrange interpolation  
        p_lag, err_lag, heat_lag = plot_interpolation(
            n_values, a, b, f, f_latex, interp_type="lagrange"
        )
        
        push!(plots_array, p_cheb, err_cheb, p_lag, err_lag)
    end
    
    # Create grid layout
    n_rows = 2 * n_funcs
    p_combined = plot(plots_array..., 
                      layout=(n_rows, 2), 
                      size=(1600, 600 * n_funcs),
                      margin=8Plots.mm)
    
    # Save if path provided
    if !isnothing(save_path)
        savefig(p_combined, save_path)
        println("Saved plot to: $save_path")
    end
    
    return p_combined
end

function main()
    # Configuration
    a, b = -1, 1
    n_values = collect(4:4:60)
    
    # Define test functions
    functions = [
        (x -> 1 / (25 * x^2 + 1), L"\frac{1}{25x^2 + 1}"),
        (x -> tanh(5 * x + 2), L"\tanh(5x + 2)"),
        (x -> cosh(sin(x)), L"\cosh(\sin(x))"),
        (x -> sin(cosh(x)), L"\sin(\cosh(x))"),
        (x -> cosh(sin(x)), L"\cosh(\sin(x))")
    ]
    
    println("\n" * "=" ^ 60)
    println("Full Comparison: All Functions")
    println("=" ^ 60)

    p_full = create_comparison_grid(n_values, a, b, functions[1:2])
    display(p_full)

    println("\nAnalysis complete! Press Enter to exit...")
    readline()

    # === Exercise 3 ===
    #=println("\n" * "=" ^ 60)=#
    #=println("Full Comparison: All Functions")=#
    #=println("=" ^ 60)=#
    #==#
    #=p_full3 = create_comparison_grid(n_values, 0, 2*π, functions[1:5])=#
    #=display(p_full3)=#
    #==#
    #=println("\nAnalysis complete! Press Enter to exit...")=#
    #=readline()=#
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
