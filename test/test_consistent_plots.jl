
push!(LOAD_PATH, joinpath(@__DIR__, "../modules"))

using ConsistentPlots
using Plots
using LaTeXStrings

# Initialize style
initialize_style()

# 1. Generic Plot
x = range(0, 10, length=100)
y = sin.(x)
p1 = plot_generic(x, y,
    title="Generic Plot Test (Sin)",
    xlabel=L"x",
    ylabel=L"\sin(x)",
    label="Sine Wave",
    lw=3,
)
save_graph("test_generic", "test_consistent_plots")

# 2. Comparison Plot
y1 = sin.(x)
y2 = cos.(x)
y3 = sin.(x) .* exp.(-x / 5)
p2 = plot_comparison(x, [y1, y2, y3],
    [L"\sin(x)", L"\cos(x)", L"\sin(x)e^{-x/5}"],
    title="Comparison Test",
    xlabel=L"x",
    ylabel="Amplitude"
)
save_graph("test_comparison", "test_consistent_plots")

# 3. Convergence Plot (Log-Log)
h = [10.0^(-k) for k in 1:10]
errors = h .^ 2 # 2nd order convergence
p3 = plot_convergence(h, errors,
    title="Convergence Test (O(h^2))",
    label=L"Error \approx h^2"
)
save_graph("test_convergence", "test_consistent_plots")

println("All test plots generated and saved to Relazione/figures/test_consistent_plots/")
