module ConsistentPlots

using Plots
using LaTeXStrings
using Measures

export initialize_style, plot_generic, plot_comparison, plot_convergence, plot_add!, save_graph

"""
Sets the global default style for Plots.jl to ensure consistency across all graphs.
"""
function initialize_style()
    default(
        fontfamily="Computer Modern",
        linewidth=2,
        framestyle=:box,
        grid=true,
        legend=:best,
        guidefontsize=12,
        tickfontsize=10,
        legendfontsize=10,
        titlefontsize=14,
        margin=5mm,
        size=(800, 500),
        dpi=300
    )
end

function plot_generic(x, y; title="", xlabel="", ylabel="", label="", kwargs...)
    plot(x, y;
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        label=label,
        kwargs...)
end

function plot_comparison(x, y_list::Vector, labels::Vector; title="", xlabel="", ylabel="", kwargs...)
    p = plot(title=title, xlabel=xlabel, ylabel=ylabel; kwargs...)

    for (i, y) in enumerate(y_list)
        lbl = length(labels) >= i ? labels[i] : "Series $i"
        plot!(p, x, y, label=lbl)
    end

    return p
end

function plot_convergence(x, errors; title="Convergence Plot", xlabel=L"h", ylabel="Error", kwargs...)
    plot(x, errors;
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        xscale=:log10,
        yscale=:log10,
        marker=:circle,
        kwargs...)
end

"""
Add a series to an existing plot `p` while inheriting the global style.
"""
function plot_add!(p, x, y; label="", kwargs...)
    plot!(p, x, y; label=label, kwargs...)
end

"""
Saves the plot on the right folder for the repor. 
Good for updating the figures in the report.
"""
function save_graph(p, filename, chapter_dir)
    path = joinpath("report", "figures", string(chapter_dir), "$(filename).pdf")
    mkpath(dirname(path))
    savefig(p, path)
    println("Saved figure to: $path")
end

end # module
