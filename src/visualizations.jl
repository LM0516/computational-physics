"""
Initializes the global default style for Plots.jl to ensure consistency 
across all graphs in the project.
"""
function plot_init()
    default(
        fontfamily="Computer Modern", # Great for LaTeX report consistency
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
    println("Plotting style initialized.")
end

"""
Creates a standard 2D line plot with optional titles and labels.
"""
function plot_generic(x, y; title="", xlabel="", ylabel="", label="", kwargs...)
    return plot(x, y;
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        label=label,
        kwargs...)
end

"""
Adds a line series to an existing plot `p`.
"""
function plot_add!(p, x, y; label="", kwargs...)
    return plot!(p, x, y; label=label, kwargs...)
end

"""
Adds a scatter series (points/markers) to an existing plot `p`.
"""
function scatter_add!(p, x, y; label="", kwargs...)
    return scatter!(p, x, y; label=label, kwargs...)
end

"""
    multi_plot(plots...; layout, kwargs...)

Combines multiple individual plot objects into a single figure (subplots).
Example: `multi_plot(p1, p2, p3; layout=(3, 1))`
"""
function multi_plot(plots...; layout=length(plots), kwargs...)
    return plot(plots...; layout=layout, kwargs...)
end

"""
Saves the plot to the specified chapter directory in the report folder.
Automatically creates the directory if it doesn't exist.
"""
function save_plot(p, filename, chapter_dir; folder="report")
    path = joinpath(folder, "figures", string(chapter_dir), "$(filename).pdf")
    mkpath(dirname(path)) # Safely creates nested directories if missing
    savefig(p, path)
    println("Saved figure to: $path")
end

"""
Saves an Animation object as a GIF in the output folder.
"""
function save_gif(anim::Animation, filename, chapter_dir; fps=15)
    path = joinpath("output", "gifs", string(chapter_dir), "$(filename).gif")
    mkpath(dirname(path))
    gif(anim, path, fps=fps) # Note: use `gif()` instead of `savefig()`
    println("Saved GIF to: $path")
end
