"""
Prepares an array `y` for logarithmic plotting.
Methods:
- `:epsilon` : Takes max(abs(y), epsilon) to prevent -Inf. (Default)
- `:missing` : Replaces y <= 0 with `missing` to break the plot line cleanly.
"""
function make_log_safe(y::AbstractArray; epsilon=1e-16, method=:epsilon)
    if method == :epsilon
        return max.(abs.(y), epsilon)
    elseif method == :missing
        return [v > 0 ? v : missing for v in y]
    else
        throw(ArgumentError("Unknown method. Use :epsilon or :missing"))
    end
end
