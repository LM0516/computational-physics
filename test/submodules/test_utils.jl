using Printf

"""
    empirical_order(ns, errors) -> order

Estimate the empirical convergence order from a sequence of discretization sizes
`ns` and the corresponding absolute errors `errors`.

The estimate is obtained by fitting the model

`error(n) ≈ C * n^{-p}`

in log-log coordinates, so the returned value is the fitted slope `p`.
This is used in the tests to verify expected first-, second-, and fourth-order
behavior without hard-coding a single mesh refinement ratio.
"""
function empirical_order(ns, errors)
    xs = log.(Float64.(ns))
    ys = log.(Float64.(errors))
    slope = sum((xs .- mean(xs)) .* (ys .- mean(ys))) / sum((xs .- mean(xs)) .^ 2)
    return -slope
end

"""
    counted_scalar_function(f) -> wrapped, calls

Wrap a function `f` with a lightweight call counter.

Returns:
- `wrapped`: a function with the same arguments and return value as `f`
- `calls`: a `Ref{Int}` whose value is incremented every time `wrapped` is called

This helper lets the tests measure computational cost in a stable way using
function evaluations instead of wall-clock timing.
"""
function counted_scalar_function(f)
    calls = Ref(0)
    wrapped = function (args...)
        calls[] += 1
        return f(args...)
    end
    return wrapped, calls
end

"""
    sample_tail(values; k=3) -> tail_values

Return the last `k` entries of `values`, or the full collection if it has fewer
than `k` elements.

The tests use this helper when only the asymptotic part of a sequence is
relevant, for example when reading off the tail of observed convergence rates.
"""
function sample_tail(values; k=3)
    n = length(values)
    first_idx = max(1, n - k + 1)
    return values[first_idx:end]
end

"""
    format_metric(value; digits=6) -> formatted

Format a numeric metric for readable test output. Floating-point values are
rounded to a fixed number of digits, while vectors are formatted elementwise.
"""
function format_metric(value; digits=6)
    if value isa AbstractVector
        return "[" * join(format_metric.(value; digits=digits), ", ") * "]"
    elseif ismissing(value)
        return "n/a"
    elseif value isa AbstractFloat
        return @sprintf("%.*f", digits, value)
    else
        return string(value)
    end
end

"""
    print_metric_table(title, rows; headers=("Metric", "Value"))

Print a compact ASCII table for test metrics.

Each row must be a two-element tuple `(name, value)`. This is intended for
human-readable summaries of convergence rates, iteration counts, and function
evaluation costs during the test run.
"""
function print_metric_table(title, rows; headers=("Metric", "Value"))
    formatted_rows = [(string(name), format_metric(value)) for (name, value) in rows]
    width1 = maximum(length.([headers[1]; first.(formatted_rows)]))
    width2 = maximum(length.([headers[2]; last.(formatted_rows)]))
    sep = "+" * repeat("-", width1 + 2) * "+" * repeat("-", width2 + 2) * "+"

    println()
    println(title)
    println(sep)
    println("| ", rpad(headers[1], width1), " | ", rpad(headers[2], width2), " |")
    println(sep)
    for (name, value) in formatted_rows
        println("| ", rpad(name, width1), " | ", rpad(value, width2), " |")
    end
    println(sep)
end

"""
    print_method_table(title, rows; headers=("Function", "Convergence", "Cost"))

Print a compact ASCII table where each row summarizes one numerical method with
separate columns for convergence information and computational cost.

Each row must be a three-element tuple `(method, convergence, cost)`.
"""
function print_method_table(title, rows; headers=("Function", "Convergence", "Cost"))
    formatted_rows = [(string(name), format_metric(conv), format_metric(cost)) for (name, conv, cost) in rows]
    width1 = maximum(length.([headers[1]; getindex.(formatted_rows, 1)]))
    width2 = maximum(length.([headers[2]; getindex.(formatted_rows, 2)]))
    width3 = maximum(length.([headers[3]; getindex.(formatted_rows, 3)]))
    sep = "+" * repeat("-", width1 + 2) * "+" * repeat("-", width2 + 2) * "+" * repeat("-", width3 + 2) * "+"

    println()
    println(title)
    println(sep)
    println("| ", rpad(headers[1], width1), " | ", rpad(headers[2], width2), " | ", rpad(headers[3], width3), " |")
    println(sep)
    for (name, conv, cost) in formatted_rows
        println("| ", rpad(name, width1), " | ", rpad(conv, width2), " | ", rpad(cost, width3), " |")
    end
    println(sep)
end
