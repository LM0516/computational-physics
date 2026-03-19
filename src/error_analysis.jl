function mclaurin_series(x::Float64, n::Int)
    sum = 0.0
    for k in 0:n
        sum += x^k / factorial(k + 1)
    end
    return sum
end

# === Double-pass variance (two loops) ===
function var_double_pass32(x::Vector{Float32})
    n = length(x)
    mean_x = sum(x) / n
    return sum((xi - mean_x)^2 for xi in x) / (n - 1)
end

function var_double_pass64(x::Vector{Float64})
    n = length(x)
    mean_x = sum(x) / n
    return sum((xi - mean_x)^2 for xi in x) / (n - 1)
end

# === Single-pass variance ===
function var_single_pass32(x::Vector{Float32})
    n = length(x)
    u = sum(xi^2 for xi in x)
    v = sum(x)
    return (u - v^2 / n) / (n - 1)
end

function var_single_pass64(x::Vector{Float64})
    n = length(x)
    u = sum(xi^2 for xi in x)
    v = sum(x)
    return (u - v^2 / n) / (n - 1)
end

function kappa_f(x, func::Function)
    return abs(x * ForwardDiff.derivative(func, x) / func(x))
end

