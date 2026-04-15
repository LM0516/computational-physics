# Composite Trapezoidal
function composite_trapezoidal(a::Real, b::Real, f::Function, m::Int)
    h = (b - a) / m

    # Check for infinite values
    fa = isfinite(f(a)) ? f(a) : 0.0
    fb = isfinite(f(b)) ? f(b) : 0.0

    # Sum of internal nodes
    internal_sum = 0.0
    for j in 1:(m-1)
        xj = a + j * h
        internal_sum += isfinite(f(xj)) ? f(xj) : 0.0
    end

    return 0.5 * h * (fa + fb) + h * internal_sum
end

safe(x) = isinf(x) || isnan(x) ? 0.0 : x

# Composite Simpson
function composite_simpson(a::Real, b::Real, f::Function, m::Int)
    h = (b - a) / (2 * m)

    S = safe(f(a)) + safe(f(b))

    for j in 1:m
        S += 4.0 * safe(f(a + (2j - 1) * h))
    end

    for j in 1:m-1
        S += 2.0 * safe(f(a + 2j * h))
    end

    return (h / 3.0) * S
end

# Gauss-Legendre Quadrature Nodes & Weights for [-1, 1]
function glq(n::Int)
    nodes = Vector{Float64}(undef, n)
    weights = Vector{Float64}(undef, n)

    function P(x, degree)
        p0, p1 = 1.0, x
        if degree == 0
            return p0
        end
        if degree == 1
            return p1
        end

        for k in 1:(degree-1)
            p0, p1 = p1, ((2k + 1) * x * p1 - k * p0) / (k + 1)
        end
        return p1
    end

    function dP(x, degree)
        val_n = P(x, degree)
        val_n_minus_1 = P(x, degree - 1)
        return degree * (x * val_n - val_n_minus_1) / (x^2 - 1.0)
    end

    for k in 1:n
        x_guess = cos((4k - 1) * π / (4n + 2))

        # Assuming newton_method is defined externally
        root, _ = newton_method(
            x -> P(x, n),
            x -> dP(x, n),
            -1.0, 1.0,
            x_init=x_guess
        )

        nodes[k] = root
        dp_val = dP(root, n)
        weights[k] = 2.0 / ((1.0 - root^2) * dp_val^2)
    end

    return nodes, weights
end

"""
    glq_integral(f, a, b, n) -> integral_solution
Gauss-Legendre integral mapped to arbitrary interval [a, b].
"""
function glq_integral(f::Function, a::Real, b::Real, n::Int)
    # Get nodes and weights for standard domain [-1, 1]
    nodes, weights = glq(n)

    # Map nodes and weights to [a, b]
    mapped_nodes = @. 0.5 * (b - a) * nodes + 0.5 * (b + a)
    mapped_weights = @. 0.5 * (b - a) * weights

    return sum(mapped_weights .* f.(mapped_nodes))
end

# Fejér's First Rule (Using Chebyshev nodes of the first kind)
function fejer_rule(f::Function, n::Int, a::Real, b::Real)
    # Roots of the Chebyshev polynomial
    θ = [(2k - 1) * π / (2n) for k in 1:n]
    nodes = cos.(θ)

    weights = Vector{Float64}(undef, n)
    for k in 1:n
        arg_sum = 0.0
        for j in 1:(n÷2)
            # Halve the last term only if n is even and we are at n/2
            bj = (2j == n) ? 1.0 : 2.0
            arg_sum += bj / (4.0 * j^2 - 1.0) * cos(2.0 * j * θ[k])
        end
        weights[k] = 2.0 / n * (1.0 - arg_sum)
    end

    # Map to [a, b]
    mapped_nodes = @. 0.5 * (b - a) * nodes + 0.5 * (b + a)
    mapped_weights = @. 0.5 * (b - a) * weights

    return sum(mapped_weights .* f.(mapped_nodes))
end

# Clenshaw-Curtis Rule (Using Chebyshev nodes of the second kind)
function clenshaw_curtis_rule(f::Function, n::Int, a::Real=(-1.0), b::Real=(1.0))
    θ = [k * π / n for k in 0:n]
    nodes = cos.(θ)

    weights = Vector{Float64}(undef, n + 1)
    for k in 0:n
        ck_val = (k == 0 || k == n) ? 1.0 : 2.0

        arg_sum = 0.0
        for j in 1:(n÷2)
            # Fix applied: check if 2j == n to properly handle odd/even splits
            bj_val = (2j == n) ? 1.0 : 2.0
            arg_sum += bj_val / (4.0 * j^2 - 1.0) * cos(2.0 * j * θ[k+1])
        end
        weights[k+1] = ck_val / n * (1.0 - arg_sum)
    end

    # Map to [a, b]
    mapped_nodes = @. 0.5 * (b - a) * nodes + 0.5 * (b + a)
    mapped_weights = @. 0.5 * (b - a) * weights

    return sum(mapped_weights .* f.(mapped_nodes))
end

# Double Exponential Quadrature (tanh-sinh rule)
function double_exponential_quadrature(g::Function, N::Int; tol::Float64=1e-15)
    t_max = 0.0
    t_step = 0.2
    max_search_t = 10.0

    for k in 1:ceil(Int, max_search_t / t_step)
        t = k * t_step

        # Check if the function has decayed below tolerance at both tails
        if abs(g(t)) < tol && abs(g(-t)) < tol
            t_max = t
            break
        end
    end

    if t_max == 0.0
        @warn "Decay tolerance not met up to t=$(max_search_t). Using t_max = $(max_search_t)."
        t_max = max_search_t
    end

    h = t_max / N
    result = 0.0

    val_zero = g(0.0)
    if isfinite(val_zero)
        result += val_zero
    end

    for k in 1:N
        t = k * h

        val_pos = g(t)
        if isfinite(val_pos)
            result += val_pos
        end

        val_neg = g(-t)
        if isfinite(val_neg)
            result += val_neg
        end
    end

    return result * h
end
