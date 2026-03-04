module NumericalIntegration
include("../modules/nonlinear_equatioins.jl")
using .NonlinearEquations
export composite_trapezoidal, composite_simpson, gauss_legendre_quadrature, fajer_rule, clenshaw_curtis_rule, double_exponential_quadrature

# FIXME: ricontrollare la funzione, non comprende tutti i casi possibili. Serve solo per risolvere l'esercizio
# in teoria dovrebbe essere in grado di trattare tutti i casi in cui la funzione non e' nel dominio e prendere
# il limite della funzione stessa.
# Questo problema potrei averlo solo per asintoti verticali visto che per applicare questa regola l'integrale
# deve essere finito.
function composite_trapezoidal(a::Real, b::Real, f::Function, m::Int)
    h = (b - a) / m

    # Check for infinite values
    fa = isfinite(f(a)) ? f(a) : 0
    fb = isfinite(f(b)) ? f(b) : 0

    # Sum of internal nodes: sum(f(a + j*h))
    internal_sum = 0.0
    for j in 1:(m-1)
        xj = a + j * h
        internal_sum += isfinite(f(xj)) ? f(xj) : 0
    end

    # Formula: h/2 * [f(a) + f(b) + 2 * internal_sum]
    # This simplifies to: h/2 * (f(a) + f(b)) + h * internal_sum
    T = 0.5h * (fa + fb) + h * internal_sum
    return T
end

function composite_simpson(a::Real, b::Real, f::Function, m::Int)
    h = (b - a) / (2m)

    # Check for infinite values
    fa = isinf(f(a)) ? 0 : f(a)
    fb = isinf(f(b)) ? 0 : f(b)

    int_sum_1 = 0.0
    int_sum_2 = 0.0

    for j in 1:(m-1)
        int_sum_1 += isinf(f(a + (2j - 1) * h)) ? 0 : f(a + (2j - 1) * h)
    end

    for j in 1:m
        int_sum_2 += isinf(f(a + 2j * h)) ? 0 : f(a + 2j * h)
    end

    S = (h / 3) * (fa + fb + 4int_sum_1 + 2int_sum_2)

    return S
end

function gauss_legendre_quadrature(n::Int; a::Real=-1.0, b::Real=1.0)
    nodes = zeros(n)
    weights = zeros(n)

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
        return degree * (x * val_n - val_n_minus_1) / (x^2 - 1)
    end

    for k in 1:n
        x_guess = cos((4k - 1) * π / (4n + 2))

        root, _ = newton_method(
            x -> P(x, n),
            x -> dP(x, n),
            a, b,
            x_init=x_guess
        )

        nodes[k] = root

        dp_val = dP(root, n)
        weights[k] = 2 / ((1 - root^2) * dp_val^2)
    end

    return nodes, weights

end

function fajer_rule(f::Function, n::Int, a::Real, b::Real)
    nodes, weights = gauss_legendre_quadrature(n, a=a, b=b)
    integral = sum(@. weights * f(nodes))
    return integral
end

function clenshaw_curtis_rule(f::Function, n::Int)
    θ = [k * π / n for k in 0:n]
    nodes = cos.(θ)

    # Calculate weights
    weights = zeros(n + 1)
    for k in 0:n
        ck_val = (k == 0 || k == n) ? 1.0 : 2.0

        arg_sum = 0.0
        for j in 1:(n÷2)
            bj_val = (j == n ÷ 2) ? 1.0 : 2.0
            arg_sum += bj_val / (4.0 * j^2 - 1.0) * cos(2.0 * j * θ[k+1])
        end
        weights[k+1] = ck_val / n * (1.0 - arg_sum)
    end

    return sum(weights .* f.(nodes))
end

function double_exponential_quadrature(g::Function, N::Int; tol::Float64=1e-15)
    t_max = 0.0
    t_step = 0.2
    max_search_t = 10.0

    for k in 1:ceil(Int, max_search_t / t_step)
        t = k * t_step

        val_pos = abs(g(t))
        val_neg = abs(g(-t))

        # Check if the function has decayed below tolerance at both tails
        if val_pos < tol && val_neg < tol
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

end
