module ODE
export euler_method, ie2, rk4_old, rk4

function euler_method(f::Function, a::Number, b::Number, n::Int, u0)
    # Determine the type of elements in the solution array
    # If u0 is a Number, we want a Vector{Float64} (or whatever type u0 promotes to with arithmetic)
    # If u0 is a Vector, we want a Vector{Vector{Float64}}

    # We can initialize the solution array with the type of u0, converted to float
    # This handles cases where u0 is Int but the solution will be Float64
    u = Vector{typeof(float(u0))}(undef, n + 1)
    u[1] = u0
    h = (b - a) / n
    t = [a + (i - 1) * h for i in 1:n+1]
    for i in 1:n
        u[i+1] = u[i] + h * f(t[i], u[i])
    end
    return t, u
end

function ie2(f::Function, a::Number, b::Number, n::Int, u0)
    u = Vector{typeof(float(u0))}(undef, n + 1)
    u[1] = u0
    h = (b - a) / n
    t = [a + (i - 1) * h for i in 1:n+1]

    for i in 1:n
        k1 = h * f(t[i], u[i])
        u[i+1] = u[i] + h * f(t[i] + 0.5 * h, u[i] + 0.5 * k1)
    end

    return t, u
end

function rk4(f::Function, a::Number, b::Number, n::Int, u0)
    u = Vector{typeof(float(u0))}(undef, n + 1)
    u[1] = u0
    h = (b - a) / n
    t = [a + (i - 1) * h for i in 1:n+1]

    for i in 1:n
        k1 = h * f(t[i], u[i])
        k2 = h * f(t[i] + 0.5 * h, u[i] + 0.5 * k1)
        k3 = h * f(t[i] + 0.5 * h, u[i] + 0.5 * k2)
        k4 = h * f(t[i] + h, u[i] + k3)
        u[i+1] = u[i] + (1 / 6) * (k1 + 2k2 + 2k3 + k4)
    end

    return t, u
end

end
