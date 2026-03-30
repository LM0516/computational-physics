using ComputationalPhysics

f = x -> (exp(x + 1) - 2 - x)
df = x -> (exp(x + 1) - 1)
a = -2
b = 2

r, x = newton_method(f, df, a, b)

q_values, C_values = conv_rate_and_asymp_const(x, r)

@show q_values[end], C_values[end]
