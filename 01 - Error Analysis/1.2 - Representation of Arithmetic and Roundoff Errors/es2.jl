# ef642.jl
#
using Plots
using Printf
using LaTeXStrings

function iteration_sequence(x0, x1, n)
  xs = zeros(typeof(x0), n+1)

  xs[1] = x0
  xs[2] = x1

  for k in 2:n
    xs[k+1] = 111 - (1130 - 3000 / xs[k-1]) / xs[k]
  end

  return xs
end 

N = 34

# iniziali in forma razionale (per esattezza)
x0_rat = BigInt(11) // 2        # 11/2
x1_rat = BigInt(61) // 11       # 61/11

# versione Float32
x0_f32 = Float32(11)/2
x1_f32 = Float32(61)/11

# versione Float64
x0_f64 = Float64(11)/2
x1_f64 = Float64(61)/11

# versione BigFloat (alta precisione)
setprecision(BigFloat, 128)     # imposta precisione per BigFloat (128 bit significativi)
x0_b = BigFloat(11)/BigFloat(2)
x1_b = BigFloat(61)/BigFloat(11)

# itera
xs_f32 = iteration_sequence(x0_f32, x1_f32, N)
xs_f64= iteration_sequence(x0_f64, x1_f64, N)
xs_b = iteration_sequence(x0_b, x1_b, N)
xs_r = iteration_sequence(x0_rat, x1_rat, N)

# stampa tabella comparativa per i primi 40 indici
println(" k   Float64                 BigFloat (128)             Rational (exact as float)")
println("--------------------------------------------------------------------------")
for k in 1:N+1
    xf32 = xs_f32[k]
    xf64 = xs_f64[k]
    xb = xs_b[k]
    xr = float(xs_r[k])   # convertiamo la razionale in Float64 solo per stampa compatta
    @printf("%2d  %22.12g  %22.12g  %25.14g  %22.12g\n", k-1, xf32, xf64, xb, xr)
end

# plotta i valori
plot(0:N, xs_f32, label="Float32", linewidth=1)
plot!(0:N, xs_f64, label="Float64", linewidth=1)
plot!(0:N, xs_b, label="BigFloat", linewidth=1)
plot!(0:N, float.(xs_r), label="Rational (exact)", linewidth=1)
hline!([6], label="convergence", linestyle=:dash, color=:black)
plot!(legend=:bottomleft)
xlabel!(L"k")
ylabel!(L"x_k")
title!("Comparison of numerical precision")

display(plot!())
readline()
