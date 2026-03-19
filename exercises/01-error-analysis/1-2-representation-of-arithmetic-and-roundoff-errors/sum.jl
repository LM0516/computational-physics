using Plots
using LaTeXStrings

function somma(a, N::Int)
  s1 = zeros(eltype(a), N)
  s2 = zeros(eltype(a), N)
  for n = 1:N
    for k = 1:n
      s1[n] += a[k]
      s2[n] += a[N-k+1]
    end
  end
  err1 = @. abs(s1 - pi^2 / 6)
  err2 = @. abs(s2 - pi^2 / 6)
  return s1, s2, err1, err2
end

N = 1000

a32 = Float32[1 / n^2 for n = 1:N]
a64 = Float64[1 / n^2 for n = 1:N]

# Single precision
no32, ro32, err32_so, err32_ro = somma(a32, N)

# Double precision
no64, ro64, err64_so, err64_ro = somma(a64, N)

# === Plot ===
p = scatter(err32_ro, label="Single normal", color=1)
scatter!(err32_so, label="Single reverse", color=2)
# plot!(no64, label="Double normal", color=2)
# plot!(ro64, label="Double reverse", color=2, ls=:dash)
display(p)
readline()
# p1 = scatter(title="Convergence errors", legend=:bottomleft, xaxis=(L"n", [1, N]), yaxis=L"\pi^2/6 S_N")

# id=1:100:N
# scatter!(id, err32_so[id], label="1 \to N", color=1)
# scatter!(id, err32_ro[id], label="1 \to N", color=1, ls=:dash)
# scatter!(id, err64_so[id], label="1 \to N", color=2)
# scatter!(id, err64_ro[id], label="1 \to N", color=2, ls=:dash)
# plot!(1:N, n->1.0/n, ls=:dash, label=L"0(1/n)", color=1)
# xlims!(1500, 7000)
# ylims!(0, 5e-4)
# vline!([4096], color=1, label=L"n = 4096")

# display(p1)
# readline()

# # === Ratio ===
# ratio = err32_so ./ err32_ro  
# minimum_ratio = minimum(ratio)
# n_min = findfirst(ratio .== minimum_ratio)
# println("Il minimo del rapporto è ", minimum_ratio, " e si trova a n = ", n_min)
# p2 = plot(ratio, xlabel="n", ylabel="Rapporto degli errori", 
#      title="Rapporto tra errore normal e reverse", 
#      yscale=:log10, label="Rapporto")

# display(p2)
# readline()
