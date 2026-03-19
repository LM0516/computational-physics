using Plots

@time begin
    N = 10000
    a = Float32[1 / n^2 for n = 1:N]

    s1 = zeros(Float32, N)
    s2 = zeros(Float32, N)

    for n = 1
        for k = 1:n
            s1[n] += a[k]
            s2[n] += a[N-k+1]
        end
    end
end

err1 = @. abs(s1 - pi^2 / 6)
err2 = @. abs(s2 - pi^2 / 6)

scatter(1:N, [err1 err2], yscale=:log10, xlabel="n", ylabel="Errore assoluto", label=["Normal precision" "Reverse precision"], title="Errore assoluto in funzione di n")
