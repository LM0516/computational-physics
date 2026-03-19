using LinearAlgebra
using BenchmarkTools
using Statistics
using Plots

include("../modules/linear_systemsV2.jl")
using .LinearSystemsV2

# --- Configuration ---
N_SIZES = [2^4, 2^5, 2^6, 2^7, 2^8] # Sizes for scaling analysis
ALGORITHMS = [
    (forward_substitution, 2, "Substitution"),
    (lu_decomposition, 3, "Decomposition"),
    (qr_mgs, 3, "Orthogonalization")
]

function automate_analysis()
    results = Dict()

    for (func, expected_k, category) in ALGORITHMS
        println("Analyzing: $func (Expected O(n^$expected_k))")
        times = Float64[]
        
        for n in N_SIZES
            # Setup specific to algorithm type
            A = rand(n, n) + n*I # Diagonally dominant for stability
            b = rand(n)
            
            if category == "Substitution"
                L = LowerTriangular(A)
                t = @belapsed $func(Matrix($L), $b)
            else
                t = @belapsed $func($A)
            end
            push!(times, t)
        end

        # Linear regression on log-log data: log(t) = k*log(n) + c
        log_n = log.(N_SIZES)
        log_t = log.(times)
        k_empirical = (sum(log_n .* log_t) - mean(log_n) * sum(log_t)) / 
                      (sum(log_n.^2) - mean(log_n) * sum(log_n))
        
        results[func] = (sizes=N_SIZES, times=times, k=k_empirical)
        println("   Empirical k: $(round(k_empirical, digits=2))\n")
    end
    return results
end

# --- Visualization ---
results = automate_analysis()
p = plot(xaxis=:log, yaxis=:log, title="Empirical Complexity Scaling", legend=:topleft)
for (func, data) in results
    plot!(data.sizes, data.times, label="$(func) (k≈$(round(data.k, digits=2)))", marker=:circle)
end
savefig("complexity_scaling.png")
