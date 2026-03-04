include("../../modules/linear_systemsV2.jl")
using .LinearSystemsV2
using LinearAlgebra

println("Modified Gram-Schmidt QR Decomposition Implementation")
println()

# === Matrix 1 ===
A1 = [2 -1
      0 1
      -2 2]

println("Matrix 1:")
display(A1)
println()
Q1, R1 = qr_mgs(Matrix{Float64}(A1))
println("\nQ =")
display(Q1)
println("\n\nR =")
display(R1)
println("\n")

# === Matrix 2 ===
A2 = [2 3 4
      4 5 10
      4 8 2]

println("\nMatrix 2:")
display(A2)
println()
Q2, R2 = qr_mgs(Matrix{Float64}(A2))
println("\nQ =")
display(Q2)
println("\n\nR =")
display(R2)
println()

# === Matrix 3 ===
A3 = [2 -1 4
      0 3 5
      7 2 -6
      1 -4 0]

println("\nMatrix 3:")
display(A3)
println()
Q3, R3 = qr_mgs(Matrix{Float64}(A3))
println("\nQ =")
display(Q3)
println("\n\nR =")
display(R3)
println()


