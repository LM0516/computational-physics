# === Tests with constant and random vectors ===
println("="^60)
println("PART 1 — Double-pass formula")
println("="^60)

test_array = [1.0, 1.0, 1.0, 1.0]

# === Test with [1, 1, 1, 1] — expected variance = 0 ===
println("\nTest array [1, 1, 1, 1] (expected variance = 0):")
test_var32 = var_double_pass(Vector{Float32}(test_array))
test_var64 = var_double_pass(Vector{Float64}(test_array))
println("  Variance (32-bit): ", test_var32)
println("  Variance (64-bit): ", test_var64)

# === Test with random vectors ===
Random.seed!(42)
random_dimension = rand(5:10)
random_value32 = rand(Float32, random_dimension)
random_value64 = rand(Float64, random_dimension)

println("\nRandom Float32 vector (length $random_dimension):")
display(random_value32)
println("  Variance (32-bit): ", var_double_pass(random_value32))

println("\nRandom Float64 vector (length $random_dimension):")
display(random_value64)
println("  Variance (64-bit): ", var_double_pass(random_value64))

# === Comparison on datasets with known variance = 1 ===
println("\n")
println("="^60)
println("PART 2 — Double-pass vs Single-pass comparison")
println("True variance = 1 for all datasets")
println("="^60)

# Sweep offsets more densely: 10^1 to 10^10
offsets = [10.0^k for k in 1:0.5:10]
datasets = [[c, c + 1.0, c + 2.0] for c in offsets]

for (c, d) in zip(offsets, datasets)
  d32 = Float32.(d)
  d64 = Float64.(d)
  println("\nDataset offset C = $c: ", d)
  println("  Double-pass 32-bit: ", var_double_pass(d32))
  println("  Double-pass 64-bit: ", var_double_pass(d64))
end
