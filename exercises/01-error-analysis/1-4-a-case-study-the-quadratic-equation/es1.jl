using ComputationalPhysics
using Random
using Plots

# === Stability Analysis ===
function generate_stability(offsets, datasets)
    true_var = 1.0
    floor_val = 1e-16  # below F64 machine epsilon, just for display

    clamp_err(e) = max(e, floor_val)

    err_sp32 = [clamp_err(abs(var_single_pass32(Float32.(d)) - true_var) / true_var) for d in datasets]
    err_dp32 = [clamp_err(abs(var_double_pass32(Float32.(d)) - true_var) / true_var) for d in datasets]
    err_sp64 = [clamp_err(abs(var_single_pass64(Float64.(d)) - true_var) / true_var) for d in datasets]
    err_dp64 = [clamp_err(abs(var_double_pass64(Float64.(d)) - true_var) / true_var) for d in datasets]

    plt = plot_generic(offsets, [err_sp32, err_dp32, err_sp64, err_dp64],
        xaxis=:log10, yaxis=:log10,
        label=["Single-Pass F32" "Double-Pass F32" "Single-Pass F64" "Double-Pass F64"],
        title="Numerical Stability: Variance Algorithms",
        xlabel="Data Offset (C)",
        ylabel="Relative Error")

    hline!([1.19e-7], color=:black, ls=:dash, label="F32 Machine Epsilon")
    hline!([2.22e-16], color=:grey, ls=:dot, label="F64 Machine Epsilon")

    return plt
end

function main()
    # === Tests with constant and random vectors ===
    println("=" ^ 60)
    println("PART 1 — Double-pass formula")
    println("=" ^ 60)

    # === Test with [1, 1, 1, 1] — expected variance = 0 ===
    println("\nTest array [1, 1, 1, 1] (expected variance = 0):")
    test_var32 = var_double_pass32(Float32[1.0, 1.0, 1.0, 1.0])
    test_var64 = var_double_pass64(Float64[1.0, 1.0, 1.0, 1.0])
    println("  Variance (32-bit): ", test_var32)
    println("  Variance (64-bit): ", test_var64)

    # === Test with random vectors ===
    Random.seed!(42)
    random_dimension = rand(5:10)
    random_value32 = rand(Float32, random_dimension)
    random_value64 = rand(Float64, random_dimension)

    println("\nRandom Float32 vector (length $random_dimension):")
    display(random_value32)
    println("  Variance (32-bit): ", var_double_pass32(random_value32))

    println("\nRandom Float64 vector (length $random_dimension):")
    display(random_value64)
    println("  Variance (64-bit): ", var_double_pass64(random_value64))

    # === Comparison on datasets with known variance = 1 ===
    println("\n")
    println("=" ^ 60)
    println("PART 2 — Double-pass vs Single-pass comparison")
    println("True variance = 1 for all datasets")
    println("=" ^ 60)

    # Sweep offsets more densely: 10^1 to 10^10
    offsets = [10.0^k for k in 1:0.5:10]
    datasets = [[c, c + 1.0, c + 2.0] for c in offsets]

    for (c, d) in zip(offsets, datasets)
        d32 = Float32.(d)
        d64 = Float64.(d)
        println("\nDataset offset C = $c: ", d)
        println("  Double-pass 32-bit: ", var_double_pass32(d32))
        println("  Double-pass 64-bit: ", var_double_pass64(d64))
        println("  Single-pass 32-bit: ", var_single_pass32(d32))
        println("  Single-pass 64-bit: ", var_single_pass64(d64))
    end

    p = generate_stability(offsets, datasets)

    save_plot(p, "single_and_double_pass", "1-4")
    display(p)

    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
