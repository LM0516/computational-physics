using ComputationalPhysics
using Plots
using LaTeXStrings

function main()
    # === Inline functions ===
    f(x) = (exp(x) - 1) / x
    f_true(x) = expm1(x) / x  # Added expm1 check for verification

    # === Condition Number Analysis ===
    x_vals = LinRange(-1, 1, 100)
    kappa_naive = kappa_f.(x_vals, f)
    max_kappa_naive = maximum(kappa_naive)
    println("=== Part (a): Condition Number Analysis ===")
    println("Maximum condition number: ", max_kappa_naive)

    # Plot condition number
    p1 = plot_generic(x_vals, kappa_naive, title="Condition Number",
            xlabel=L"x", ylabel=L"\kappa_f(x)")

    # === Naive algorithm ===
    x_test = [10.0^(-k) for k in 3:16]
    f_naive = f.(x_test)
    println("\n=== Part (b): Naive Algorithm Results ===")
    for (x, fx) in zip(x_test, f_naive)
        println("f($x) = $fx")
    end

    # === Maclaurin series ===
    n = 5 # number of terms in the series
    f_mclaurin_test = [mclaurin_series(x, n) for x in x_test]
    println("\n=== Part (c): Maclaurin Series Results (n=$n) ===")
    for (x, fx) in zip(x_test, f_mclaurin_test)
        println("f($x) = $fx")
    end

    # === Test multiple values of n ===
    n_values = [3, 5, 7, 9, 11]
    println("\n=== Part (c): Maclaurin Series - Testing different n values ===")
    for n in n_values
        println("\nResults for n=$n:")
        f_mc = [mclaurin_series(x, n) for x in x_test]
        for (x, fx) in zip(x_test[1:3], f_mc[1:3])  # Show first 3 for brevity
            println("  f($x) = $fx")
        end
    end

    # Use a reasonable n for final comparison
    n = 15  # More terms for better accuracy at small x
    # Compute on full range for visualization
    f_mclaurin_vals = [mclaurin_series(x, n) for x in x_vals]
    p2 = plot_generic(x_vals, f_mclaurin_vals, title="Maclaurin Series Results",
            xlabel=L"x", ylabel=L"f(x)", label="Maclaurin Series (n=$n)")

    # === Comparison ===
    f_true_test = f_true.(x_test) # Compute the ground truth

    println("\n=== Part (d): Comparison ===")
    # Formatted to include the expm1 check and compare Maclaurin against it
    println("x\t\tNaive\t\t\tMaclaurin\t\texpm1\t\t\tDiff(Macl - expm1)")
    for (x, fn, fm, ft) in zip(x_test, f_naive, f_mclaurin_test, f_true_test)
        println("$x\t$fn\t$fm\t$ft\t$(abs(fm - ft))")
    end

    # Comparison plot
    p3 = plot_generic(x_test, f_naive, xscale=:log10, label="Naive", 
            xlabel=L"x", ylabel=L"f(x)", title="Comparison of Algorithms", linewidth=1)
    plot_add!(p3, x_test, f_mclaurin_test, label="Maclaurin (n=$n)", linewidth=1)
    # Add the expm1 truth as a dashed line
    plot_add!(p3, x_test, f_true_test, label="expm1 (Truth)", linestyle=:dash, linewidth=1)

    # === Show the results ===
    p = multi_plot(p1, p2, p3, layout=(1, 3), size=(1200, 400))
    save_plot(p, "series_comparison", "1-4")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
