using ComputationalPhysics
using Plots
using LaTeXStrings

function condition_number_analysis(x_vals, f)
    # Corrected to loop over the array of x values using your library's kappa_f
    κ_f = [kappa_f(x, f) for x in x_vals]
    max_κ = maximum(κ_f)
    println("Maximum condition number: ", max_κ)
    return κ_f
end

function find_optimal_n(x::Float64; tol=eps(Float64))
    # To find the best n I have to find the next term in the series that is smaller
    # than the precision output of the floating-point type.
    # In this case I have a max x value of 1 so I have to insert that value
    # in this function. I'm expecting to have n = 16 to maintain full precision

    n = 1
    while true
        # The next term in your specific Maclaurin series
        next_term = x^n / factorial(n + 1)
        if abs(next_term) < tol
            return n - 1
        end
        n += 1
        if n > 20
            break
        end
    end
    return n
end

function main()
    # === Inline functions ===
    f(x) = (exp(x) - 1) / x
    f_true(x) = expm1(x) / x

    # === Naive algorithm ===
    x_test = [10.0^(-k) for k in 3:16]
    f_naive = f.(x_test)
    f_true_test = f_true.(x_test)

    println("\n=== Naive Algorithm Results ===")
    for (x, fx) in zip(x_test, f_naive)
        println("f($x) = $fx")
    end

    # === Maclaurin series ===
    println("\n=== Find optimal Maclaurin order ===")
    best_n = find_optimal_n(1.0)
    println("Best n for x=1 is: ", best_n)

    f_mclaurin_test = [mclaurin_series(Float64(x), best_n) for x in x_test]
    println("\n=== Maclaurin Series Results (n=$best_n) ===")
    for (x, fx) in zip(x_test, f_mclaurin_test)
        println("f($x) = $fx")
    end

    # === Test multiple values of n ===
    n_values = [3, 5, 7, 9, 11]
    println("\n=== Maclaurin Series - Testing different n values ===")
    for n in n_values
        println("\nResults for n=$n:")
        f_mc = [mclaurin_series(Float64(x), n) for x in x_test]
        for (x, fx) in zip(x_test[1:3], f_mc[1:3])  # Show first 3 for brevity
            println("  f($x) = $fx")
        end
    end

    # === Condition Number Analysis ===
    x_vals = range(-1.0, 1.0, length=100)

    println("\n=== Condition Number Analysis ===")
    κ_naive = condition_number_analysis(x_vals, f)
    κ_true = condition_number_analysis(x_vals, f_true)

    p1 = plot_generic(x_vals, κ_naive,
        xlabel=L"x", ylabel=L"\kappa_f(x)", label="Naive results")
    plot_add!(p1, x_vals, κ_true, label="expm1 results", linestyle=:dash)

    # === Comparison ===
    println("\n=== Comparison ===")
    println("Naive VS Maclaurin (n=$best_n): ")
    @show f_mclaurin_test - f_naive
    println("\nexpm1 VS Maclaurin (n=$best_n): ")
    @show f_true_test - f_mclaurin_test

    p2 = plot_generic(x_test, f_naive, xscale=:log10, label="Naive",
        xlabel=L"x", ylabel=L"f(x)", linewidth=1)
    plot_add!(p2, x_test, f_mclaurin_test, label="Maclaurin (n=$best_n)", linewidth=1)
    plot_add!(p2, x_test, f_true_test, label="expm1 (Truth)", linestyle=:dash, linewidth=1)

    # === Show the results ===
    p = multi_plot(p1, p2, layout=(1, 2), size=(1200, 400))
    save_plot(p, "series_comparison", "1-4")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
