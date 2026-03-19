using ComputationalPhysics
using Plots
using LaTeXStrings

function somma(a, N::Int)
  s1 = zeros(eltype(a), N)  # Normal ordering: 1 → N
  s2 = zeros(eltype(a), N)  # Reverse ordering: N → 1

  # Normal Ordering
  s1[1] = a[1]
  for n = 2:N
    s1[n] = s1[n-1] + a[n]
  end

  # Reverse Ordering
  s2[1] = a[N]
  for n = 2:N
    s2[n] = s2[n-1] + a[N-n+1]
  end

  exact = π^2 / 6
  err1 = @. abs(s1 - exact)
  err2 = @. abs(s2 - exact)

  return s1, s2, err1, err2
end

function main()
    N = 10000

    # Single precision
    a32 = Float32[1 / n^2 for n = 1:N]
    s32_normal, s32_reverse, err32_normal, err32_reverse = somma(a32, N)

    println("=== SINGLE PRECISION ===")
    println("(a) Normal ordering (1→N):  S($N) = ", s32_normal[end])
    println("(b) Reverse ordering (N→1): S($N) = ", s32_reverse[end])
    println("    Exact value: π²/6 = ", π^2 / 6)
    println("    Error (normal):  ", err32_normal[end])
    println("    Error (reverse): ", err32_reverse[end])

    # Double precision
    a64 = Float64[1 / n^2 for n = 1:N]
    s64_normal, s64_reverse, err64_normal, err64_reverse = somma(a64, N)

    println("\n=== DOUBLE PRECISION ===")
    println("(d) Normal ordering (1→N):  S($N) = ", s64_normal[end])
    println("(d) Reverse ordering (N→1): S($N) = ", s64_reverse[end])
    println("    Error (normal):  ", err64_normal[end])
    println("    Error (reverse): ", err64_reverse[end])

    # Plot convergence for single precision
    p1 = plot_generic(1:N, err32_normal,
        title="Single Precision",
        xlabel=L"N", ylabel=L"|S(N) - \pi^2/6|",
        label="Normal ordering (1→N)", yscale=:log10, color=1)
    plot_add!(p1, 1:N, err32_reverse, label="Reverse ordering (N→1)", color=2)
    plot_add!(p1, 1:N, n -> 1.0 / n, label=L"O(1/n)", ls=:dash, color=:black, lw=1)

    # Plot convergence for double precision
    p2 = plot_generic(1:N, err64_normal,
    title="Double Precision",
    xlabel=L"N", ylabel=L"|S(N) - \pi^2/6|",
    label="Normal ordering (1→N)", yscale=:log10, color=1)
    plot_add!(p2, 1:N, err64_reverse, label="Reverse ordering (N→1)", color=2)
    plot_add!(p2, 1:N, n -> 1.0 / n, label=L"O(1/n)", ls=:dash, color=:black, lw=1)

    # Combine and save
    combined = multi_plot(p1, p2, plot_title="Convergence Error", layout=(1, 2), size=(1000, 800))
    save_plot(combined, "single_and_double_precision", "1-2")
    display(combined)

    readline()

end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
