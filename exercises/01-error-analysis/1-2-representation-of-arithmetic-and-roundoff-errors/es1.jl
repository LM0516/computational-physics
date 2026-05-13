using ComputationalPhysics
using Plots
using LaTeXStrings

# Calculates the sum for a specific N using type T (Float32 or Float64)
function compute_sums(N::Int, ::Type{T}) where {T<:AbstractFloat}
  # (a) Normal ordering: 1 → N
  s_normal = zero(T)
  for n in 1:N
    s_normal += T(1) / (T(n)^2)
  end

  # (b) Reverse ordering: N → 1
  s_reverse = zero(T)
  for n in N:-1:1
    s_reverse += T(1) / (T(n)^2)
  end

  return s_normal, s_reverse
end

function main()
  # Generate an array of N values to study convergence
  N_vals = collect(100:1000:1000000)
  count_N = length(N_vals)

  exact = π^2 / 6

  # Preallocate error arrays
  err32_normal = Vector{Float32}(undef, count_N)
  err32_reverse = Vector{Float32}(undef, count_N)
  err64_normal = Vector{Float64}(undef, count_N)
  err64_reverse = Vector{Float64}(undef, count_N)

  println("Calculating sums for various N")

  for (i, N) in enumerate(N_vals)
    # (c) Single precision (Float32)
    s32_n, s32_r = compute_sums(N, Float32)
    err32_normal[i] = absolute_error(s32_n, Float32(exact))
    err32_reverse[i] = absolute_error(s32_r, Float32(exact))

    # (d) Double precision (Float64)
    s64_n, s64_r = compute_sums(N, Float64)
    err64_normal[i] = absolute_error(s64_n, Float64(exact))
    err64_reverse[i] = absolute_error(s64_r, Float64(exact))
  end

  # Print the final values for the largest N
  println("\n=== SINGLE PRECISION (N = $(N_vals[end])) ===")
  println("Error (normal):  ", err32_normal[end])
  println("Error (reverse): ", err32_reverse[end])

  println("\n=== DOUBLE PRECISION (N = $(N_vals[end])) ===")
  println("Error (normal):  ", err64_normal[end])
  println("Error (reverse): ", err64_reverse[end])

  # Plot convergence for Single Precision
  p1 = plot_generic(N_vals, err32_normal,
    title="Single Precision (Float32)",
    xlabel=L"N", ylabel=L"|S(N) - \pi^2/6|",
    label="Normal (1→N)", xscale=:log10, yscale=:log10)
  plot_add!(p1, N_vals, err32_reverse, label="Reverse (N→1)")
  plot_add!(p1, N_vals, n -> 1.0 ./ n, label=L"O(1/n)", ls=:dash)

  # Plot convergence for Double Precision
  p2 = plot_generic(N_vals, err64_normal,
    title="Double Precision (Float64)",
    xlabel=L"N", ylabel=L"|S(N) - \pi^2/6|",
    label="Normal (1→N)", xscale=:log10, yscale=:log10)
  plot_add!(p2, N_vals, err64_reverse, label="Reverse (N→1)")
  plot_add!(p2, N_vals, n -> 1.0 ./ n, label=L"O(1/n)", ls=:dash)

  # Combine and display
  combined = multi_plot(p1, p2, layout=(1, 2), size=(1000, 800))
  save_plot(combined, "single_and_double_precision", "1-2")
  display(combined)
  readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
  plot_init()
  main()
end
