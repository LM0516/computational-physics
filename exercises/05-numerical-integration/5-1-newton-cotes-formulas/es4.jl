using ComputationalPhysics
using Plots
using LaTeXStrings

function solutions(a::Real, b::Real, f::Function, nodes, label, exact; rule::String="trapezoidal")
  trapezoidal = Vector{Float64}(undef, length(nodes))
  simpson = Vector{Float64}(undef, length(nodes))
  n_nodes = length(nodes)

  # Compute solutions 
  for j in 1:n_nodes
    trapezoidal[j] = composite_trapezoidal(a, b, f, Int(nodes[j]))
  end

  for j in 1:n_nodes
    simpson[j] = composite_simpson(a, b, f, Int(nodes[j]))
  end

  # Plot solutions
  p = scatter_generic(nodes, simpson, label=label)
  hline!([exact], label="Exact solution (Simpson rule): $(round(exact, digits=4))", linestyle=:dash)
  xlabel!("Nodes")
  ylabel!("Solution")

  return trapezoidal, simpson, p
end

function errors(s::Vector{Float64}, exact::Float64, nodes)
  println("-"^40)

  absolute_errors = @. abs(s - exact)

  # Plot absolute errors
  absolute_plots = plot(make_log_safe(nodes), make_log_safe(absolute_errors), label="Absolute errors", yscale=:log10, xscale=:log10)
  xlabel!("Nodes")
  ylabel!("Absolute error")

  return absolute_plots
end

function main()
  k = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  nodes = @. 10 * 2^k
  rule = "trapezoidal"

  println("=== $rule rule for integration ===")
  println("-"^40)

  println("=== First integral: x log(1 + x) ===")
  f1 = x -> x * log(1 + x)
  a1 = 0.0
  b1 = 1.0
  exact1 = 1 / 4
  trapezoidal1, simpson1, p1 = solutions(a1, b1, f1, nodes, L"x log(1 + x)", exact1, rule=rule)
  println("Trapezoidal solution: $(trapezoidal1[end])")
  println("Simpson solution: $(simpson1[end])")
  ap1 = errors(simpson1, exact1, nodes)

  println("=== Second integral: x^2 tan^-1(x) ===")
  f2 = x -> x^2 * atan(x)
  a2 = 0.0
  b2 = 1.0
  exact2 = (π - 2 + 2log(2)) / 12
  trapezoidal2, simpson2, p2 = solutions(a2, b2, f2, nodes, L"x^2 tan^{-1}(x)", exact2, rule=rule)
  println("Trapezoidal solution: $(trapezoidal2[end])")
  println("Simpson solution: $(simpson2[end])")
  ap2 = errors(simpson2, exact1, nodes)

  println("=== Third integral: e^x cos(x) ===")
  f3 = x -> exp(x) * cos(x)
  a3 = 0.0
  b3 = π / 2
  exact3 = (exp(π / 2) - 1) / 2
  trapezoidal3, simpson3, p3 = solutions(a3, b3, f3, nodes, L"e^x cos(x)", exact3, rule=rule)
  println("Trapezoidal solution: $(trapezoidal3[end])")
  println("Simpson solution: $(simpson3[end])")
  ap3 = errors(simpson3, exact1, nodes)

  println("=== Fourth integral: (tan^-1 (√2 + x^2))/(1 + x^2)(√2 + x^2) ===")
  f4 = x -> atan(sqrt(2 + x^2)) / ((1 + x^2) * sqrt(2 + x^2))
  a4 = 0.0
  b4 = 1.0
  exact4 = 5π^2 / 96
  trapezoidal4, simpson4, p4 = solutions(a4, b4, f4, nodes, L"\frac{\tan^{-1} (\sqrt{2 + x^2})}{(1 + x^2)\sqrt{2 + x^2}}", exact4, rule=rule)
  println("Trapezoidal solution: $(trapezoidal4[end])")
  println("Simpson solution: $(simpson4[end])")
  ap4 = errors(simpson4, exact1, nodes)

  println("=== Fifth integral: √x log(x) ===")
  f5 = x -> sqrt(x) * log(x)
  a5 = 0.0
  b5 = 1.0
  exact5 = -4 / 9
  trapezoidal5, simpson5, p5 = solutions(a5, b5, f5, nodes, L"\sqrt{x} log(x)", exact5, rule=rule)
  println("Trapezoidal solution: $(trapezoidal5[end])")
  println("Simpson solution: $(simpson5)")
  ap5 = errors(simpson5, exact1, nodes)

  println("=== Sixth integral: √(1 - x^2) ===")
  f6 = x -> sqrt(1 - x^2)
  a6 = 0.0
  b6 = 1.0
  exact6 = π / 4
  trapezoidal6, simpson6, p6 = solutions(a6, b6, f6, nodes, L"\sqrt{1 - x^2}", exact6, rule=rule)
  println("Trapezoidal solution: $(trapezoidal6[end])")
  println("Simpson solution: $(simpson6[end])")
  ap6 = errors(simpson6, exact1, nodes)

  # Show the results
  p_results = multi_plot(p1, p2, p3, p4, p5, p6, layout=(2, 3), size=(1200, 800))
  save_plot(p_results, "newton-cotes-comparison", "5-1")
  display(p_results)
  readline()

  # Show the absolute errors
  p_absolute = multi_plot(ap1, ap2, ap3, ap4, ap5, ap6, layout=(2, 3), size=(1200, 800))
  save_plot(p_absolute, "newton-cotes-comparison-absolute-error", "5-1")
  display(p_absolute)
  readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
  plot_init()
  main()
end

