using ComputationalPhysics
using Plots
using LaTeXStrings

function investment(initial_capital, years::Int, T::Type)
  capital = zeros(T, years)

  capital[1] = T(initial_capital) - 1

  for i = 2:years
    capital[i] = capital[i-1] * i - 1
  end

  return capital
end

function main()
    years = 30
    initial_capital = 1.71828182845904523536028747135

    final_capital_32 = investment(initial_capital, years, Float32)
    final_capital_64 = investment(initial_capital, years, Float64)
    final_capital_BF = investment(initial_capital, years, BigFloat)

    p = plot_generic(1:length(final_capital_32), final_capital_32, label="Float32")
    plot_add!(p, 1:length(final_capital_64), final_capital_64, label="Float64")
    plot_add!(p, 1:length(final_capital_BF), final_capital_BF, label="BigFloat")
    display(p)
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_init()
    main()
end
