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

years = 30
initial_capital = 1.71828182845904523536028747135

final_capital_32 = investment(initial_capital, years, Float32)
final_capital_64 = investment(initial_capital, years, Float64)
final_capital_BF = investment(initial_capital, years, BigFloat)

# yearly_yeld = (final_capital[years] / initial_capital)^(1 / years) - 1
# yeld = (final_capital[years] - initial_capital) / initial_capital * 100

# println("Initial capital: ", initial_capital)
# println("Final capital: ", final_capital[years])
# println("All time yield: ", yeld, "%")
# println("Yearly yield: ", yearly_yeld * 100, "%")

p = plot(final_capital_32, label="Float32")
plot!(final_capital_64, label="Float64")
plot!(final_capital_BF, label="BigFloat")
display(p)
readline()
