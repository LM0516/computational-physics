module Interpolations
using Distributions
using LinearAlgebra
export barycentric_lagrange, chi_square, p_value

const ∏ = prod
const ∑ = sum

"""
Evaluate the Lagrange interpolating polynomial at a given point `x_eval` using the
barycentric interpolation formula.
"""
function barycentric_lagrange(x_eval::Float64, x_nodes::Vector{Float64}, f::Function; type="lagrange")
  n = length(x_nodes)

  # Compute function values at nodes
  f_values = [f(x) for x in x_nodes]

  # Compute barycentric weights λ
  λ = zeros(n)
  if type == "lagrange"
    for j in 1:n
      # Product of (x[j] - x[i]) for all i ≠ j
      λ[j] = 1 / ∏(x_nodes[j] - x_nodes[i] for i in 1:n if i != j)
    end
  elseif type == "chebyshev"
    for j in 1:n
      if j == 1 || j == n
        λ[j] = (-1)^j * 1 / 2
      else
        λ[j] = (-1)^j
      end
    end
  end

  # Check if x_eval is one of the nodes
  for j in 1:n
    if x_eval ≈ x_nodes[j]
      return f_values[j]
    end
  end

  # Compute interpolation using barycentric formula
  numerator = ∑(λ[j] * f_values[j] / (x_eval - x_nodes[j]) for j in 1:n)
  denominator = ∑(λ[j] / (x_eval - x_nodes[j]) for j in 1:n)

  p = numerator / denominator
  return p
end

"""
Chi-square test for data fit.

Interpretation of reduced chi-square:
- ≈ 0 → good fit
- >> 1 → poor fit (model doesn't describe data well)
"""
function chi_square(o::Vector{Float64}, e::Vector{Float64})
    return ∑((o .- e).^2 ./ e)
end

"""
p-value test for data fit.

Interpretation:
- p > 0.05 → fit is acceptable (fail to reject null hypothesis that model fits data)
- p < 0.05 → fit is poor (model likely doesn't describe data well)
- p < 0.01 → fit is very poor
"""
function p_value(chi2::Float64, dof::Int)
    dist = Chisq(dof)
    p = 1 - cdf(dist, chi2)
    return p
end

# TODO: Implement a t-test (z-test) to compare real results to calculated results.
function t_test()
end

end # module Interpolations
