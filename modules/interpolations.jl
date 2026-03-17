module Interpolations
using Distributions
using LinearAlgebra
export barycentric_lagrange, chi_square, p_value, fit_goodness

const ∏ = prod
const ∑ = sum

# NOTE: For better performance look 'dispatch', there is an if that is called to
# many times.
# TODO: Check memory allocations
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
function chi_square(o::AbstractVector, e::AbstractVector)
    # eachindex(o, e) ensures both vectors have the exact same length, 
    # throwing an error automatically if they don't.
    return ∑(eachindex(o, e); init=0.0) do i
        ei = e[i]
        # Only compute if the expected value is not near zero
        abs(ei) > 1e-10 ? (o[i] - ei)^2 / ei : 0.0
    end
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

"""
Combined version of `chi_square` and `p_value` functions.
"""
function fit_goodness(y::Vector{Float64}, y_true::Vector{Float64}, coeffs::Vector{Float64})
    dof = length(y) - length(coeffs)
    chi2 = chi_square(y, y_true)
    chi2r = chi2 / dof 
    p = p_value(chi2, dof)

    println("Chi-square: ", chi2)
    println("Degrees of freedom: ", dof)
    println("Reduced chi-square: ", chi2r)
    println("P-value: ", p)
    println()
    return chi2, dof, chi2r, p
end

end # module Interpolations
