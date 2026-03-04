module Interpolations
using LinearAlgebra
export barycentric_lagrange

const ∏ = prod
const ∑ = sum

"""
    barycentric_lagrange(x_eval::Float64, x_nodes::Vector{Float64}, f::Function) -> Float64

Evaluate the **Lagrange interpolating polynomial** at a given point `x_eval` using the 
**barycentric interpolation formula**.

# Arguments
- `x_eval::Float64`:  
  The point at which the interpolating polynomial is to be evaluated.
  
- `x_nodes::Vector{Float64}`:  
  A vector containing the interpolation nodes (distinct x-values).

- `f::Function`:  
  The function to be interpolated. It should accept a single `Float64` argument.

# Returns
- `p::Float64`:  
  The interpolated value `p = P(x_eval)` of the Lagrange polynomial that passes through 
  the points `(x_nodes[i], f(x_nodes[i]))`.
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

end # module Interpolations
