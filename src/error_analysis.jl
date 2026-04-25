"""
#### 1. Type Annotations (`::`)

The syntax:

```julia
x::AbstractVector{T}
```

means that `x` must be a vector-like object whose elements are of type `T`.

* `AbstractVector` allows any vector-like container (e.g. `Vector`, views, etc.)
* `{T}` represents the element type

---

#### 2. Type Parameters (`where {T}`)

The syntax:

```julia
where {T<:AbstractFloat}
```

introduces a type parameter `T` with a constraint.

* `T` is a placeholder for a concrete type (e.g. `Float32`, `Float64`)
* `<:` means “is a subtype of”
* `AbstractFloat` restricts `T` to floating-point types

---

#### 3. Combined Meaning

```julia
function f(x::AbstractVector{T}) where {T<:AbstractFloat}
```

This means:

> The function accepts any vector whose elements are floating-point numbers, 
and refers to that element type as `T`.

---

#### 4. Why Use This Pattern

* **Generic code**: works with `Float32`, `Float64`, `BigFloat`, etc.
* **Performance**: Julia compiles a specialized version for each `T`
* **Flexibility**: supports custom array types and views
* **Type stability**: ensures consistent output types

---

#### 5. Using `T` Inside the Function

The type parameter `T` can be used to keep computations consistent:

```julia
nT = T(n)        # convert integer to same type as elements
one(T)           # multiplicative identity of type T
zero(T)          # additive identity of type T
```

This avoids mixing numeric types (e.g. `Int` and `Float64`).
"""


function mclaurin_series(x::Float64, n::Int)
    sum = 0.0
    for k in 0:n
        sum += x^k / factorial(k + 1)
    end
    return sum
end

"""
    var_double_pass(values) -> variance

Compute the double-pass variance (two loops).
"""
function var_double_pass(x::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(x)
    mean_x = sum(x) / T(n)

    return sum((xi - mean_x)^2 for xi in x) / (n - 1)
end

"""
    var_double_pass(values) -> variance

Compute the single-pass variance (one loop).
"""
function var_single_pass(x::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(x)
    u = sum(xi^2 for xi in x)
    v = sum(x)
    nT = T(n)

    return (u - v^2 / nT) / (nT - one(T))
end

"""
    κ_f(x, f(x)) -> condition_number

Compute the condition number of a function.
"""
function kappa_f(x, func::Function)
    return abs(x * ForwardDiff.derivative(func, x) / func(x))
end

