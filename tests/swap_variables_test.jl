using BenchmarkTools

function three_var(x1, x2, iter)
    for i in 1:iter
        x_new = x1 + x2
        x1 = x2
        x2 = x_new
    end
    return x2
end

function two_var(x1, x2, iter)
    for _ in 1:iter
        # This creates a tuple and unpacks it
        x1, x2 = x2, x1 + x2
    end
    return x2
end

x1 = 1
x2 = 2
iter = 100000 # Increased iterations to make the work measurable

# We use '$' to interpolate variables so we don't benchmark 
# the cost of looking up global variables.
println("--- Three Variables ---")
@btime three_var($x1, $x2, $iter)

println("\n--- Two Variables ---")
@btime two_var($x1, $x2, $iter)
