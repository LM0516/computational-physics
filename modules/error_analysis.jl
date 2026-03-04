module ErrorAnalysis

function mclaurin_series(x::Float64, n::Int)
    sum = 0.0
    for k in 0:n
        sum += x^k / factorial(k + 1)
    end
    return sum
end

end
