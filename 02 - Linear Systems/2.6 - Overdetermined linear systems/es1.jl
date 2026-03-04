include("../../modules/linear_systems.jl")
using .LinearSystems
using LinearAlgebra

function main()
    A = [
        2 -1
        0 1
        -2 2
    ]
    b = [
        1
        -5
        6
    ]
    x = least_squares(A, b)
    println(x)
    println(A \ b)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
