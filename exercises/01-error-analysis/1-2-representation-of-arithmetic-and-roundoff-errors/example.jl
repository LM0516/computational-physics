println(eps(Float64))
@show floatmin(),floatmax();

e = eps()/2
println((1.0 + e) - 1.0)
println(1.0 + (e - 1.0))
