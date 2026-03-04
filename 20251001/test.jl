y = zeros(4)

x = 1:4

println(x)

for i = 1:4
    y[i] =  exp(i)
end

println("Risultato: ", y)
