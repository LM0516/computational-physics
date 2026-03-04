files = [
    "01 - Error Analysis/1.2 - Representation of Arithmetic and Roundoff Errors/es1.jl",
    #="01 - Error Analysis/1.2 - Representation of Arithmetic and Roundoff Errors/es2.jl",=#
    #="01 - Error Analysis/1.2 - Representation of Arithmetic and Roundoff Errors/es3.jl",=#
    "01 - Error Analysis/1.4 - A case study: the quadratic equation/es1.jl",
    "01 - Error Analysis/1.4 - A case study: the quadratic equation/es2.jl",
    "02 - Linear Systems/2.1 - Triangular systems/es1.jl",
    "02 - Linear Systems/2.3 - LU decomposition/es2.jl",
    "02 - Linear Systems/2.3 - LU decomposition/es3.jl",
    # "02 - Linear Systems/2.3 - LU decomposition/es4.jl",
    "02 - Linear Systems/2.3 - LU decomposition/es6.jl",
    "02 - Linear Systems/2.5 - Conditioning of linear systems/es1.jl",
    "02 - Linear Systems/2.6 - Overdetermined linear systems/es1.jl",
    "02 - Linear Systems/2.6 - Overdetermined linear systems/es2.jl",
    "02 - Linear Systems/2.6 - Overdetermined linear systems/es3.jl",
    "02 - Linear Systems/2.7 - QR decomposition/es2.jl",
    "02 - Linear Systems/2.7 - QR decomposition/es3.jl",
    "02 - Linear Systems/2.7 - QR decomposition/es4.jl",
    "03 - Interpolations/3.1 - Polynomial interpolation/es1.jl",
    "03 - Interpolations/3.2 - Lagrande-Waring interpolation/es2.jl",
    "03 - Interpolations/3.4 - Error on interpolation/es2.jl",
    "03 - Interpolations/3.4 - Error on interpolation/es3.jl",
    "03 - Interpolations/3.4 - Error on interpolation/es4.jl",
    "03 - Interpolations/3.4 - Error on interpolation/es5.jl",
    "04 - Roots of Nonlinear Equations/4.2 - Bisection method/es2.jl",
    "04 - Roots of Nonlinear Equations/4.3 - Newton method/es2.jl",
    "04 - Roots of Nonlinear Equations/4.3 - Newton method/es3.jl",
    "04 - Roots of Nonlinear Equations/4.3 - Newton method/es4.jl",
    "04 - Roots of Nonlinear Equations/4.4 - Interpolation methods/es2.jl",
    "04 - Roots of Nonlinear Equations/4.4 - Interpolation methods/es3.jl",
    "05 - Numerical Integration/5.1 - Newton-Cotes formulas/es2 backup 1.jl",
    "05 - Numerical Integration/5.1 - Newton-Cotes formulas/es2.jl",
    "05 - Numerical Integration/5.2 - Free-nodes integration/es1.jl",
    "05 - Numerical Integration/5.3 - Advanced topics in integration/es1.jl",
    "05 - Numerical Integration/5.3 - Advanced topics in integration/es3.jl",
    "05 - Numerical Integration/5.3 - Advanced topics in integration/es4.jl",
    "06 - Ordinary Differential Equations/6.2 - Basics of one-step methods/es2.jl",
    "06 - Ordinary Differential Equations/6.2 - Basics of one-step methods/es3.jl",
    "06 - Ordinary Differential Equations/6.3 - Runge-Kutta methods/es1.jl",
    "06 - Ordinary Differential Equations/6.3 - Runge-Kutta methods/es2.jl",
    "06 - Ordinary Differential Equations/6.3 - Runge-Kutta methods/es3.jl",
    "06 - Ordinary Differential Equations/6.5 - Some interesting equation/es1.jl",
    "06 - Ordinary Differential Equations/6.5 - Some interesting equation/es2.jl",
    "06 - Ordinary Differential Equations/6.5 - Some interesting equation/es3.jl",
    "06 - Ordinary Differential Equations/6.5 - Some interesting equation/es5.jl"
]

for file in files
    println("="^60)
    println("Running $file")
    include(file)
    println("Press Enter to run nex file...")
    readline()
end
