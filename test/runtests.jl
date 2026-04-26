using Test
using ComputationalPhysics
using Random
using LinearAlgebra
using SparseArrays
using Statistics

Random.seed!(1234)

include(joinpath(@__DIR__, "submodules", "test_utils.jl"))
include(joinpath(@__DIR__, "submodules", "error_analysis_tests.jl"))
include(joinpath(@__DIR__, "submodules", "interpolations_tests.jl"))
include(joinpath(@__DIR__, "submodules", "linear_systems_tests.jl"))
include(joinpath(@__DIR__, "submodules", "nonlinear_equations_tests.jl"))
include(joinpath(@__DIR__, "submodules", "numerical_integration_tests.jl"))
include(joinpath(@__DIR__, "submodules", "ode_tests.jl"))
include(joinpath(@__DIR__, "submodules", "schrodinger_tests.jl"))
