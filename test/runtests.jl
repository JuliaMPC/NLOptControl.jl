using NLOptControl
using Polynomials
using Base.Test


include("tests.jl")

# first test
actual, approx = test1();
@test_approx_eq_eps(actual,approx,1e-3)


#@testset "Basic Single Interval Tests" begin

#end
