using NLOptControl
using Base.Test

# constants
const tol = 5e-2
const big_tol = 0.5 # can reduce if the number of points are increased
const integrationConfigs = [:lgrExplicit,:lgrImplicit,:trapezoidal,:bkwEuler]
include("ocp.jl")


include("mpc.jl")
