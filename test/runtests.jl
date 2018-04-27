using NLOptControl
using Base.Test

# NOTE: for some reason if a test fails, then later tests are effected.
# This behaviour can be very perplexing.
# If a test fails ignore all remaining errors they may be erroneous!
# Either fix the broken test or put it at the bottom of the stack.
const tol = 5e-2
const big_tol = 0.5 # can reduce if the number of points are increased
const integrationConfigs = [:lgrExplicit,:lgrImplicit,:trapezoidal,:bkwEuler]
include("ocp.jl")
include("mpc.jl")
