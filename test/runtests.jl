using NLOptControl
using Base.Test

# constants
const tol=5e-2;
const integrationConfigs=[:lgrImplicit,:lgrExplicit,:trapezoidal,:bkwEuler]

include("ocp.jl")

#using NBInclude
#nbinclude(joinpath(dirname(@__FILE__), "..", "examples", "tutorial.ipynb"))
