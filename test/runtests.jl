using NLOptControl
using JuMP
using Base.Test

# constants
const tol=1e-2;
const integrationConfigs=[[:ps,:lgrExplicit],[:tm,:bkwEuler],[:tm,:trapezoidal]]

include(Pkg.dir("NLOptControl/test/ocp.jl"))
