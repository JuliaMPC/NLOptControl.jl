using NLOptControl
using Test
import Statistics: mean

const tol = 5e-2
const big_tol = 0.5 # can reduce if the number of points are increased
#const integrationConfigs = [:lgrExplicit,:lgrImplicit,:trapezoidal,:bkwEuler]
const integrationConfigs = [:lgrExplicit,:trapezoidal,:bkwEuler]
const integrationConfigs_no_trap = [:lgrExplicit,:bkwEuler] #  TODO figure out why this is broken with trapezoidal: using new names for states and controls; also naming them x and u to check for bugs

@testset "OCP" begin
    include("ocp.jl")
end
@testset "MPC" begin
    include("mpc.jl")
end
