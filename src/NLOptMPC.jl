module NLOptMPC

using JuMP
using OrdinaryDiffEq
using DiffEqBase
using Parameters

include("NLOptBase.jl")
using .NLOptBase

export MPC

########################################################################################
# MPC structs
# MPC = Model-Predictive Control
########################################################################################

# ? Initial Plant
@with_kw mutable struct IP
    control::Control    = Control() #
    state::State        = State()   #
end

# ? Expected Plant
@with_kw mutable struct EP
    control::Control    = Control() #
    state::State        = State()   #
end

# Model-Predictive Control (MPC) Variables
@with_kw mutable struct MPCvariables{ T <: Number }
    t::T                        = convert(T,0.0)    # current simulation time (s)
    tp::Union{T,JuMP.Variable}  = convert(T,0.0)    # prediction time (if finalTimeDV == true -> this is not known before optimization)
    tex::T                      = convert(T,0.5)    # execution horizon time
    t0Actual::T                 = convert(T,0.0)    # actual initial time # TODO: ?
    t0::T                       = convert(T,0.0)    # mpc initial time # TODO: ?
    tf::T                       = convert(T,0.0)    # mpc final time # TODO: ?
    t0Param::Any                = Any               # parameter for mpc t0  # TODO: ? # ! was Any
    evalNum::Int                = 1                 # parameter for keeping track of number of MPC evaluations
    goal::Vector{T}             = Vector{T}()       # goal location w.r.t OCP
    goalTol::Vector{T}          = Vector{T}()       # tolerance on goal location
    initOptNum::Int             = 3                 # number of initial optimization
    previousSolutionNum::Int    = 3                 # number of times the previous solution should be used
end

# Model-Predictive Control (MPC)
@with_kw mutable struct MPC{ T <: Number }
    v::MPCvariables{T}  = MPCvariables{T}() #
    ip::IP              = IP()              #
    ep::EP              = EP()              #
end

end # module
