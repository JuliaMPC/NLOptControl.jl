module NLOptMPC

using JuMP
using OrdinaryDiffEq
using DiffEqBase

include("NLOptBase.jl")
using .NLOptBase

export MPC

########################################################################################
# MPC structs
# MPC = Model-Predictive Control
########################################################################################

# ? Initial Plant
mutable struct IP
    control::Control
    state::State
    IP() = new(
        Control(),
        State()
    )
end

# ? Expected Plant
mutable struct EP
    control::Control
    state::State
    EP() = new(
        Control(),
        State()
    )
end

# Model-Predictive Control (MPC) Variables
mutable struct MPCvariables{ T <: Number }
    t::T                        # current simulation time (s)
    tp::Union{T,JuMP.Variable}  # prediction time (if finalTimeDV == true -> this is not known before optimization)
    tex::T                      # execution horizon time
    t0Actual::T                 # actual initial time # TODO: ?
    t0::T                       # mpc initial time # TODO: ?
    tf::T                       # mpc final time # TODO: ?
    t0Param::Any                # parameter for mpc t0  # TODO: ? # ! was Any
    evalNum::Int                # parameter for keeping track of number of MPC evaluations
    goal::Vector{T}             # goal location w.r.t OCP
    goalTol::Vector{T}          # tolerance on goal location
    initOptNum::Int             # number of initial optimization
    previousSolutionNum::Int    # number of times the previous solution should be used
    MPCvariables(T::DataType=Float64) = new{T}(
        convert(T,0.0), # t
        convert(T,0.0), # tp (might be a JuMP.Variable) # ! was Any
        convert(T,0.5), # tex
        convert(T,0.0), # t0Actual
        convert(T,0.0), # t0
        convert(T,0.0), # tf
        Any,            # t0Param
        1,              # evalNum
        Vector{T}(),    # goal
        Vector{T}(),    # goalTol
        3,              # initOptNum
        3               # previousSolutionNum
    )
end

# Model-Predictive Control (MPC)
mutable struct MPC{ T <: Number }
    v::MPCvariables{T}
    ip::IP
    ep::EP
    MPC(T::DataType=Float64) = new{T}(
        MPCvariables(T),
        IP(),
        EP()
    )
end

end # module
