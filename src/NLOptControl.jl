isdefined(Base, :__precompile__) && __precompile__()

module NLOptControl
#TODO  enable setvalue() functionality

using JuMP
import JuMP: setRHS, 
             getvalue, 
             setvalue, 
             @NLexpression, 
             @NLobjective, 
             @NLparameter, 
             @NLconstraint, 
             internalmodel
using Ipopt
using FastGaussQuadrature
using DataFrames
using CSV
using Interpolations
import Printf
import LinearAlgebra

include("NLOptBase.jl")
using .NLOptBase
export  State,
        Control,
        Constraint,
        Results,
        Settings

include("NLOptMPC.jl")
using .NLOptMPC
export MPC

# Optimal Control Problem Structure
mutable struct OCP{T <: Number}
  # general properties
  state::State                      # state data
  control::Control                  # control data
  tf::T                             # final time
  t0::T                             # initial time # TODO: consider getting rid of this or replacing it with `n.mpc.v.t0Param`
  tV::Vector{T}                     # vector for use with time varying constraints

  # Boundary conditions
  X0::Vector{T}                     # initial state conditions
  X0_tol::Vector{T}                 # initial state tolerance
  x0s::Vector{JuMP.Variable}        # initial state variables
  XF::Vector{T}                     # final state conditions
  XF_tol::Vector{T}                 # final state tolerance
  xFs::Vector{JuMP.Variable}        # final state variables

  # Constant bounds on state variables
  XL::Vector{T}                     # Constant lower bound on state variables
  XU::Vector{T}                     # Constant upper bound on state variables

  # Variables for linear bounds on state variables
  mXL::Vector{T}                    # slope on XL -> time always starts at zero
  mXU::Vector{T}                    # slope on XU -> time always starts at zero
  XL_var::JuMP.Variable             # time varying lower bounds on states # ! NOTE: not used currently
  XU_var::JuMP.Variable             # time varying upper bounds on states # ! NOTE: not used currently

  # Constant bounds on control variables
  CL::Vector{T}                     # Constant lower bound on control variables
  CU::Vector{T}                     # Constant upper bound on control variables

  # Pseudospectral method data
  Nck::Vector{Int}                  # number of collocation points per interval
  Nck_cum::Vector{Int}              # cumulative number of points per interval
  Nck_full::Vector{Int}             # [0;cumsum(n.ocp.Nck+1)]
  Ni::Int                           # number of intervals
  tau::Matrix{T}                    # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
  ts::Matrix{T}                     # time scaled based off of tau
  w::Matrix{T}                      # weights
  ws::Matrix{T}                     # scaled weights
  DMatrix::Vector{Matrix{T}}        # differention matrix
  IMatrix::Vector{Matrix{T}}        # integration matrix

  # tm method data
  N::Int                            # number of points in discretization
  dt::Vector{T}                     # array of dts

  mdl::JuMP.Model                   # JuMP model
  params                            # parameters for the models
  DXexpr                            # ? DX expression
  NLcon                             # ! NOTE: not used currently

  # Scaling factors
  XS::Vector{T}                    # scaling factors on states
  CS::Vector{T}                    # scaling factors on controls

    OCP(T::DataType=Float64) = new{T}(
        State(),                              # state data
        Control(),                            # control data
        convert(T, 0),                        # final time
        @NLparameter(JuMP.Model(), x == 0),   # initial time # TODO: consider getting ride of this! or replace it with n.mpc.v.t0Param
        Vector{T}(),                          # vector for use with time varying constraints
        Vector{T}(),                          # initial state conditions
        Vector{T}(),                          # initial state tolerance
        Vector{JuMP.Variable}(),              # initial state variables
        Vector{T}(),                          # final state conditions
        Vector{T}(),                          # final state tolerance
        Vector{JuMP.Variable}(),              # final state variables
        Vector{T}(),                          # Constant lower bound on state variables
        Vector{T}(),                          # Constant upper bound on state variables
        Vector{T}(),                          # slope on XL -> time always starts at zero
        Vector{T}(),                          # slope on XU -> time always starts at zero
        Matrix{T}(undef, 1, 1),               # time varying lower bounds on states # ! NOTE: not used currently
        Matrix{T}(undef, 1, 1),               # time varying upper bounds on states # ! NOTE: not used currently
        Vector{T}(),                          # Constant lower bound on control variables
        Vector{T}(),                          # Constant upper bound on control variables
        Vector{Int}(),                        # number of collocation points per interval
        Vector{Int}(),                        # cumulative number of points per interval
        Vector{Int}(),                        # [0;cumsum(n.ocp.Nck+1)]
        Int(0),                               # number of intervals
        Matrix{T}(undef, 1, 1),               # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
        Matrix{T}(undef, 1, 1),               # time scaled based off of tau
        Matrix{T}(undef, 1, 1),               # weights
        Matrix{T}(undef, 1, 1),               # scaled weights
        Vector{Matrix{T}}(),                  # differention matrix
        Vector{Matrix{T}}(),                  # integration matrix
        Int(0),                               # number of points in discretization
        Vector{T}(),                          # array of dts
        JuMP.Model(),                         # JuMP model
        Any[],                                # Jump model parameters
        Any[],                                # ? DX expression
        Any[],                                # ! Unused - Nonlinear conditions
        Vector{T}(),                          # scaling factors on states
        Vector{T}(),                          # scaling factors on controls
    )
end

mutable struct OCPFlags
  defined::Bool  # a bool to tell if define() has been called
  OCPFlags() = new(
    false
)
end

mutable struct MPCFlags
    defined::Bool
    goalReached::Bool
    simFailed::Vector{Union{Bool, Symbol}}   # a bool to indicate that the simulation failed and a symbol to indicate why
    ipDefined::Bool
    epDefined::Bool
    MPCFlags() = new(
    false,
    false,
    [false, :NaN],
    false,
    false
)
end

mutable struct Flags
    ocp::OCPFlags
    mpc::MPCFlags
    Flags() = new(
    OCPFlags(),
    MPCFlags()
)
end

abstract type AbstractNLOpt end
mutable struct NLOpt{T <: Number} <: AbstractNLOpt
    # major data structs
    ocp::OCP{T}
    mpc::MPC{T}
    s::Settings
    r::Results{T}
    f::Flags
    NLOpt(T::DataType=Float64) = new{T}(
    OCP(T),
    MPC(T),
    Settings(),
    Results(T),
    Flags()
)
end

export  NLOpt,
        OCP


include("NLOptBaseutils.jl")
export  resultsDir!,
        evalConstraints!,
        postProcess!,
        optimize!,
        interpolateLagrange!,
        interpolateLinear!,
        interpolate_lagrange,
        opt2dfs!
# TODO: make a hamiltonian function

include("NLOptMPCutils.jl")
export  defineMPC!,
        initOpt!,
        defineIP!,
        mapNames!,
        simIPlant!,
        updateX0!,
        currentIPState,
        goalReached!,
        simMPC!,
        plant2dfs!,
        predictX0!

include("utils.jl")
export  defineTolerances!,
        linearStateTolerances!,
        newConstraint!,
        evalMaxDualInf,
        states!,
        controls!,
        minDF,
        maxDF,
        savePlantData!,
        saveData,
        linearSpline,
        saveOptData,
        integrate!

include("setup.jl");
export  define,
        defineSolver!,
        configure!,
        dynamics!,
        constraints!,
        NLExpr

include("ps.jl");
include("diffeq.jl")

include("PrettyPlots/PrettyPlots.jl")
using .PrettyPlots

export  minDF,
        maxDF,
        plotSettings,
        _pretty_defaults,
        currentSettings,
        statePlot,
        controlPlot,
        costatesPlot,
        costatesPlots,
        allPlots,
        adjust_axis,
        mpcPlot,
        tPlot,
        optPlot,
        xlims!,
        ylims!,
        plot

# JuMP functions
export @NLexpression,
       @NLobjective,
       @NLparameter,
       @NLconstraint,
       setvalue,
       getvalue,
       setRHS

end # module
