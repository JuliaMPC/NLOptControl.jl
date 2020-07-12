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
export @NLexpression,
       @NLobjective,
       @NLparameter,
       @NLconstraint,
       setvalue,
       getvalue,
       setRHS

using Ipopt
using FastGaussQuadrature
using DataFrames
using CSV
using Interpolations
using Parameters
import Printf
import LinearAlgebra

include("math.jl")

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
@with_kw mutable struct OCP{T <: Number}
    # general properties
    state::State                        = State()                               # state data
    control::Control                    = Control()                             # control data
    tf::T                               = convert(T, 0)                         # final time
    t0::JuMP.JuMPTypes                  = @NLparameter(JuMP.Model(), x == 0)    # initial time # TODO: consider getting rid of this or replacing it with `n.mpc.v.t0Param`
    tV::Vector{T}                       = Vector{T}()                           # vector for use with time varying constraints

    # Boundary conditions
    X0::Vector{T}                       = Vector{T}()                           # initial state conditions
    X0_tol::Vector{T}                   = Vector{T}()                           # initial state tolerance
    x0s::Vector{JuMP.JuMPTypes}         = Vector{JuMP.JuMPTypes}()              # initial state variables
    XF::Vector{T}                       = Vector{T}()                           # final state conditions
    XF_tol::Vector{T}                   = Vector{T}()                           # final state tolerance
    xFs::Vector{JuMP.JuMPTypes}         = Vector{JuMP.JuMPTypes}()              # final state variables

    # Constant bounds on state variables
    XL::Vector{T}                       = Vector{T}()                           # Constant lower bound on state variables
    XU::Vector{T}                       = Vector{T}()                           # Constant upper bound on state variables

    # Variables for linear bounds on state variables
    mXL::Vector{T}                      = Vector{T}()                           # slope on XL -> time always starts at zero
    mXU::Vector{T}                      = Vector{T}()                           # slope on XU -> time always starts at zero
    XL_var::Matrix{T}                   = Matrix{T}(undef,0,0)                  # time varying lower bounds on states # ! NOTE: not used currently - was JuMP.Variable
    XU_var::Matrix{T}                   = Matrix{T}(undef,0,0)                  # time varying upper bounds on states # ! NOTE: not used currently - was JuMP.Variable

    # Constant bounds on control variables
    CL::Vector{T}                       = Vector{T}()                           # Constant lower bound on control variables
    CU::Vector{T}                       = Vector{T}()                           # Constant upper bound on control variables

    # Pseudospectral method data
    Nck::Vector{Int}                    = Vector{Int}()                         # number of collocation points per interval
    Nck_cum::Vector{Int}                = Vector{Int}()                         # cumulative number of points per interval
    Nck_full::Vector{Int}               = Vector{Int}()                         # [0;cumsum(n.ocp.Nck+1)]
    Ni::Int                             = Int(0)                                # number of intervals
    tau::Matrix{T}                      = Matrix{T}(undef,0,0)                  # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
    ts::Matrix{T}                       = Matrix{T}(undef,0,0)                  # time scaled based off of tau
    w::Matrix{T}                        = Matrix{T}(undef,0,0)                  # weights
    ws::Matrix{T}                       = Matrix{T}(undef,0,0)                  # scaled weights
    DMatrix::Vector{Matrix{T}}          = Vector{Matrix{T}}()                   # differention matrix
    IMatrix::Vector{Matrix{T}}          = Vector{Matrix{T}}()                   # integration matrix

    # tm method data
    N::Int                              = 0                                     # number of points in discretization
    dt::Vector{T}                       = Vector{T}()                           # array of dts

    mdl::JuMP.Model                     = JuMP.Model()                          # JuMP model
    params                              = Any[]                                 # parameters for the models
    DXexpr                              = Any[]                                 # ? DX expression
    NLcon                               = Any[]                                 # ! NOTE: not used currently

    # Scaling factors
    XS::Vector{T}                       = Vector{T}()                           # scaling factors on states
    CS::Vector{T}                       = Vector{T}()                           # scaling factors on controls
end

@with_kw mutable struct OCPFlags
  defined::Bool = false # a bool to tell if define() has been called
end

@with_kw mutable struct MPCFlags
    defined::Bool                           = false         #
    goalReached::Bool                       = false         #
    simFailed::Vector{Union{Bool, Symbol}}  = [false, :NaN] # a bool to indicate that the simulation failed and a symbol to indicate why
    ipDefined::Bool                         = false         #
    epDefined::Bool                         = false         #
end

@with_kw mutable struct Flags
    ocp::OCPFlags   = OCPFlags()    #
    mpc::MPCFlags   = MPCFlags()    #
end

@with_kw mutable struct NLOpt{T <: Number}
    # major data structs
    ocp::OCP{T}     = OCP{T}()      #
    mpc::MPC{T}     = MPC{T}()      #
    s::Settings     = Settings{T}() #
    r::Results{T}   = Results{T}()  #
    f::Flags        = Flags()       #
end
NLOpt() = NLOpt{Float64}()

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

end # module
