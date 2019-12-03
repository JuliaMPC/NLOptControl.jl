isdefined(Base, :__precompile__) && __precompile__()

module NLOptControl
#TODO  enable setvalue() functionality

using JuMP
import JuMP.set_normalized_rhs, JuMP.getvalue, JuMP.setvalue, JuMP.@NLexpression, JuMP.@NLobjective, JuMP.@NLparameter, JuMP.@NLconstraint, JuMP.backend
using Ipopt
# using KNITRO
using FastGaussQuadrature
using DataFrames # https://discourse.julialang.org/t/dataframes-0-11-released/7296
using CSV
using Interpolations

include("Base.jl")
using .Base

# include("MPC_Module.jl")
# using .MPC_Module

################################################################################
# Model Classfield
################################################################################
struct OCP
  # general properties
  state::State              # state data
  control::Control          # control data
  tf::Any                   # final time
  t0::Any                   # initial time TODO consider getting ride of this! or replace it with n.mpc.v.t0Param
  tV::Any                   # vector for use with time varying constraints

  # boundary constraits
  X0::Array{Float64,1}      # initial state conditions
  X0_tol::Array{Float64,1}  # initial state tolerance
  x0s::Array{JuMP.Variable,1}
  XF::Array{Float64,1}      # final state conditions
  XF_tol::Array{Float64,1}  # final state tolerance
  xFs::Array{JuMP.Variable,1}

  # constant bounds on state variables
  XL::Array{Float64,1}
  XU::Array{Float64,1}

  # variables for linear bounds on state variables
  mXL::Array{Any,1}           # slope on XL -> time always starts at zero
  mXU::Array{Any,1}           # slope on XU -> time always starts at zero
  XL_var::Any                 # time varying lower bounds on states NOTE Not used currently
  XU_var::Any                 # time varying upper bounds on states NOTE Not used currently

  # constant bounds on control variables
  CL::Array{Float64,1}
  CU::Array{Float64,1}

  # ps method data
  Nck::Array{Int64,1}             # number of collocation points per interval
  Nck_cum::Array{Int64,1}         # cumulative number of points per interval
  Nck_full::Array{Int64,1}        # [0;cumsum(n.ocp.Nck+1)]
  Ni::Int64                       # number of intervals
  tau::Array{Array{Float64,1},1}  # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
  ts::Array{Array{Float64,1},1}   # time scaled based off of tau
  w::Array{Array{Float64,1},1}    # weights
  ws::Array{Array{Any,1},1}       # scaled weights
  DMatrix::Array{Array{Any,2},1}  # differention matrix
  IMatrix::Array{Array{Any,2},1}  # integration matrix

  # tm method data
  N::Int64                      # number of points in discretization
  dt::Array{Any,1}              # array of dts

  mdl::JuMP.Model               # JuMP model
  params                        # parameters for the models
  DXexpr
  NLcon# NOTE Not used currently

  # scaling factors
  XS::Array{Float64,1}           # scaling factors on states
  CS::Array{Float64,1}           # scaling factors on controls
end

# Default Constructor
function OCP()
OCP(
      State(),            # state data
      Control(),          # control data
      Any,                # final time
      0.0,                # initial time
      Any,                # optional vector for use with time varying constraints
      Float64[],          # initial state conditions
      Float64[],          # tolerances on inital state constraint
      Vector{Any}[],
      Float64[],          # final state conditions
      Vector{Any}[],
      Float64[],          # tolerances on final state constraint
      Float64[],          # XL
      Float64[],          # XU
      Any[],              # slopes on XL -> time always starts at zero
      Any[],              # slopes on XU -> time always starts at zero
      Any,                # time varying lower bounds on states
      Any,                # time varying upper bounds on states
      Float64[],          # CL
      Float64[],          # CU
      Int[],              # number of collocation points per interval
      Int[],              # Nck_cum
      Int[],              # Nck_full
      0,                  # number of intervals
      Vector{Float64}[],  # tau
      Vector{Any}[],      # ts
      Vector{Float64}[],  # weights
      Vector{Any}[],      # scaled weights
      Matrix{Any}[],      # DMatrix
      Matrix{Any}[],      # IMatrix
      0,                  # number of points in discretization
      Any[],              # array of dts
      JuMP.Model(),       # JuMP model
      Any[],
      Any[],
      Any[],
      Float64[],
      Float64[]
      )
end

struct OCPFlags
  defined::Bool  # a bool to tell if define() has been called
end

function OCPFlags()
  OCPFlags(
   false
  )
end

struct MPCFlags
 defined::Bool
 goalReached::Bool
 simFailed::Array{Any,1}   # a bool to indicate that the simulation failed and a symbol to indicate why
 ipDefined::Bool
 epDefined::Bool
end

function MPCFlags()
 MPCFlags(
  false,
  false,
  [false, NaN],
  false,
  false
  )
end

struct Flags
  ocp::OCPFlags
  mpc::MPCFlags
end

function Flags()
  Flags(
   OCPFlags(),
   MPCFlags()
  )
end

abstract type AbstractNLOpt end
struct NLOpt <: AbstractNLOpt
  # major data structs
  ocp::OCP
  mpc::MPC
  s::Settings
  r::Results
  f::Flags
end

# Default Constructor
function NLOpt()
NLOpt(
  OCP(),
  MPC(),
  Settings(),
  Results(),
  Flags()
    )
end

# scripts
include("utils.jl")
include("setup.jl");
include("ps.jl");
include("diffeq.jl")

include("PrettyPlots/PrettyPlots.jl")
using .PrettyPlots

export
       # Base functions  TODO make a hamiltonian function
       evalConstraints!,
       postProcess!,
       optimize!,
       interpolateLagrange!,
       interpolateLinear!,
       interpolate_lagrange,
       opt2dfs!,

       # setup functions
       define,
       configure!,
       dynamics!,
       constraints!,
       NLExpr,

       # extra functions
       defineSolver!,

       # math functions
       integrate!,

       # optimization related functions - utils.jl
       defineTolerances!,
       linearStateTolerances!,

       # data processing  - utils.jl
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

       # MPC_Module.jl
       MPC,
       defineMPC!,
       initOpt!,
       defineIP!,
       mapNames!,
       simIPlant!,
       updateX0!,
       currentIPState,
       goalReached!,
       simMPC!,
       plant2dfs!,
       predictX0!,

       # Objects
       NLOpt,

       # results
       resultsDir!,  # function to make a results folder

       #JuMP macros and functions
       @NLexpression,
       @NLobjective,
       @NLparameter,
       @NLconstraint,
       setvalue,
       getvalue,
       set_normalized_rhs,


       # PrettyPlots
        minDF,
        maxDF,
        plotSettings,
        _pretty_defaults,
        currentSettings,

        # NLOptControl plots
        statePlot,
        controlPlot,
        costatesPlot,
        costatesPlots,
        allPlots,
        adjust_axis,

        # MPC plots
        mpcPlot,
        tPlot,
        optPlot,

        # Plots.jl exported functions
        xlims!,
        ylims!,
        plot
end # module
