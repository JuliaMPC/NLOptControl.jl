isdefined(Base, :__precompile__) && __precompile__()

module NLOptControl

using JuMP
import JuMP.setRHS, JuMP.getvalue, JuMP.setvalue, JuMP.@NLexpression, JuMP.@NLobjective, JuMP.@NLparameter, JuMP.@NLconstraint, JuMP.internalmodel
using Ipopt
using KNITRO
using FastGaussQuadrature
using DataFrames # https://discourse.julialang.org/t/dataframes-0-11-released/7296
using CSV
using Interpolations

include("Base.jl")
using .Base

include("MPC_Module.jl")
using .MPC_Module

################################################################################
# Model Classfield
################################################################################
type OCP
  # general properties
  state::State              # state data
  control::Control          # control data
  tf::Any                   # final time
  t0::Any                   # initial time TODO consider getting ride of this! or replace it with n.mpc.v.t0Param
  tV::Any                   # vector for use with time varying constraints

  # boundary constraits
  X0::Array{Float64,1}      # initial state conditions
  X0_tol::Array{Float64,1}  # initial state tolerance
  XF::Array{Float64,1}      # final state conditions
  XF_tol::Array{Float64,1}  # final state tolerance

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
      Float64[],          # final state conditions
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
      Any[]
      )
end

type OCPFlags
  defined::Bool  # a bool to tell if define() has been called
end

function OCPFlags()
  OCPFlags(
   false
  )
end

type MPCFlags
 defined::Bool
 goalReached::Bool
 ipDefined::Bool
 epDefined::Bool
end

function MPCFlags()
 MPCFlags(
  false,
  false,
  false,
  false
 )
end

type Flags
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
type NLOpt <: AbstractNLOpt
  # major data types
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
include("utils.jl");
include("setup.jl");
include("ps.jl");
include("diffeq.jl")

export
       # Base functions  TODO make a hamiltonian function
       evalConstraints!,
       postProcess!,
       optimize!,
       interpolateLagrange!,
       interpolateLinear!,
       interpolate_lagrange,

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

       # MPC_Module.jl
       MPC,
       defineMPC!,
       initOpt!,
       defineIP!,
       mapNames!,
       simIPlant!,
       updateX0!,
       currentIPState,
       updateOCPState!,
       goalReached!,
       simMPC!,

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
       setRHS,

       # TMP
       State,
       Control,
       initState,
       plant2dfs!,
       scale_tau,
       intervals,
       interpolate_lagrange

end # module
