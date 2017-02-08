module NLOptControl

using Media, Dierckx, Parameters, Interpolations, FastGaussQuadrature, Polynomials, JuMP, SymPy, VehicleModels
# To copy a particular piece of code (or function) in some location
macro def(name, definition)
  return quote
    macro $name()
      esc($(Expr(:quote,definition)))
    end
  end
end

# Model Class
abstract AbstractNLOpt
type NLOpt <: AbstractNLOpt
  # general properties
  stateEquations
  numStates::Int64          # number of states
  state
  numControls::Int64        # number of controls
  control
  numPoints::Array{Int64,1} # number of dv discretization within each interval
  numStatePoints::Int64     # number of dvs per state
  numControlPoints::Int64   # numer of dvs per control
  lengthDV::Int64           # total number of dv discretizations per variables
  tf::Any                   # final time
  t0::Any                   # initial time
  tf_max::Any               # maximum final time

  # boundary constraits
  X0::Array{Float64,1}      # initial state conditions
  XF::Array{Float64,1}      # final state conditions

  # linear bounds on variables
  XL::Array{Float64,1}
  XU::Array{Float64,1}
  CL::Array{Float64,1}
  CU::Array{Float64,1}

  # ps method data
  Nck::Array{Int64,1}           # number of collocation points per interval
  Ni::Int64                     # number of intervals
  τ::Array{Array{Float64,1},1}  # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
  ts::Array{Array{Any,1},1}     # time scaled based off of τ
  ω::Array{Array{Float64,1},1}  # weights
  ωₛ::Array{Array{Any,1},1}     # scaled weights
  DMatrix::Array{Array{Any,2},1}# differention matrix
  IMatrix::Array{Array{Any,2},1}# integration matrix

  # tm method data
  N::Int64                      # number of points in discretization
  dt::Array{Any,1}              # array of dts

  # options
  finalTimeDV::Bool
  integrationMethod::Symbol
  integrationScheme::Symbol
end

# Default Constructor
function NLOpt()
NLOpt(Any,                # state equations
      0,                  # number of states
      Symbol[],
      0,                  # number of controls
      Symbol[],
      Int[],              # number of dv discretization within each interval
      0,                  # number of dvs per state
      0,                  # number of dvs per control
      0,                  # total number of dv discretizations per variables
      Any,                # final time
      Any,                # initial time
      Any,                # maximum final time
      Float64[],          # initial state conditions
      Float64[],          # final state conditions
      Float64[],          # XL
      Float64[],          # XU
      Float64[],          # CL
      Float64[],          # CU
      Int[],              # number of collocation points per interval
      0,                  # number of intervals
      Vector{Float64}[],  # τ
      Vector{Any}[],      # ts
      Vector{Float64}[],  # weights
      Vector{Any}[],      # scaled weights
      Matrix{Any}[],      # DMatrix
      Matrix{Any}[],      # IMatrix
      0,                  # number of points in discretization
      Any[],              # array of dts
      false,              # finalTimeDV
      :ts,                # integrtionMethod
      :bkwEuler           # integrationScheme
    );
end

# Result Class
abstract AbstractNLOpt
type Result <: AbstractNLOpt
  t_ctr  # time vector for control
  t_st   # time vector for state
  x      # JuMP states
  u      # JuMP controls
  X      # states
  U      # controls
  x0_con # handle for intial state constraints
  xf_con # handle for final state constraints
  dyn_con # dynamics constraints
  constraint  # constraint handels and data
  eval_num::Int64 # number of times optimization has been run
  status::Vector{Symbol} # optimization status
  t_solve::Vector{Float64} # solve time for optimization
  obj_val::Vector{Float64} # objective function value
end

# Default Constructor
function Result()
Result( Vector{Any}[], # time vector for control
        Vector{Any}[], # time vector for state
        Matrix{Any}[], # JuMP states
        Matrix{Any}[], # JuMP controls
        Matrix{Any}[], # states
        Matrix{Any}[], # controls
        nothing,       # handle for intial state constraints
        nothing,       # handle for final state constraints
        nothing,       # dynamics constraint
        nothing,       # constraint data
        0, # number of times optimization has been run
        Symbol[], # optimization status
        Float64[], # solve time for optimization
        Float64[]  # objective function value
      );
end

# Settings Class
abstract AbstractNLOpt
type Settings <: AbstractNLOpt
  # plotting
  lw1::Float64 # line width 1
  lw2::Float64 # line width 2
end

# Default Constructor
function Settings()
Settings(8, # line width 1
         3  # line width 2
        );
end

# scripts
include("utils.jl");
include("setup.jl")
include("ps.jl");
include("ocp.jl")
include("post_processing.jl")

       # Functions
export
       # setup functions
       NLOpt, define, configure,
       OCPdef,

       # ps functions
       LGR_matrices,
       scale_tau, scale_w, createIntervals,
       lagrange_basis_poly, interpolate_lagrange,
       polyDiff,
       D_matrix,

       # math functions
       integrate,

       # optimization functions
       optimize,

       # data processing
       postProcess,
       newConstraintData,
       evalConstraints,
       stateNames,
       controlNames,

       # Objects
       NLOpt,
       Result,
       Settings
end
