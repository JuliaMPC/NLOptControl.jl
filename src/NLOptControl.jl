module NLOptControl

using Media, Dierckx, Parameters, Interpolations, FastGaussQuadrature, Polynomials, JuMP

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
  numControls::Int64        # number of controls
  numPoints::Array{Int64,1} # number of dv discretization within each interval
  lengthDV::Int64           # total number of dv discretizations per variables
  tf::Any                   # final time

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
      0,                  # number of controls
      Int[],              # number of dv discretization within each interval
      0,                  # total number of dv discretizations per variables
      Any,                # final time
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
      0,                  # number of points in discretization
      Any[],              # array of dts
      false,              # finalTimeDV
      :ts,                # integrtionMethod
      :bkwEuler           # integrationScheme
    );
end

# scripts
include("utils.jl");
include("setup.jl")
include("LGR.jl");
include("ocp.jl")

       # Functions
export
       # setup functions
       NLOpt, define, configure,
       OCPdef,

       # ps functions
       LGR, lgr_diff, LGR_matrices,
       scale_tau, scale_w, create_intervals,
       lagrange_basis_poly, interpolate_lagrange,
       polyDiff,
       D_matrix,
       lepoly, poldif,

       # math functions
       integrate,

       # Objects
       NLOpt
end
