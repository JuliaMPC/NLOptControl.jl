isdefined(Base, :__precompile__) && __precompile__()

module NLOptControl

using JuMP
import JuMP.setRHS, JuMP.getvalue, JuMP.setvalue, JuMP.@NLexpression, JuMP.@NLobjective, JuMP.@NLparameter, JuMP.@NLconstraint, JuMP.internalmodel
using Ipopt
using KNITRO
using FastGaussQuadrature
using DataFrames
using Ranges

include("Base.jl")
using .Base

include("MPC_Module.jl")
using .MPC_Module

################################################################################
# Constants
################################################################################

const _Ipopt_defaults=Dict(
   :print_level                =>0,
   :warm_start_init_point      =>"yes",
   :tol                        =>1e-8,
   :max_iter                   =>3000,
   :max_cpu_time               =>1e6,
   :dual_inf_tol               =>1.,
   :constr_viol_tol            =>0.0001,
   :compl_inf_tol              =>0.0001,
   :acceptable_tol             =>1e-6,
   :acceptable_constr_viol_tol =>0.01,
   :acceptable_dual_inf_tol    =>1e-10,
   :acceptable_compl_inf_tol   =>0.01,
   :acceptable_obj_change_tol  =>1e20,
   :diverging_iterates_tol     =>1e20
)

const _Ipopt_MPC=Dict(
   :print_level                =>0,
   :warm_start_init_point      =>"yes",
   :tol                        =>5e-1,
   :max_iter                   =>500,
   :max_cpu_time               =>0.47,
   :dual_inf_tol               =>5.,
   :constr_viol_tol            =>1e-1,
   :compl_inf_tol              =>1e-1,
   :acceptable_tol             =>1e-2,
   :acceptable_constr_viol_tol =>0.01,
   :acceptable_dual_inf_tol    =>1e10,
   :acceptable_compl_inf_tol   =>0.01,
   :acceptable_obj_change_tol  =>1e20,
   :diverging_iterates_tol     =>1e20
)

const _KNITRO_defaults=Dict(
  :outlev                       =>1,
  :maxit                        =>0,
  :maxtime_real                 =>1.0e8,
  :infeastol                    =>1.0e-8,
  :feastol                      =>1.0e-6,
  :feastol_abs                  =>1e-3,
  :opttol                       =>1.0e-6,
  :opttol_abs                   =>1.0e-3,
  :algorithm                    =>0,
  :bar_initpt                   =>0,
  :bar_murule                   =>0,
  :bar_penaltycons              =>0,
  :bar_penaltyrule              =>0,
  :bar_switchrule               =>0,
  :linesearch                   =>0,
  :linsolver                    =>0
)

const _KNITRO_MPC=Dict(
  :outlev                       =>0,
  :maxit                        =>500,
  :maxtime_real                 =>0.47,
  :infeastol                    =>1e-2,
  :feastol                      =>1.0e20,
  :feastol_abs                  =>7e-2,
  :opttol                       =>1.0e20,
  :opttol_abs                   =>5e-1,
  :algorithm                    =>1,
  :bar_initpt                   =>3,
  :bar_murule                   =>4,
  :bar_penaltycons              =>1,
  :bar_penaltyrule              =>2,
  :bar_switchrule               =>2,
  :linesearch                   =>1,
  :linsolver                    =>2
)

################################################################################
# Basic Types
################################################################################
############################### control ########################################
type Control
  name::Vector{Any}
  description::Vector{Any}
end
function Control()
  Control([],
          []);
end

# ############################# state  ##########################################
type State
  name::Vector{Any}
  description::Vector{Any}
end
function State()
  State([],
        []);
end

############################## constraint ######################################
type Constraint
  name::Vector{Any}
  handle::Vector{Any}
  value::Vector{Any}
  nums  # range of indecies in g(x)
end
function Constraint()
  Constraint([],
             [],
             [],
             []);
end

############################### solver  ########################################
type Solver
    name
    settings
end

function Solver()
       Solver(:Ipopt,
              _Ipopt_defaults);
end

# Result Class
type Result
  t_ctr                       # time vector for control
  t_st                        # time vector for state
  x                           # JuMP states
  u                           # JuMP controls
  X                           # states
  U                           # controls
  CS                          # costates
  t_polyPts                   # time sample points for polynomials
  X_polyPts                   # state evaluated using Lagrange polynomial
  CS_polyPts                  # costate evaluated using Lagrange polynomial
  U_polyPts                   # control evaluated using Lagrane polynomial
  t_pts                       # vector time sample points for polynomials
  X_pts                       # vector state evaluated using Lagrange polynomial
  U_pts                       # vector control evaluated using Lagrane polynomial
  CS_pts                      # vector costate evaluated using Lagrange polynomial
  x0_con                      # handle for intial state constraints
  xf_con                      # handle for final state constraints
  dyn_con                     # dynamics constraints
  constraint::Constraint      # constraint handles and data
  eval_num                    # number of times optimization has been run
  iter_nums                   # mics. data, perhaps an iteration number for a higher level algorithm
  status                      # optimization status
  t_solve                     # solve time for optimization
  obj_val                     # objective function value
  dfs                         # results in DataFrame for plotting
  dfs_opt                     # optimization information in DataFrame for plotting
  dfs_plant                   # plant data
  dfs_con                     # constraint data
  results_dir                 # string that defines results folder
  main_dir                    # string that defines main folder
end

# Default Constructor
function Result()
Result( Vector{Any}[], # time vector for control
        Vector{Any}[], # time vector for state
        Matrix{Any}[], # JuMP states
        Matrix{Any}[], # JuMP controls
        Matrix{Any}[], # states
        Matrix{Any}[], # controls
        [],            # costates
        [],            # time sample points for polynomials
        [],            # state evaluated using Lagrange polynomial
        [],            # costate evaluated using Lagrange polynomial
        [],            # control evaluated using Lagrane polynomial
        Vector{Any}[], # vector time sample points for polynomials
        Vector{Any}[], # vector state evaluated using Lagrange polynomial
        Vector{Any}[], # vector control evaluated using Lagrane polynomial
        Vector{Any}[], # vector costate evaluated using Lagrange polynomial
        nothing,       # handle for intial state constraints
        nothing,       # handle for final state constraints
        nothing,       # dynamics constraint
        Constraint(),  # constraint data
        1,             # current evaluation number
        [],            # mics. data, perhaps an iteration number for a higher level algorithm
        Symbol,        # optimization status
        Float64,       # solve time for optimization
        Float64,       # objective function value
        [],            # results in DataFrame for plotting
        [],            # optimization information in DataFrame for plotting
        [],            # plant data
        [],            # constraint data
        string(pwd(),"/results/"),  # string that defines results folder
        pwd()                       # string that defines main folder
      );
end

# Settings Class
type Settings   # options
  solver::Solver                # solver information
  finalTimeDV::Bool
  integrationMethod::Symbol
  integrationScheme::Symbol
  MPC::Bool                     # bool for doing MPC
  save::Bool                    # bool for saving data
  reset::Bool                   # bool for reseting data
  evalConstraints::Bool         # bool for evaluating duals of the constraints
  evalCostates::Bool            # bool for evaluating costates
  tf_max::Any                   # maximum final time
  numInterpPts::Int64           # number of points to sample polynomial running through collocation points
end

# Default Constructor NOTE currently not using these, they get overwritten
function Settings()
        Settings(
         Solver(),           # default solver
         false,              # finalTimeDV
         :ts,                # integrationMethod
         :bkwEuler,          # integrationScheme
         false,              # bool for doing MPC
         true,               # bool for saving data
         false,              # bool for reseting data
         false,              # bool for evaluating duals of the constraints
         false,              # bool for evaluating costates
         400.0,              # maximum final time
         250                 # number of points to sample polynomial running through collocation points
                );
end

################################################################################
# Model Class
################################################################################
abstract type AbstractNLOpt end
type NLOpt <: AbstractNLOpt
  # general properties
  numStates::Int64          # number of states
  state::State              # state data
  numControls::Int64        # number of controls
  control::Control          # control data
  numStatePoints::Int64     # number of dvs per state
  numControlPoints::Int64   # numer of dvs per control
  lengthDV::Int64           # total number of dv discretizations per variables
  tf::Any                   # final time
  t0::Any                   # initial time TODO consider getting ride of this! or replace it with n.mpc.t0_param
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
  XL_var::Any                 # time varying lower bounds on states
  XU_var::Any                 # time varying upper bounds on states

  # constant bounds on control variables
  CL::Array{Float64,1}
  CU::Array{Float64,1}

  # ps method data
  Nck::Array{Int64,1}             # number of collocation points per interval
  Nck_cum::Array{Int64,1}         # cumulative number of points per interval
  Nck_full::Array{Int64,1}        # [0;cumsum(n.Nck+1)]
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

  # major data types
  mpc::MPC                      # mpc data
  mdl::JuMP.Model               # JuMP model
  s::Settings                   # settings
  r::Result                     # results
  params                        # parameters for the models
  DXexpr
  NLcon

  # problem state
  define::Bool                  # a bool to tell if define!() has been called
end

# Default Constructor
function NLOpt()
NLOpt(
      0,                  # number of states
      State(),            # state data
      0,                  # number of controls
      Control(),          # control data
      0,                  # number of dvs per state
      0,                  # number of dvs per control
      0,                  # total number of dv discretizations per variables
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
      MPC(),              # mpc data
      JuMP.Model(),       # JuMP model
      Settings(),
      Result(),
      Any[],
      Any[],
      Any[],
      false);
end

# scripts
include("utils.jl");
include("setup.jl");
include("ps.jl");
include("diffeq.jl")

export
       # Base functions  TODO make a costate function and hamiltonian function
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

       # MPC_Module.jl
       autonomousControl!,
       initializeMPC!,
       driveStraight!,
       updateX0!,
       simPlant!,
       simModel,
       MPC,

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
       setRHS
end # module
