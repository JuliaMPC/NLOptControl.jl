using JuMP
using OrdinaryDiffEq
using DiffEqBase
using DataFrames
using Interpolations
using CSV
using Parameters
using Printf

# These functions are required for NLOptMPC.jl and PrettyPlots.jl (resultsDir!)
export  State,
        Control,
        Constraint,
        Results,
        Settings,
        _Ipopt_defaults,
        _Ipopt_MPC,
        simulationModes,
        MPC

# IPOPT Settings
# :print_level      : Print level
# :constr_viol_tol  : Absolute tolerance on the constraint violation.
# :max_iter         : Maximum number of iterations.
# :max_cpu_time     : A limit on CPU seconds that Ipopt can use to solve one problem
const _Ipopt_defaults=Dict{Symbol,Any}(
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

# IPOPT Settings for Model Predictive Control
const _Ipopt_MPC = Dict{Symbol,Any}(
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

# Simulation Modes
const simulationModes = [ :OCP , :IP , :IPEP , :EP]

################################
# Optimal Control Common Types #
################################

# Control
@with_kw mutable struct Control
    name::Vector{Symbol}                = Vector{Symbol}[]          #
    description::Vector{AbstractString} = Vector{AbstractString}[]  #
    num::Int                            = Int(0)                    #
    pts::Int                            = Int(0)                    #
end

# State
@with_kw mutable struct State
    name::Vector{Symbol}                 = Vector{Symbol}()          #
    description::Vector{AbstractString}  = Vector{AbstractString}()  #
    num::Int                             = 0                         #
    pts::Int                             = 0                         #
    model::Any                    = Any           #
end

# Constraint
@with_kw mutable struct Constraint{ T <: Number }
    name::Vector{Symbol}                    = Vector{Symbol}()                  #
    handle::Vector{Any}      = Vector{Any}()      #
    value::Vector{T}                        = Vector{T}()                       #
    nums::Vector{Any}                       = Vector{Any}()                     # range of indices in g(x)
end

# Solver
@with_kw mutable struct Solver
    name::Symbol                = :Ipopt            #
    settings::Dict{Symbol,Any}  = _Ipopt_defaults   #
end

# Optimal Control Problem (OCP) Results
@with_kw mutable struct OCPResults{ T <: Number }
    tctr::Vector{Any}                           = Vector{T}()                           # Time vector for control
    tst::Vector{Any}                            = Vector{T}()                           # Time vector for state
    x           = Matrix{Any}[]      # JuMP states
    u           = Matrix{Any}[]      # JuMP controls
    X                        = Matrix{T}[]                   # States
    U                        =  Matrix{T}[]                 # Controls
    X0                              = Vector{T}[]                           # Initial states for OCP
    CS                               = []                           # Costates
    tpolyPts                         = []                           # Time sample points for polynomials  (NOTE these interpolated solutions were developed for calculating error, between them and a known Optimal solution)
    XpolyPts                = []                   # State evaluated using Lagrange/Linear polynomial
    CSpolyPts                = []               # Costate evaluated using Lagrange/Linear polynomial
    UpolyPts                 = []                  # Control evaluated using Lagrane/Linear polynomial
    AlltpolyPts                      = []                           # Time sample points for polynomials
    AllXpolyPts              = []                  # State evaluated using Lagrange/Linear polynomial
    AllCSpolyPts             = []                   # Costate evaluated using Lagrange/Linear polynomial
    AllUpolyPts              = []                   # Control evaluated using Lagrane/Linear polynomial
    tpts::Vector{T}                             = Vector{T}()                           # Vector time sample points
    Xpts                    = Vector{T}[]                   # Vector state sample points
    Upts                     = Vector{T}[]                  # Vector control sample points
    CSpts                    = Vector{T}[]                  # Vector costate sample points
    x0Con          = nothing          # Handle for initial state constraints
    x0sCon          = nothing           # ? Unsure what this is yet (slack variable constraints?)
    xfCon           =  nothing          # Handle for final state constraints
    xfsCon          = nothing     # ? Unsure what this is yet (slack variable constraints?)
    dynCon = nothing   # Dynamics constraints
    constraint::Constraint{T}                   = Constraint{T}()                       # Constraint handles and data
    evalNum::Int                                = 1                                     # Number of times optimization has been run
    iterNum                                     = []                                    # Mics. data, perhaps an iteration number for a higher level algorithm
    status::Symbol                              = :undef                                # Optimization status
    tSolve::T                                   = convert(T,0)                          # Solve time for optimization
    objVal::T                                   = convert(T,0)                          # Objective function value
    dfs::Vector{DataFrame}                      = Vector{DataFrame}()                   # Results in DataFrame for plotting
    dfsOpt::DataFrame                           = DataFrame()                           # Optimization information in DataFrame for plotting
    dfsCon::DataFrame                           = DataFrame()                           # Constraint data
end

# Optimal Control Problem (OCP) Settings
@with_kw mutable struct OCPSettings{ T <: Number }
    solver::Solver              = Solver()          # solver information
    finalTimeDV::Bool           = false             #
    integrationMethod::Symbol   = :ts               #
    integrationScheme::Symbol   = :bkwEuler         #
    save::Bool                  = true              # bool for saving data
    reset::Bool                 = false             # bool for reseting data
    evalConstraints::Bool       = false             # bool for evaluating duals of the constraints
    evalCostates::Bool          = false             # bool for evaluating costates
    tfMin::T                    = convert(T,   0.0) # minimum final time
    tfMax::T                    = convert(T, 400.0) # maximum final time # ? Should choose something else as default?
    tfOptimal::Union{Bool, T}   = false             # known optimal final time
    numInterpPts::Int           = 250               # number of points to sample polynomial running through collocation points
    cacheOnly::Bool             = false             # bool for only caching the results when using optimize!()
    linearInterpolation::Bool   = false             # bool for using linear interpolation even if integrationMethod ==:ps
    interpolationOn::Bool       = false             # bool to indicate if user wants solution interpolated for them
    x0slackVariables::Bool      = false             #
    xFslackVariables::Bool      = false             #
end

# Optimal Control Problem (OCP) Flags
@with_kw mutable struct OCPFlags
  defined::Bool = false # a bool to tell if define() has been called
end

# Optimal Control Problem (OCP)
@with_kw mutable struct OCP{T <: Number}

    # General properties
    state::State                            = State()                               # state data
    control::Control                        = Control()                             # control data
    tf             = Any                         # final time
    t0::JuMP.JuMPTypes                      = @NLparameter(JuMP.Model(), x == 0)    # initial time # TODO: consider getting rid of this or replacing it with `n.mpc.v.t0Param`
    tV                           = Any                          # vector for use with time varying constraints

    # Boundary conditions
    X0::Vector{T}                           = Vector{T}()                           # initial state conditions
    X0_tol::Vector{T}                       = Vector{T}()                           # initial state tolerance
    x0s::Vector{JuMP.JuMPTypes}             = Vector{JuMP.JuMPTypes}()              # initial state variables
    XF::Vector{T}                           = Vector{T}()                           # final state conditions
    XF_tol::Vector{T}                       = Vector{T}()                           # final state tolerance
    xFs::Vector{JuMP.JuMPTypes}             = Vector{JuMP.JuMPTypes}()              # final state variables

    # Constant bounds on state variables
    XL::Vector{T}                           = Vector{T}()                           # Constant lower bound on state variables
    XU::Vector{T}                           = Vector{T}()                           # Constant upper bound on state variables

    # Variables for linear bounds on state variables
    mXL::Vector{Bool}                       = Vector{Bool}()                        # slope on XL -> time always starts at zero
    mXU::Vector{Bool}                       = Vector{Bool}()                        # slope on XU -> time always starts at zero
    XL_var = Any[]     # time varying lower bounds on states # ! NOTE: not used currently - was JuMP.Variable
    XU_var = Any[]     # time varying upper bounds on states # ! NOTE: not used currently - was JuMP.Variable

    # Constant bounds on control variables
    CL::Vector{T}                           = Vector{T}()                           # Constant lower bound on control variables
    CU::Vector{T}                           = Vector{T}()                           # Constant upper bound on control variables

    # Pseudospectral method data
    Nck::Vector{Int}                        = Vector{Int}()                         # number of collocation points per interval
    Nck_cum::Vector{Int}                    = Vector{Int}()                         # cumulative number of points per interval
    Nck_full::Vector{Int}                   = Vector{Int}()                         # [0;cumsum(n.ocp.Nck+1)]
    Ni::Int                                 = Int(0)                                # number of intervals
    tau::Array{Array{T,1},1}                = Array{Array{T,1},1}()                 # Node points ---> Nc increasing and distinct numbers ∈ [-1,1]
    ts::Array{Array{T,1},1} = Array{Array{T,1},1}()                                 # time scaled based off of tau
    w::Array{Array{T,1},1}  = Array{Array{T,1},1}()                                 # weights
    ws::Array{Array{T,1},1} = Array{Array{T,1},1}()                                 # scaled weights
    DMatrix::Vector{Matrix{T}}              = Vector{Matrix{T}}()                   # differention matrix
    IMatrix::Vector{Matrix{T}}              = Vector{Matrix{T}}()                   # integration matrix

    # tm method data
    N::Int                              = 0                                     # number of points in discretization
    dt::Array{Any,1}                    = Array{Any,1}()       # array of dts

    mdl::JuMP.Model                     = JuMP.Model()                          # JuMP model
    params                              = Any[]                                 # parameters for the models
    DXexpr                              = Any[]                                 # ? DX expression
    NLcon                               = Any[]                                 # ! NOTE: not used currently

    # Scaling factors
    XS::Vector{T}                       = Vector{T}()                           # scaling factors on states
    CS::Vector{T}                       = Vector{T}()                           # scaling factors on controls
end

##################################
# Model-Predictive Control Types #
##################################

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
    tp::Any  = Any   # prediction time (if finalTimeDV == true -> this is not known before optimization)
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

# Model-Predictive Contrl (MPC) Settings
@with_kw mutable struct MPCSettings
    on::Bool                    = false     #
    mode::Symbol                = :OCP      #
    predictX0::Bool             = false     #
    fixedTp::Bool               = true      #
    IPKnown::Bool               = true      #
    saveMode::Symbol            = :all      #
    maxSim::Int                 = 100       # maximum number of total MPC updates
    shiftX0::Bool               = true      # a bool to indicate that a check on the feasibility of the linear constraints for the initial conditions should be performed
    lastOptimal::Bool           = true      # a bool to indicate that the previous solution should be used if the current solution in not optimal
    printLevel::Int             = 2         #
    expandGoal::Bool            = true      # bool to indicate if the goal needs to be expanded
    enlargeGoalTolFactor::Int   = 2         # scaling factor to enlare the goal
    onlyOptimal::Bool           = false     #
end

# Plant Results
# TODO: Figure out what these "Any" types are
@with_kw mutable struct PlantResults{ T <: Number }
  plant::DataFrame            = DataFrame()             # plant
  X0p::Vector{Any}            = Vector{Any}()           # X0p
  X0a::Vector{Any}            = Vector{Any}()           # X0a
  X0e::Vector{Any}            = Vector{Any}()           # X0e
  e::Matrix{Any}              = Matrix{Any}(undef,0,0)  # e
  dfsplant::Vector{DataFrame} = Vector{DataFrame}()     # dfsplant
  dfsplantPts::DataFrame      = DataFrame()             # dfsplantPts
  dfsX0p::DataFrame           = DataFrame()             # dfsX0p
  dfsX0a::DataFrame           = DataFrame()             # dfsX0a
  dfsX0e::DataFrame           = DataFrame()             # dfsX0e
  dfse::DataFrame             = DataFrame()             # dfse
end

# Model-Predictive Control (MPC) Flags
@with_kw mutable struct MPCFlags
    defined::Bool                           = false         #
    goalReached::Bool                       = false         #
    simFailed::Vector{Union{Bool, Symbol}}  = [false, :NaN] # a bool to indicate that the simulation failed and a symbol to indicate why
    ipDefined::Bool                         = false         #
    epDefined::Bool                         = false         #
end

# Model-Predictive Control (MPC)
@with_kw mutable struct MPC{ T <: Number }
    v::MPCvariables{T}  = MPCvariables{T}() #
    ip::IP              = IP()              #
    ep::EP              = EP()              #
end

##################
# Combined Types #
##################

# Combined Settings
@with_kw mutable struct Settings{ T <: Number }
    ocp::OCPSettings{T} = OCPSettings{T}()  #
    mpc::MPCSettings    = MPCSettings()     #
end

# Combined Results
@with_kw mutable struct Results{T <: Number}
    ocp::OCPResults{T}          = OCPResults{T}()           #
    ip::PlantResults{T}         = PlantResults{T}()         #
    ep::PlantResults{T}         = PlantResults{T}()         #
    resultsDir::AbstractString  = joinpath(pwd(),"results") # string that defines results folder
    mainDir::AbstractString     = pwd()                     # string that defines main folder
end

# Combined Flags
@with_kw mutable struct Flags
    ocp::OCPFlags   = OCPFlags()    #
    mpc::MPCFlags   = MPCFlags()    #
end

# Nonlinear Optimal Control (NLOpt) Problem
@with_kw mutable struct NLOpt{T <: Number}
    # major data structs
    ocp::OCP{T}     = OCP{T}()      #
    mpc::MPC{T}     = MPC{T}()      #
    s::Settings     = Settings{T}() #
    r::Results{T}   = Results{T}()  #
    f::Flags        = Flags()       #
end

# Defaualt Number subtype is Float64
NLOpt() = NLOpt{Float64}()
