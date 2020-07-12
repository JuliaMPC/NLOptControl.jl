module NLOptBase

using JuMP
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
        simulationModes

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
    model::JuMP.Model                    = JuMP.Model()              #
end

# Constraint
@with_kw mutable struct Constraint{ T <: Number }
    name::Vector{Symbol}                    = Vector{Symbol}()                  #
    handle::Vector{JuMP.ConstraintRef}      = Vector{JuMP.ConstraintRef}()      #
    value::Vector{T}                        = Vector{T}()                       #
    nums::Vector{Any}                       = Vector{Any}()                     # range of indices in g(x)
end

# Solver
@with_kw mutable struct Solver
    name::Symbol                = :Ipopt            #
    settings::Dict{Symbol,Any}  = _Ipopt_defaults   #
end

# Plant Results
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

# Optimal Control Problem (OCP) Results
@with_kw mutable struct OCPResults{ T <: Number }
    tctr::Vector{Any}                   = Vector{T}()                           # Time vector for control
    tst::Vector{Any}                    = Vector{T}()                           # Time vector for state
    x::Matrix{JuMP.JuMPTypes}           = Matrix{JuMP.JuMPTypes}(undef,0,0)     # JuMP states
    u::Matrix{JuMP.JuMPTypes}           = Matrix{JuMP.JuMPTypes}(undef,0,0)     # JuMP controls
    X::Matrix{T}                        = Matrix{T}(undef,0,0)                  # States
    U::Matrix{T}                        = Matrix{T}(undef,0,0)                  # Controls
    X0::Vector{T}                       = Vector{T}()                           # Initial states for OCP
    CS::Vector{T}                       = Vector{T}()                           # Costates
    tpolyPts::Vector{T}                 = Vector{T}()                           # Time sample points for polynomials  (NOTE these interpolated solutions were developed for calculating error, between them and a known Optimal solution)
    XpolyPts::Matrix{T}                 = Matrix{T}(undef,0,0)                  # State evaluated using Lagrange/Linear polynomial
    CSpolyPts::Matrix{T}                = Matrix{T}(undef,0,0)                  # Costate evaluated using Lagrange/Linear polynomial
    UpolyPts::Matrix{T}                 = Matrix{T}(undef,0,0)                  # Control evaluated using Lagrane/Linear polynomial
    AlltpolyPts::Vector{T}              = Vector{T}()                           # Time sample points for polynomials
    AllXpolyPts::Matrix{T}              = Matrix{T}(undef,0,0)                  # State evaluated using Lagrange/Linear polynomial
    AllCSpolyPts::Matrix{T}             = Matrix{T}(undef,0,0)                  # Costate evaluated using Lagrange/Linear polynomial
    AllUpolyPts::Matrix{T}              = Matrix{T}(undef,0,0)                  # Control evaluated using Lagrane/Linear polynomial
    tpts::Vector{T}                     = Vector{T}()                           # Vector time sample points
    Xpts::Matrix{T}                     = Matrix{T}(undef,0,0)                  # Vector state sample points
    Upts::Matrix{T}                     = Matrix{T}(undef,0,0)                  # Vector control sample points
    CSpts::Matrix{T}                    = Matrix{T}(undef,0,0)                  # Vector costate sample points
    x0Con::Array{JuMP.ConstraintRef}    = Array{JuMP.ConstraintRef}(undef)      # Handle for initial state constraints
    x0sCon::Array{JuMP.ConstraintRef}   = Array{JuMP.ConstraintRef}(undef)      # ? Unsure what this is yet (slack variable constraints?)
    xfCon::Array{JuMP.ConstraintRef}    = Array{JuMP.ConstraintRef}(undef)      # Handle for final state constraints
    xfsCon::Array{JuMP.ConstraintRef}   = Array{JuMP.ConstraintRef}(undef)      # ? Unsure what this is yet (slack variable constraints?)
    dynCon::Array{JuMP.ConstraintRef}   = Array{JuMP.ConstraintRef}(undef)      # Dynamics constraints
    constraint::Constraint{T}           = Constraint{T}()                       # Constraint handles and data
    evalNum::Int                        = 1                                     # Number of times optimization has been run
    iterNum                             = []                                    # Mics. data, perhaps an iteration number for a higher level algorithm
    status::Symbol                      = :undef                                # Optimization status
    tSolve::T                           = convert(T,0)                          # Solve time for optimization
    objVal::T                           = convert(T,0)                          # Objective function value
    dfs::Vector{DataFrame}              = Vector{DataFrame}()                   # Results in DataFrame for plotting
    dfsOpt::DataFrame                   = DataFrame()                           # Optimization information in DataFrame for plotting
    dfsCon::DataFrame                   = DataFrame()                           # Constraint data
end

# Combined Results
@with_kw mutable struct Results{T <: Number}
    ocp::OCPResults{T}          = OCPResults{T}()           #
    ip::PlantResults{T}         = PlantResults{T}()         #
    ep::PlantResults{T}         = PlantResults{T}()         #
    resultsDir::AbstractString  = joinpath(pwd(),"results") # string that defines results folder
    mainDir::AbstractString     = pwd()                     # string that defines main folder
end

# Model-Predictive Contrl (MPC) Settings Class
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

# Optimal Control Problem (OCP) Settings Class
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

# Combined Settings
@with_kw mutable struct Settings{ T <: Number }
    ocp::OCPSettings{T} = OCPSettings{T}()  #
    mpc::MPCSettings    = MPCSettings()     #
end

end # module
