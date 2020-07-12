module NLOptBase

using JuMP
using DataFrames
using Interpolations
using CSV
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
mutable struct Control
    name::Vector{Symbol}
    description::Vector{AbstractString}
    num::Int
    pts::Int
    Control() = new(
        Vector{Symbol}[],
        Vector{AbstractString}[],
        Int(0),
        Int(0)
    )
end

# State
mutable struct State
    name::Vector{Symbol}
    description::Vector{AbstractString}
    num::Int
    pts::Int
    model::JuMP.Model
    State() = new(
        Vector{Symbol}(),
        Vector{AbstractString}(),
        0,
        0,
        JuMP.Model()
    )
end

# Constraint
mutable struct Constraint{ T <: Number }
    name::Vector{Symbol}
    handle::Vector{JuMP.ConstraintRef}
    value::Vector{T}
    nums::Vector{Any} # range of indices in g(x)
    Constraint(T::DataType=Float64) = new{T}(
        Vector{Symbol}(),
        Vector{JuMP.ConstraintRef}(),
        Vector{T}(),
        Vector{Any}()
    )
end

# Solver
mutable struct Solver
    name::Symbol
    settings::Dict{Symbol,Any}
    Solver() = new(:Ipopt, _Ipopt_defaults)
end

# Plant Results
mutable struct PlantResults{ T <: Number }
  plant::DataFrame
  X0p::Vector{Any}
  X0a::Vector{Any} # closest one in time to X0e
  X0e::Vector{Any}
  e::Matrix{Any}
  dfsplant::Vector{DataFrame}
  dfsplantPts::DataFrame # plant data extracted into a single DataFrame
  dfsX0p::DataFrame
  dfsX0a::DataFrame
  dfsX0e::DataFrame
  dfse::DataFrame
  PlantResults(T::DataType=Float64) = new{T}(
      DataFrame(),             # plant
      Vector{Any}(),           # X0p
      Vector{Any}(),           # X0a
      Vector{Any}(),           # X0e
      Matrix{Any}(undef,0,0),  # e
      Vector{DataFrame}(),     # dfsplant
      DataFrame(),             # dfsplantPts
      DataFrame(),             # dfsX0p
      DataFrame(),             # dfsX0a
      DataFrame(),             # dfsX0e
      DataFrame()              # dfse
  )
end

# Optimal Control Problem (OCP) Results
mutable struct OCPResults{ T <: Number }
    tctr::Vector{Any}          # time vector for control
    tst::Vector{Any}           # time vector for state
    x::Matrix{JuMP.Variable}   # JuMP states
    u::Matrix{JuMP.Variable}   # JuMP controls
    X::Matrix{T}               # states
    U::Matrix{T}               # controls
    X0::Vector{T}              # initial states for OCP
    CS::Vector{T}              # costates
    tpolyPts::Vector{T}        # time sample points for polynomials  (NOTE these interpolated solutions were developed for calculating error, between them and a known Optimal solution)
    XpolyPts::Matrix{T}        # state evaluated using Lagrange/Linear polynomial
    CSpolyPts::Matrix{T}       # costate evaluated using Lagrange/Linear polynomial
    UpolyPts::Matrix{T}        # control evaluated using Lagrane/Linear polynomial
    AlltpolyPts::Vector{T}     # time sample points for polynomials
    AllXpolyPts::Matrix{T}     # state evaluated using Lagrange/Linear polynomial
    AllCSpolyPts::Matrix{T}    # costate evaluated using Lagrange/Linear polynomial
    AllUpolyPts::Matrix{T}     # control evaluated using Lagrane/Linear polynomial
    tpts::Vector{T}            # vector time sample points
    Xpts::Matrix{T}            # vector state sample points
    Upts::Matrix{T}            # vector control sample points
    CSpts::Matrix{T}           # vector costate sample points
    x0Con                      # handle for initial state constraints
    x0sCon
    xfCon                      # handle for final state constraints
    xfsCon
    dynCon                     # dynamics constraints
    constraint::Constraint{T}  # constraint handles and data
    evalNum::Int               # number of times optimization has been run
    iterNum                    # mics. data, perhaps an iteration number for a higher level algorithm
    status::Symbol             # optimization status
    tSolve::T                  # solve time for optimization
    objVal::T                  # objective function value
    dfs::Vector{DataFrame}     # results in DataFrame for plotting
    dfsOpt::DataFrame          # optimization information in DataFrame for plotting
    dfsCon                     # constraint data
    OCPResults(T::DataType=Float64) = new{T}(
        Vector{T}(),                        # time vector for control
        Vector{T}(),                        # time vector for state
        Matrix{JuMP.Variable}(undef,0,0),   # JuMP states
        Matrix{JuMP.Variable}(undef,0,0),   # JuMP controls
        Matrix{T}(undef,0,0),               # states
        Matrix{T}(undef,0,0),               # controls
        Vector{T}(),                        # initial states for OCP
        Vector{T}(),                        # costates
        Vector{T}(),                        # time sample points for polynomials
        Matrix{T}(undef,0,0),               # state evaluated using Lagrange/Linear polynomial
        Matrix{T}(undef,0,0),               # costate evaluated using Lagrange/Linear polynomial
        Matrix{T}(undef,0,0),               # control evaluated using Lagrane/Linear polynomial
        Vector{T}(),                        # time sample points for polynomials
        Matrix{T}(undef,0,0),               # state evaluated using Lagrange/Linear polynomial
        Matrix{T}(undef,0,0),               # costate evaluated using Lagrange/Linear polynomial
        Matrix{T}(undef,0,0),               # control evaluated using Lagrane/Linear polynomial
        Vector{T}(),                        # vector time sample points
        Matrix{T}(undef,0,0),               # vector state sample points
        Matrix{T}(undef,0,0),               # vector control sample points
        Matrix{T}(undef,0,0),               # vector costate sample points
        nothing,
        nothing,
        nothing,                            # handle for final state constraints
        nothing,
        nothing,                            # dynamics constraint
        Constraint(T), # constraint data
        1,                                  # current evaluation number
        [],                                 # mics. data, perhaps an iteration number for a higher level algorithm
        :undef,                             # optimization status
        parse(T,0),                         # solve time for optimization
        parse(T,0),                         # objective function value
        Vector{DataFrame}(),                # results in DataFrame for plotting
        DataFrame(),                        # optimization information in DataFrame for plotting
        []                                  # constraint data
    )
end

# Combined Results
mutable struct Results{T <: Number}
    ocp::OCPResults{T}
    ip::PlantResults{T}
    ep::PlantResults{T}
    resultsDir::AbstractString  # string that defines results folder
    mainDir::AbstractString     # string that defines main folder
    Results(T::DataType=Float64) = new{T}(
        OCPResults(T),
        PlantResults(T),
        PlantResults(T),
        joinpath(pwd(),"results"),  # string that defines results folder
        pwd()                       # string that defines main folder
    )
end

# Model-Predictive Contrl (MPC) Settings Class
mutable struct MPCSettings
    on::Bool
    mode::Symbol
    predictX0::Bool
    fixedTp::Bool
    IPKnown::Bool
    saveMode::Symbol
    maxSim::Int                 # maximum number of total MPC updates
    shiftX0::Bool               # a bool to indicate that a check on the feasibility of the linear constraints for the initial conditions should be performed
    lastOptimal::Bool           # a bool to indicate that the previous solution should be used if the current solution in not optimal
    printLevel::Int
    expandGoal::Bool            # bool to indicate if the goal needs to be expanded
    enlargeGoalTolFactor::Int   # scaling factor to enlare the goal
    onlyOptimal::Bool
    MPCSettings() = new(
        false,
        :OCP,
        false,
        true,
        true,
        :all,
        100,
        true,
        true,
        2,
        true,
        2,
        false
    )
end

# Optimal Control Problem (OCP) Settings Class
mutable struct OCPSettings
    solver::Solver              # solver information
    finalTimeDV::Bool
    integrationMethod::Symbol
    integrationScheme::Symbol
    save::Bool                  # bool for saving data
    reset::Bool                 # bool for reseting data
    evalConstraints::Bool       # bool for evaluating duals of the constraints
    evalCostates::Bool          # bool for evaluating costates
    tfMin::Any                  # minimum final time
    tfMax::Any                  # maximum final time
    tfOptimal::Any              # known optimal final time
    numInterpPts::Int           # number of points to sample polynomial running through collocation points
    cacheOnly::Bool             # bool for only caching the results when using optimize!()
    linearInterpolation::Bool   # bool for using linear interpolation even if integrationMethod ==:ps
    interpolationOn::Bool       # bool to indicate if user wants solution interpolated for them
    x0slackVariables::Bool
    xFslackVariables::Bool
    OCPSettings() = new(
        Solver(),           # default solver
        false,              # finalTimeDV
        :ts,                # integrationMethod
        :bkwEuler,          # integrationScheme
        true,               # bool for saving data
        false,              # bool for reseting data
        false,              # bool for evaluating duals of the constraints
        false,              # bool for evaluating costates
        0.001,              # minimum final time
        400.0,              # maximum final time
        false,              # known optimal final time
        250,                # number of points to sample polynomial running through collocation points
        false,              # bool for only caching the results when using optimize!()
        false,              # bool for using linear interpolation even if integrationMethod ==:ps
        false,              # bool to indicate if user wants solution interpolated for them
        false,
        false
    )
end

# Combined Settings
mutable struct Settings
    ocp::OCPSettings
    mpc::MPCSettings
    Settings() = new(
        OCPSettings(),
        MPCSettings()
    )
end

end # module
