module MPC_Module

using JuMP

include("Base.jl")
using .Base

export
      autonomousControl!,
      initializeMPC!,
      updateX0!,
      simPlant!,
      simModel,
      MPC


const simulationModes = [:OCP,:InternalPlant,:InternalExternalPlant,:ExternalPlant]

"""
initializeMPC(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""


#function configureMPC!(n::NLOpt; kwargs... )
function configureMPC!(n; kwargs... )
  kw = Dict(kwargs)

  # simulationMode
  if !haskey(kw,:simulationMode)
    n.mpc.simulationMode = :OCP
  else
    simulationMode = get(kw,:simulationMode,0)
      if simulationMode in simulationModes
        n.mpc.simulationMode = get(kw,:simulationMode,0)
      else
        error(simulationMode," is not in: ", simulationModes)
      end
  end

  if isequal(n.mpc.simulationMode,:OCP)

  elseif isequal(n.mpc.simulationMode,:internalPlant)

  elseif

end
















type MPC
  # models
  plantEquations
  modelEquations

  # constants
  tp::Any              # predication time (if finalTimeDV == true -> this is not known before optimization)
  tex::Float64         # execution horizon time
  max_iter::Int64      # maximum number of iterations
  numPlantControls::Int64      # number of states in plant
  numPlantStates::Int64        # number of controls in plant

  # variables
  goal_reached::Bool           # flag to indicate that goal has been reached
  t0_actual                    # actual initial time
  t0::Float64                  # mpc initial time
  tf::Float64                  # mpc final time
  t0_param::Any                # parameter for mpc t0
  X0p::Array{Float64,1}        # predicted controller initial states
  X0                           # array of all controller initial states (n.X0 leads because we are in simulation)
  X0p_plant::Array{Float64,1}  # predicted initial states
  X0_plant                     # array of all actual initial states (n.X0 leads because we are in simulation)

  # options
  PredictX0::Bool
  FixedTp::Bool
  Mode
  InternalPlantKnown::Bool    # bool to indicate if the internal plant model is assumed to be known
  ExternalPlant::Bool         # bool to indicate if there is also an external plant model
end

function MPC()
  MPC(
      Any,
      Any,
      Any,    # might be a variable tp
      0.0,
      0,
      0,      # number of states in plant
      0,      # number of controls in plant
      false,
      0.0,
      0.0,
      0.0,
      Any,
      Float64[],          # predicted initial state conditions
      [],
      Float64[],          # predicted actual initial state conditions
      [],
      true,
      true,
      true,
      false);
end

########################################################################################
# mpc functions
########################################################################################

"""
initializeMPC(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function initializeMPC!(n;
                        numPlantStates::Int64=0,
                        numPlantControls::Int64=0,
                        X0_plant=fill(NaN,numPlantStates,),
                        InternalPlantKnown::Bool=true,
                        FixedTp::Bool=true,
                        PredictX0::Bool=true,
                        tp::Any=NaN,tex::Float64=0.5,
                        max_iter::Int64=10)
  if n.mpc.t0_param!=Any
    error("\n initializeMPC!() must be called before OCPdef!() \n")
  end
  # validate input
  if numPlantControls <= 0
      error("numPlantControls must be > 0","\n");
  end
  if numPlantStates <= 0
      error("numPlantStates must be > 0","\n");
  end
  if length(X0_plant) != numPlantStates
    error(string("\n Length of X0_plant must match number of plant states \n"));
  end

  n.s.MPC = true
  n.mpc::MPC = MPC()
  n.mpc.InternalPlantKnown = InternalPlantKnown
  n.mpc.numPlantControls = numPlantControls
  n.mpc.numPlantStates = numPlantStates
  n.mpc.X0_plant = [X0_plant]
  n.mpc.FixedTp = FixedTp
  n.mpc.PredictX0 = PredictX0
  n.mpc.tp = tp
  n.mpc.tex =tex
  n.mpc.max_iter = max_iter
  n.mpc.t0 = 0.0
  n.mpc.t0_actual = 0.0
  n.mpc.tf = tex
  return nothing
end

"""
# for simulating the model of the plant given control commands
simModel(n,X0,t,U,t0,tf)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simModel(n,X0,t,U,t0,tf)
  n.mpc.modelEquations(n,X0,t,U,t0,tf)
end

"""
updateStates!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/17/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function updateStates!(n)

  # update states based off of n.mpc.X0p
  for st in 1:n.numStates   # update states based off of n.mpc.X0p
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(n.r.x0_con[st,1], (n.mpc.X0p[st]+n.X0_tol[st]));
      JuMP.setRHS(n.r.x0_con[st,2],-(n.mpc.X0p[st]-n.X0_tol[st]));
    else
      JuMP.setRHS(n.r.x0_con[st],n.mpc.X0p[st]);
    end
  end
  return nothing
end

"""
updateX0!(n,X0;(:userUpdate=>true))    # user defined update of X0
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateX0!(n,args...;kwargs...)
  kw = Dict(kwargs);
  # check to see how the initial states are being updated
  if !haskey(kw,:userUpdate); kw_ = Dict(:userUpdate => false); userUpdate = get(kw_,:userUpdate,0);
  else; userUpdate = get(kw,:userUpdate,0);
  end

  if userUpdate
    X0=args[1];
    if length(X0)!=n.numStates
      error(string("\n Length of X0 must match number of states \n"));
    end
    n.X0 = X0
  else # update using the current location of plant
    for st in 1:n.numStates
      n.X0[st] = n.r.dfs_plant[end][n.state.name[st]][end]
    end
  end
  updateStates!(n)
  append!(n.mpc.X0,[copy(n.X0)])
  return nothing
end

"""
simPlant(n)
# consider X0=n.mpc.X0[(n.r.eval_num)] when plant and controller are different
# NOTE if previous solution was Infeasible it will just pass the begining of the last Optimal solution again
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/09/2018 \n
--------------------------------------------------------------------------------------\n
"""
function simPlant!(n;X0=n.X0,t=n.r.t_ctr+n.mpc.t0,U=n.r.U,t0=n.mpc.t0_actual,tf=n.r.eval_num*n.mpc.tex)
  sol = n.mpc.plantEquations(n,X0,t,U,t0,tf)
  plant2dfs!(n,sol)  #TODO consider passing U, t0 etc.. and only run this if saving data
  return nothing
end

"""
# TODO eventually the "plant" will be different from the "model"
t0p=predictX0(n)
# also shifts the initial time, but this is accounted for in updateMPC()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function predictX0!(n)
  if n.mpc.tf==n.mpc.t0
    error("n.mpc.tf==n.mpc.t0")
  end

  if n.mpc.FixedTp
    t0p = n.mpc.tex
  elseif n.r.eval_num == 0
    t0p = n.t0
  else
    t0p = r.dfs_opt[end][:t_solve][1]
  end

  if n.r.eval_num != 0
    # based off of "current X0". Even though we may have the next X0 we should not (i.e.look at @show length(n.mpc.X0)). It is because it is a simulation (in reality they would be running in parallel)
    sol = simModel(n,n.mpc.X0[(n.r.eval_num)],n.r.t_ctr+n.mpc.t0,n.r.U,n.mpc.t0,n.mpc.t0+t0p)
    n.mpc.X0p = sol(n.mpc.t0+t0p)[:]
  else
    n.mpc.X0p = n.X0  # assuming vehicle did not move
  end

  return t0p
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function mpcUpdate!(n)
  if n.mpc.PredictX0      # predict where X0 will be when optimized signal is actually sent to the vehicle
    t0p = predictX0!(n)   # predicted start time -> important for time varying constraints
  else
    t0p = 0
    n.mpc.X0p = n.mpc.X0[end]  # current "known" plant states  TODO make sure the plant is never simulated ahead
  end
  n.mpc.t0 = copy(n.mpc.t0_actual+t0p)      # there are two different time scales-> the plant leads by t0p
  setvalue(n.mpc.t0_param,copy(n.mpc.t0))   # update for time varying constraints
  n.mpc.tf = copy(n.mpc.t0+n.mpc.tex)
  if n.mpc.t0_param!=Any
    setvalue(n.mpc.t0_param,copy(n.mpc.t0))
  end
  return nothing
end

"""
status = autonomousControl!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/1/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function autonomousControl!(n)
  if n.r.eval_num != 0
    mpcUpdate!(n)
    updateStates!(n)
  end
  optimize!(n)
  return n.r.status
end


end # module
