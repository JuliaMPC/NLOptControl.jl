########################################################################################
# mpc functions
########################################################################################

"""
defineMPC!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function defineMPC!(n;
                        numPlantStates::Int64=0,
                        numPlantControls::Int64=0,
                        X0_plant=fill(NaN,numPlantStates,),
                        InternalPlantKnown::Bool=true,
                        FixedTp::Bool=true,
                        PredictX0::Bool=true,
                        tp::Any=NaN,tex::Float64=0.5,
                        max_iter::Int64=10)
  if n.mpc.v.t0Param!=Any
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

  n.s.mpc.on = true
  n.mpc::MPC = MPC()
  n.mpc.InternalPlantKnown = InternalPlantKnown
  n.mpc.numPlantControls = numPlantControls
  n.mpc.numPlantStates = numPlantStates
  n.mpc.X0_plant = [X0_plant]
  n.mpc.FixedTp = FixedTp
  n.mpc.PredictX0 = PredictX0
  n.mpc.v.tp = tp
  n.mpc.v.tex =tex
  n.mpc.v.evalNum = max_iter
  n.mpc.v.t0 = 0.0
  n.mpc.v.t0Actual = 0.0
  n.mpc.tf = tex
  return nothing
end

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

  else
  end

end





"""
updateStates!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/17/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function updateStates!(n)

  # NOTE assuming known  IP
  # update states based off of n.r.ip.X0p
  for st in 1:n.ocp.state.num   # update states based off of n.r.ip.X0p
    if any(!isnan(n.ocp.X0_tol[st]))
      JuMP.setRHS(n.r.ocp.x0Con[st,1], (n.r.ip.X0p[st]+n.ocp.X0_tol[st]));
      JuMP.setRHS(n.r.ocp.x0Con[st,2],-(n.r.ip.X0p[st]-n.ocp.X0_tol[st]));
    else
      JuMP.setRHS(n.r.ocp.x0Con[st],n.r.ip.X0p[st]);
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
    if length(X0)!=n.ocp.state.num
      error(string("\n Length of X0 must match number of states \n"));
    end
    n.ocp.X0 = X0
  else # update using the current location of plant
    for st in 1:n.ocp.state.num
      n.ocp.X0[st] = n.r.ip.dfsplant[end][n.ocp.state.name[st]][end]
    end
  end
  updateStates!(n)
  append!(n.r.ip.X0a,[copy(n.ocp.X0)])
  return nothing
end

"""
simPlant(n)
# consider X0=n.r.ip.X0a[(n.r.ocp.evalNum)] when plant and controller are different
# NOTE if previous solution was Infeasible it will just pass the begining of the last Optimal solution again
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/09/2018 \n
--------------------------------------------------------------------------------------\n
"""
function simPlant!(n;X0=n.ocp.X0,t=n.r.ocp.tctr+n.mpc.v.t0,U=n.r.ocp.U,t0=n.mpc.v.t0Actual,tf=n.r.ocp.evalNum*n.mpc.v.tex)
  sol = n.mpc.ip.state.model(n,X0,t,U,t0,tf)
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
  if n.mpc.tf==n.mpc.v.t0
    error("n.mpc.tf==n.mpc.v.t0")
  end

  if n.mpc.FixedTp
    t0p = n.mpc.v.tex
  elseif n.r.ocp.evalNum == 0
    t0p = n.t0
  else
    t0p = r.dfs_opt[end][:t_solve][1]
  end

  if n.r.ocp.evalNum != 0
    # based off of "current X0". Even though we may have the next X0 we should not (i.e.look at @show length(n.r.ip.X0a)). It is because it is a simulation (in reality they would be running in parallel)
    sol = simModel(n,n.r.ip.X0a[(n.r.ocp.evalNum)],n.r.ocp.tctr+n.mpc.v.t0,n.r.ocp.U,n.mpc.v.t0,n.mpc.v.t0+t0p)
    n.r.ip.X0p = sol(n.mpc.v.t0+t0p)[:]
  else
    n.r.ip.X0p = n.ocp.X0  # assuming vehicle did not move
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
    n.r.ip.X0p = n.r.ip.X0a[end]  # current "known" plant states  TODO make sure the plant is never simulated ahead
  end
  n.mpc.v.t0 = copy(n.mpc.v.t0Actual+t0p)      # there are two different time scales-> the plant leads by t0p
  setvalue(n.mpc.v.t0Param,copy(n.mpc.t0))   # update for time varying constraints
  n.mpc.tf = copy(n.mpc.t0+n.mpc.v.tex)
  if n.mpc.v.t0Param!=Any
    setvalue(n.mpc.v.t0Param,copy(n.mpc.t0))
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
  if n.r.ocp.evalNum != 0
    mpcUpdate!(n)
    updateStates!(n)
  end
  optimize!(n)
  return n.r.ocp.status
end

"""
mapNames!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/9/2018, Last Modified: 4/9/2018 \n
--------------------------------------------------------------------------------------\n
"""
function mapNames!(n)
  if isequal((n.mpc.s.simulationMode,:IP)
    s1 = n.ocp.state.name
    c1 = n.ocp.control.name
    s2 = n.mpc.stateIP.name
    c2 = n.mpc.controlIP.name
  elseif isequal((n.mpc.s.simulationMode,:EP)
    error(":EP function not ready")
  else
    error("mode must be either :IP or :EP")
  end

  m = []
  # go through all states in OCP
  idxOCP = 1
  for var in s1
    # go through all states in IP
    idxIP = find(var.==s2)
    if !isempty(idxIP)
      push!(m, [var; :stOCP; idxOCP; :stIP; idxIP[1]])
    end

    # go through all controls in IP
    idxIP = find(var.==c2)
    if !isempty(idxIP)
      push!(m, [var; :stOCP; idxOCP; :ctrIP; idxIP[1]])
    end
    idxOCP = idxOCP + 1
  end

  # go through all controls in OCP
  idxOCP = 1
  for var in c1
    # go through all states in IP
    idxIP = find(var.==s2)
    if !isempty(idxIP)
      push!(m, [var; :ctrOCP; idxOCP; :stIP; idxIP[1]])
    end

    # go through all controls in IP
    idxIP = find(var.==c2)
    if !isempty(idxIP)
      push!(m, [var; :ctrOCP; idxOCP; :ctrIP; idxIP[1]])
    end
    idxOCP = idxOCP + 1
  end

  if isequal(mode,:IP)
    n.mpc.mIP = m
  elseif isequal(mode,:EP)
    error(":EP function not ready")
  else
    error("mode must be either :IP or :EP")
  end

  return nothing
end

using VehicleModels

# Bryson Denham

This problem can be found [here](http://www.gpops2.com/Examples/Bryson-Denham.html).


# OCP
using NLOptControl, PrettyPlots
n=define(numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,NaN],XU=[1/9,NaN]);

states!(n,[:x1,:x2])
controls!(n,[:u1])

dx=[:(x2[j]),:(u1[j])]
dynamics!(n,dx)
configure!(n;(:Nck=>[100]),(:finalTimeDV=>true));
obj=integrate!(n,:(0.5*u1[j]^2));
@NLobjective(n.ocp.mdl,Min,obj);
optimize!(n);

# Define IP (maybe IP can just be OCP)
defineMPC!(n)

sn = [:x1,:x2]
cn = [:u1]
X0a = []
defineModel!(n)



# OCP
#[:x,:y,:psi,:ux]
#[:sa,:ax]

# IP
#[:x,:y,:v,:r,:psi,:sa,:ux,:ax]
#[:sr,:jx]


##### funcs.jl
function simIPlant!(n)
  X0 = currentIPState(n)[1]
  t = n.r.ocp.tctr
  U = n.r.ocp.U  # NOTE this is OK for the :OCP case
  t0 = n.mpc.v.t
  tf = n.mpc.v.t + n.mpc.v.tex
  sol, U = n.mpc.ip.state.model(X0,t,U,t0,tf)
  plant2dfs!(n,sol,U)

  return nothing
end

function currentIPState(n)
  if isempty(n.r.ip.plant)
    error("there is no data in n.r.ip.plant")
  end

  # even though may have solution for plant ahead of time
  # can only get the state up to n.mpc.v.t
  idx = find((n.mpc.v.t - n.r.ip.plant[:t]) .>= 0)
  if isempty(idx)
    error("(n.mpc.v.t - n.r.ip.plant[:t]) .>= 0) is empty.")
  else
    X0 = [zeros(n.mpc.ip.state.num),n.mpc.v.t]
    for st in 1:n.mpc.ip.state.num
      X0[1][st] = n.r.ip.plant[n.mpc.ip.state.name[st]][idx[end]]
    end
  end
  return X0
end

function updateOCPState!(n)
  # update states with n.ocp.X0
  for st in 1:n.ocp.state.num
    if any(!isnan(n.ocp.X0_tol[st]))
      JuMP.setRHS(n.r.ocp.x0Con[st,1], (n.ocp.X0[st]+n.ocp.X0_tol[st]));
      JuMP.setRHS(n.r.ocp.x0Con[st,2],-(n.ocp.X0[st]-n.ocp.X0_tol[st]));
    else
      JuMP.setRHS(n.r.ocp.x0Con[st],n.ocp.X0[st]);
    end
  end
  return nothing
end


function goal!(n)
  # TODO deal with NaNs
  if isequal(n.s.mpc.mode,:OCP)
    X = currentIPState(n)[1]
  else
    #TODO
  end

  if all((abs.(X - n.mpc.v.goal) .<= n.mpc.v.goalTol))
    println("Goal Attained! \n")
    n.f.mpc.goalReached = true
  end

 return n.f.mpc.goalReached
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/08/2018, Last Modified: 4/08/2018 \n
--------------------------------------------------------------------------------------\n
"""
function initOpt!(n)
  n.s.ocp.save = false
  n.s.mpc.on = false
  n.s.ocp.evalConstraints = false
  n.s.ocp.cacheOnly = true

  if n.s.ocp.save
   warn("saving initial optimization results where functions where cached!")
  end
  for k in 1:3 # initial optimization (s)  TODO make this a parameter
   status = optimize!(n);
   if status==:Optimal; break; end
  end
  # defineSolver!(n,solverConfig(c)) # modifying solver settings NOTE currently not in use
  n.s.ocp.save = true  # NOTE set to false if running in parallel to save time
  n.s.ocp.cacheOnly = false
  n.s.ocp.evalConstraints = false # NOTE set to true to investigate infeasibilities
  return nothing
end

function simMPC!(n)
  for ii = 1:n.s.mpc.maxSim
   println("Running model for the: ",n.mpc.v.evalNum + 1," time")
    #############################
    # (A) and (B) in "parallel"
    #############################
    # (A) solve OCP
    if !n.s.mpc.predictX0 #  use the current known plant state to update OCP
      # X0p is simply the current known location of the plant
      push!(n.r.ip.X0p,currentIPState(n))

      # need to map n.r.ip.X0p to n.X0 (states may be different)
      # NOTE for the :OCP mode this is OK
      if isequal(n.s.mpc.mode,:OCP)
        n.ocp.X0 = n.r.ip.X0p[end][1] # the n.ocp. structure is for running things
        push!(n.r.ocp.X0,n.ocp.X0)  # NOTE this may be for saving data
      else
        error("TODO")
      end
    else
      error("TODO")
    end
    updateOCPState!(n)
    optimize!(n)

    # (B) simulate plant
    simIPlant!(n) # the plant simulation time will lead the actual time

    # advance the actual time
    n.mpc.v.t = n.mpc.v.t + n.mpc.v.tex
    n.mpc.v.evalNum = n.mpc.v.evalNum + 1

    #NOTE want to try and avoid having time not start at 0 in OCP
    setvalue(n.ocp.t0,copy(n.mpc.v.t))

    goal!(n)     # check to see if the goal is in range
    if n.f.mpc.goalReached
      break
    end
  end
end
