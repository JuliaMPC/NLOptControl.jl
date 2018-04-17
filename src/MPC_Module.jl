module MPC_Module

using JuMP
using OrdinaryDiffEq
using DiffEqBase

include("Base.jl")
using .Base

export
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
     simMPC!

########################################################################################
# MPC types
########################################################################################

type IP
 control::Control
 state::State
end

function IP()
 IP(
  Control(),
  State()
  )
end

type EP
 control::Control
 state::State
end

function EP()
 EP(
  Control(),
  State()
 )
end

type MPCvariables
 # variables
 t::Float64           # current simulation time (s)
 tp::Any              # prediction time (if finalTimeDV == true -> this is not known before optimization)
 tex::Float64         # execution horizon time
 t0Actual                    # actual initial time TODO ?
 t0::Float64                  # mpc initial time TODO ?
 tf::Float64                  # mpc final time TODO ?
 t0Param::Any        # parameter for mpc t0  TODO ?
 evalNum::Int64       # parameter for keeping track of number of MPC evaluations
 goal                 # goal location w.r.t OCP
 goalTol             # tolerance on goal location
 initOptNum::Int64  # number of initial optimization
 previousSolutionNum::Int64  # number of times the previous solution should be used
end

function MPCvariables()
 MPCvariables(
              0.0,    # t
              Any,    # tp (might be a variable)
              0.5,    # tex
              0.0,
              0.0,
              0.0,
              Any,
              0,
              [],
              [],
              3,
              3
 )
end

type MPC
 v::MPCvariables
 ip::IP
 ep::EP
end

function MPC()
 MPC(
     MPCvariables(),
     IP(),
     EP()
     )
end

########################################################################################
# MPC functions
########################################################################################
"""
defineMPC!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/8/2018 \n
--------------------------------------------------------------------------------------\n
"""
function defineMPC!(n;
                   mode::Symbol=:OCP,
                   predictX0::Bool=true,
                   fixedTp::Bool=true,
                   tp::Any=Any,
                   tex::Float64=0.5,
                   IPKnown::Bool=true,
                   saveMode::Symbol=:all,
                   maxSim::Int64=100,
                   goal=n.ocp.XF,
                   goalTol=0.1*abs.(n.ocp.X0 - n.ocp.XF),
                   lastOptimal::Bool=true,
                   printLevel::Int64=2)
 n.s.mpc.on = true
 n.mpc::MPC = MPC()
 n.s.mpc.mode = mode
 n.s.mpc.predictX0 = predictX0
 n.s.mpc.fixedTp = fixedTp
 n.mpc.v.tp = tp
 n.mpc.v.tex = tex
 n.s.mpc.IPKnown = IPKnown
 n.s.mpc.saveMode = saveMode
 n.s.mpc.maxSim = maxSim
 n.mpc.v.goal = goal
 n.mpc.v.goalTol = goalTol
 n.s.mpc.lastOptimal = lastOptimal
 n.s.mpc.printLevel = printLevel
 n.f.mpc.defined = true
 return nothing
end

"""
# TODO consider letting user pass options
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
  for k in 1:n.mpc.v.initOptNum # initial optimization (s)
   status = optimize!(n)
   if status==:Optimal; break; end
  end
  # defineSolver!(n,solverConfig(c)) # modifying solver settings NOTE currently not in use
  n.s.ocp.save = true  # NOTE set to false if running in parallel to save time
  n.s.ocp.cacheOnly = false
  n.s.ocp.evalConstraints = true # NOTE set to true to investigate infeasibilities
  return nothing
end

"""
# add a mode that solves as quickly as possible
# consider using the IP always.
defineModel!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/12/2018, Last Modified: 4/12/2018 \n
--------------------------------------------------------------------------------------\n
"""
function defineIP!(n,model;stateNames=[],controlNames=[],X0a=[])

   if isequal(n.s.mpc.mode,:OCP) # this function is called automatically for this mode
    if !isempty(stateNames)
     error("stateNames are set automatically for :mode == :OCP and cannot be provided.")
    end
    if !isempty(controlNames)
     error("controlNames are set automatically for :mode == :OCP and cannot be provided.")
    end
    if !isempty(X0a)
     error("X0a is set automatically for :mode == :OCP and cannot be provided.")
    end
    n.r.ip.X0a = copy(n.ocp.X0)  # NEED to append time
    n.mpc.ip.state.model = model
    n.mpc.ip.state.name = n.ocp.state.name
    n.mpc.ip.state.description = n.ocp.state.description
    n.mpc.ip.state.num = n.ocp.state.num
    n.mpc.ip.state.pts = n.ocp.state.pts

    n.mpc.ip.control.name = n.ocp.control.name
    n.mpc.ip.control.description = n.ocp.control.description
    n.mpc.ip.control.num = n.ocp.control.num
    n.mpc.ip.control.pts = n.ocp.control.pts

    # add X0 t0 plant dfs
    n.r.ip.plant[:t] = n.mpc.v.t0
    for st in 1:n.mpc.ip.state.num
      n.r.ip.plant[n.mpc.ip.state.name[st]] = copy(n.ocp.X0)[st]
    end
    for ctr in 1:n.mpc.ip.control.num
      n.r.ip.plant[n.mpc.ip.control.name[ctr]] = 0
    end

   elseif isequal(n.s.mpc.mode,:IP)
    if isempty(stateNames)
     error("unless :mode == :OCP the stateNames must be provided.")
    end
    if isempty(controlNames)
     error("unless :mode == :OCP the controlNames must be provided.")
    end
    if isempty(X0a)
     error("unless :mode == :OCP X0a must be provided.")
    end
    if isempty(model)
     error("A model needs to be passed for the IP mode.")
    else
    if isequal(length(X0a),length(stateNames))
      error(string("\n Length of X0a must match length(stateNames) \n"))
    end

     n.mpc.ip.state::State = State() # reset
     n.mpc.ip.state.num = length(stateNames)
     for i in 1:n.mpc.ip.state.num
       if stateNames[i]==:xxx
         error("xxx is OFF limits for a state name; please choose something else. \n")
       end
       push!(n.mpc.ip.state.name,stateNames[i])
     end

     n.mpc.ip.control::Control = Control() # reset
     n.mpc.ip.control.num = length(controlNames)
     for i in 1:n.mpc.ip.control.num
       if controlNames[i]==:xxx
         error("xxx is OFF limits for a control name; please choose something else. \n")
       end
       push!(n.mpc.ip.control.name,controlNames[i])
     end
     n.mpc.r.ip.X0a = X0a
     n.mpc.ip.state.model = model # TODO validate type of model
    end
   elseif isequal(n.s.mpc.mode,:EP)
    error("not setup for :EP")
   else
    error("n.mpc.s.mode = ",n.s.mpc.mode," not defined." )
   end

   # consider calling mapNames
   return nothing
end


"""
mapNames!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/9/2018, Last Modified: 4/12/2018 \n
--------------------------------------------------------------------------------------\n
"""
function mapNames!(n)
  if isequal(n.s.mpc.mode,:IP)
    s1 = n.ocp.state.name
    c1 = n.ocp.control.name
    s2 = n.mpc.stateIP.name
    c2 = n.mpc.controlIP.name
  elseif isequal(n.s.mpc.mode,:EP)
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

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/09/2018 \n
--------------------------------------------------------------------------------------\n
"""
function simIPlant!(n)
  X0 = currentIPState(n)[1]
  t = n.r.ocp.tctr
  U = n.r.ocp.U  # NOTE this is OK for the :OCP case
  t0 = n.mpc.v.t
  tf = n.mpc.v.t + n.mpc.v.tex
  sol, U = n.mpc.ip.state.model(n,X0,t,U,t0,tf)
  plant2dfs!(n,sol,U)
  return nothing
end

"""
updateX0!(n,X0;(:userUpdate=>true))    # user defined update of X0
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 4/12/2018 \n
--------------------------------------------------------------------------------------\n
"""
function updateX0!(n,args...;kwargs...)
  kw = Dict(kwargs)

  if isequal(mode,:OCP) # push the vehicle along
   for st in 1:n.ocp.state.num
     n.ocp.X0[st] = n.r.ip.dfsplant[end][n.ocp.state.name[st]][end]
   end
  elseif isequal(mode,:IP)
    error(":IP function not ready")
  elseif isequal(mode,:EP)
     error(":EP function not ready")
     X0 = args[1]
     if length(X0)!=n.ocp.state.num
       error(string("\n Length of X0 must match number of states \n"));
     end
     n.ocp.X0 = X0
  else
    error("mode must be :OCP, :IP, or :EP")
  end
  updateStates!(n)
  append!(n.r.ip.X0a,[copy(n.ocp.X0)])
  return nothing
end


"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/08/2018, Last Modified: 4/08/2018 \n
--------------------------------------------------------------------------------------\n
"""
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

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/08/2018, Last Modified: 4/08/2018 \n
--------------------------------------------------------------------------------------\n
"""
function updateOCPState!(n)
  if n.s.mpc.shiftX0 # TODO consider saving linear shifting occurances
    for st in 1:n.ocp.state.num
      if n.ocp.X0[st] < n.ocp.XL[st]
        n.ocp.X0[st] = n.ocp.XL[st]
      end
      if n.ocp.X0[st] > n.ocp.XU[st]
        n.ocp.X0[st] = n.ocp.XU[st]
      end
    end
  end
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


"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/08/2018, Last Modified: 4/08/2018 \n
--------------------------------------------------------------------------------------\n
"""
function goalReached!(n)
  # TODO deal with NaNs
  if isequal(n.s.mpc.mode,:OCP)
    X = currentIPState(n)[1]
  else
    #TODO
  end
  if all((abs.(X - n.mpc.v.goal) .<= n.mpc.v.goalTol))
   if isequal(n.s.mpc.printLevel,2)
    println("Goal Attained! \n")
   end
    n.f.mpc.goalReached = true
  end
 return n.f.mpc.goalReached
end
# if the vehicle is very close to the goal sometimes the optimization returns with a small final time
# and it can even be negative (due to tolerances in NLP solver). If this is the case, the goal is slightly
# expanded from the previous check and one final check is performed otherwise the run is failed
#if getvalue(n.ocp.tf) < 0.01
#  if ((n.r.ip.dfplant[end][:x][end]-c["goal"]["x"])^2 + (n.r..ip.dfplant[end][:y][end]-c["goal"]["yVal"])^2)^0.5 < 2*c["goal"]["tol"]
#  println("Expanded Goal Attained! \n"); n.f.mpc.goal_reached=true;
#  break;
#  else
#  warn("Expanded Goal Not Attained! -> stopping simulation! \n"); break;
#  end
#elseif getvalue(n.ocp.tf) < 0.5 # if the vehicle is near the goal => tf may be less then 0.5 s
#  tf = (n.r.evalNum-1)*n.mpc.v.tex + getvalue(n.ocp.tf)
#else
#  tf = (n.r.evalNum)*n.mpc.v.tex
#end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/06/2018, Last Modified: 4/08/2018 \n
--------------------------------------------------------------------------------------\n
"""
function simMPC!(n;updateFunction::Any=[])
  for ii = 1:n.s.mpc.maxSim
    if isequal(n.s.mpc.printLevel,2)
     println("Running model for the: ",n.mpc.v.evalNum + 1," time")
    end
    #############################
    # (A) and (B) in "parallel"
    #############################
    # (A) solve OCP
    if !isequal(typeof(updateFunction),Array{Any,1})
      updateFunction(n)
    end

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

    setvalue(n.ocp.t0,copy(n.mpc.v.t))

    # check to see if the goal has been reached
    if goalReached!(n); break; end
  end
end

end # module
