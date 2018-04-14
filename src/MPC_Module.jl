module MPC_Module

using JuMP
using OrdinaryDiffEq
using DiffEqBase

include("Base.jl")
using .Base

export
      MPC,
      defineMPC!,
      defineIP!
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
 goal                 # goal location
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
              []
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
                   IPKnown::Bool=true,
                   saveMode::Symbol=:all,
                   maxSim::Int64=100,
                   goal=n.ocp.XF)
  n.s.mpc.on = true
  n.mpc::MPC = MPC()
  n.s.mpc.mode = mode
  n.s.mpc.predictX0 = predictX0
  n.s.mpc.fixedTp = fixedTp
  n.s.mpc.IPKnown = IPKnown
  n.s.mpc.saveMode = saveMode
  n.s.mpc.maxSim = maxSim
  n.mpc.v.goal = goal
  n.f.mpc.defined = true

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
    n.r.ip.X0a = copy(n.ocp.X0)
    n.mpc.ip.state.model = model
    n.mpc.ip.state.name = n.ocp.state.name
    n.mpc.ip.state.description = n.ocp.state.description
    n.mpc.ip.state.num = n.ocp.state.num
    n.mpc.ip.state.pts = n.ocp.state.pts

    n.mpc.ip.control.name = n.ocp.control.name
    n.mpc.ip.control.description = n.ocp.control.description
    n.mpc.ip.control.num = n.ocp.control.num
    n.mpc.ip.control.pts = n.ocp.control.pts

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
Date Create: 4/12/2018, Last Modified: 4/12/2018 \n
--------------------------------------------------------------------------------------\n
"""
function goalReached!(n)
 if ((n.r.ip.dfsplant[end][:x][end]-c["goal"]["x"])^2 + (n.r.ip.dfsplant[end][:y][end]-c["goal"]["yVal"])^2)^0.5 < c["goal"]["tol"]
   println("Goal Attained! \n"); n.f.mpc.goalReached=true;
 end
 return n.f.mpc.goalReached
end
"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/12/2018, Last Modified: 4/12/2018 \n
--------------------------------------------------------------------------------------\n
"""
function simMpc(c)

 n = initializeAutonomousControl(c)

 for ii = 1:n.mpc.s.maxSim
  println("Running model for the: ",n.mpc.v.evalNum + 1," time")
  updateAutoParams!(n,c)                 # update model parameters
  status = autonomousControl!(n)         # rerun optimization

  if n.r.ocp.status==:Optimal || n.r.ocp.status==:Suboptimal || n.r.ocp.status==:UserLimit
    println("Passing Optimized Signals to 3DOF Vehicle Model");
  elseif n.r.ocp.status==:Infeasible
    println("\n FINISH:Pass PREVIOUSLY Optimized Signals to 3DOF Vehicle Model \n"); break;
  else
    println("\n There status is nor Optimal or Infeaible -> create logic for this case!\n"); break;
  end

  n.mpc.t0_actual = (n.r.ocp.evalNum-1)*n.mpc.v.tex  # external so that it can be updated easily in PathFollowing

  # if the vehicle is very close to the goal sometimes the optimization returns with a small final time
  # and it can even be negative (due to tolerances in NLP solver). If this is the case, the goal is slightly
  # expanded from the previous check and one final check is performed otherwise the run is failed
  if getvalue(n.ocp.tf) < 0.01
    if ((n.r.ip.dfsplant[end][:x][end]-c["goal"]["x"])^2 + (n.r.ip.dfsplant[end][:y][end]-c["goal"]["yVal"])^2)^0.5 < 2*c["goal"]["tol"]
    println("Expanded Goal Attained! \n"); n.f.mpc.goalReached=true;
    break;
    else
    warn("Expanded Goal Not Attained! -> stopping simulation! \n"); break;
    end
  elseif getvalue(n.ocp.tf) < 0.5 # if the vehicle is near the goal => tf may be less then 0.5 s
    tf = (n.r.ocp.evalNum-1)*n.mpc.v.tex + getvalue(n.ocp.tf)
  else
    tf = (n.r.ocp.evalNum)*n.mpc.v.tex
  end

  if isequal(c["misc"]["model"],:ThreeDOFv2)
    U = n.r.ocp.U # TODO change to v1 for plant sim
  elseif isequal(c["misc"]["model"],:KinematicBicycle)
    #U = hcat(n.r.ocp.U[:,1],n.r.ocp.X[:,4])# TODO change to v1 for plant sim
    U = n.r.ocp.U
  end

  simPlant!(n;tf=tf,U=U)
  updateX0!(n)
  if n.r.ocp.evalNum==n.mpc.v.evalNum
    warn(" \n This is the last itteration! \n i.e. the maximum number of iterations has been reached while closing the loop; consider increasing (max_iteration) \n")
  end
  if ((n.r.ip.dfsplant[end][:x][end]-c["goal"]["x"])^2 + (n.r.ip.dfsplant[end][:y][end]-c["goal"]["yVal"])^2)^0.5 < c["goal"]["tol"]
    println("Goal Attained! \n"); n.f.mpc.goalReached=true;
    break;
  end
  if checkCrash(n,c,c["misc"]["sm2"];(:plant=>true))
    warn(" \n The vehicle crashed -> stopping simulation! \n"); break;
  end
 end
 return n
end

end # module
