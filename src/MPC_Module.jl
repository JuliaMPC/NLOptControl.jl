module MPC_Module

using JuMP
using Ranges

include("Base.jl")
using .Base

export
      autonomousControl!,
      initializeMPC!,
      driveStraight!,
      updateX0!,
      simPlant!,
      simModel,
      MPC

type MPC
  # models
  plantEquations
  modelEquations

  # constants
  tp::Any          # predication time (if finalTimeDV == true -> this is not known before optimization)
  tex::Float64         # execution horizon time
  max_iter::Int64      # maximum number of iterations

  # variables
  goal_reached::Bool           # flag  to indicate that goal has been reached
  t0_actual                    # actual intial time
  t0::Float64                  # mpc inital time
  tf::Float64                  # mpc final time
  t0_param::Any                # parameter for mpc t0
  X0p::Array{Float64,1}        # predicted initial states
  X0                           # array of all actual intial states (n.X0 leads because we are in simulation)

  # options
  PredictX0::Bool
  FixedTp::Bool
end

function MPC()
  MPC(
      Any,
      Any,
      Any,  # might be a variable tp
      0.0,
      0,
      false,
      0.0,
      0.0,
      0.0,
      Any,
      Float64[],          # predicted intial state conditions
      [],
      true,
      true);
end

########################################################################################
# mpc functions
########################################################################################

"""
initializeMPC(n)

initializeMPC!(n;FixedTp=c.m.FixedTp,PredictX0=c.m.PredictX0,tp=c.m.tp,tex=copy(c.m.tex),max_iter=c.m.mpc_max_iter);

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function initializeMPC!(n;FixedTp::Bool=true,PredictX0::Bool=true,tp::Any=NaN,tex::Float64=0.5,max_iter::Int64=10)
  if n.mpc.t0_param!=Any
    error("\n initializeMPC!() must be called before OCPdef!() \n")
  end
  n.s.MPC=true;
  n.mpc::MPC=MPC();
  n.mpc.FixedTp=FixedTp;
  n.mpc.PredictX0=PredictX0;
  n.mpc.tp=tp;
  n.mpc.tex=tex;
  n.mpc.max_iter=max_iter;
  n.mpc.t0=0.0;
  n.mpc.t0_actual=0.0;
  n.mpc.tf=tex;
  return nothing
end
#initializeMPC!(n,r,c)=initializeMPC!(n,r;FixedTp=c.m.FixedTp,PredictX0=c.m.PredictX0,tp=c.m.tp,tex=c.m.tex,max_iter=c.m.mpc_max_iter);

"""
# for simulating the model of the plant given control commands
simModel(n,X0,t,U,t0,tf)

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simModel(n,X0,t,U,t0,tf)
  n.mpc.modelEquations(n,X0,t,U,t0,tf);
end

"""
updateStates!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/17/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateStates!(n)
  for st in 1:n.numStates   # update states
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(n.r.x0_con[st,1], (n.X0[st]+n.X0_tol[st]));
      JuMP.setRHS(n.r.x0_con[st,2],-(n.X0[st]-n.X0_tol[st]));
   else
     JuMP.setRHS(n.r.x0_con[st],n.X0[st]);
   end
  end
  return nothing
end

"""
# NOTE: this function is called inside simPlant()
updateX0!(n,X0;(:userUpdate=>true))    # user defined update of X0
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateX0!(n,args...;kwargs...)
  kw = Dict(kwargs);
  # check to see how the intial states are being updated
  if !haskey(kw,:userUpdate); kw_ = Dict(:userUpdate => false); userUpdate = get(kw_,:userUpdate,0);
  else; userUpdate = get(kw,:userUpdate,0);
  end

  if userUpdate
    X0=args[1];
    if length(X0)!=n.numStates
      error(string("\n Length of X0 must match number of states \n"));
    end
    n.X0=X0;
  else # update using the current location of plant
    for st in 1:n.numStates
      n.X0[st]=n.r.dfs_plant[end][n.state.name[st]][end];
    end
  end
  updateStates!(n)
  append!(n.mpc.X0,[copy(n.X0)])
  return nothing
end

"""
simPlant(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 6/27/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simPlant!(n;X0=n.X0,t=n.r.t_ctr+n.mpc.t0,U=n.r.U,t0=n.mpc.t0_actual,tf=n.r.eval_num*n.mpc.tex)
  sol=n.mpc.plantEquations(n,X0,t,U,t0,tf);
  plant2dfs!(n,sol);  #TODO consider passing U, t0 etc..
  t0p=updateX0!(n);
  return nothing
end

"""
# TODO eventually the "plant" will be different from the "model"
t0p=predictX0(n)
# also shifts the intial time, but this is accounted for in updateMPC()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function predictX0!(n)
  if n.mpc.tf==n.mpc.t0
    error("n.mpc.tf==n.mpc.t0")
  end

  if n.mpc.FixedTp || r.eval_num==1
    t0p=n.mpc.tex;
  else
    t0p=r.dfs_opt[end][:t_solve][1]
  end
  if n.r.eval_num != 1
    # based off of "current X0". Even though we may have the next X0 we should not (i.e.look at @show length(n.mpc.X0)). It is because it is a simulation (in reality they would be running in parallel)
    sol=simModel(n,n.mpc.X0[n.r.eval_num],n.r.t_ctr+n.mpc.t0,n.r.U,n.mpc.t0,n.mpc.t0+t0p)
    n.mpc.X0p=sol(n.mpc.t0+t0p)[:];
  else # using this with driveStraight!()
    n.mpc.X0p = n.X0  # NOTE assuming the vehicle did not move
  end

  return t0p
end
"""
# this function currently only works when the "actual model"==simPlant()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/9/2017, Last Modified: 6/27/2017 \n
--------------------------------------------------------------------------------------\n
"""
function driveStraight!(n;t0::Float64=n.mpc.t0,tf::Float64=n.mpc.tf)
  # add these signals to r so that they can be used for predictions during optimization
  n.r.U=0*Matrix{Float64}(n.numControlPoints,n.numControls);
  n.r.t_ctr=Vector(Ranges.linspace(t0,tf,n.numControlPoints)); # gives a bunch of points
  postProcess!(n;(:Init=>true));                               # to make solutions line up

  # simulate the "actual vehicle" response
  simPlant!(n;X0=n.X0,t=n.r.t_ctr+n.mpc.t0,U=n.r.U,t0=t0,tf=tf)

  # update the inital time for the optimization
  #setvalue(n.mpc.t0_param,copy(tf));
  mpcUpdate!(n)

  return nothing
end


"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 8/1/2017 \n
--------------------------------------------------------------------------------------\n
"""
function mpcUpdate!(n)
  if n.mpc.PredictX0    # predict where X0 will be when optimized signal is actually sent to the vehicle
    t0p=predictX0!(n);  # predicted start time -> important for time varying constraints
  else
    t0p=0; n.mpc.X0p=n.mpc.X0[end];  # current "known" plant states  TODO make sure the plant is never simulated ahead
  end
  n.mpc.t0=copy(n.mpc.t0_actual+t0p);      # there are two different time scales-> the plant leads by t0p
  setvalue(n.mpc.t0_param,copy(n.mpc.t0)); # update for time varying constraints
  n.mpc.tf=copy(n.mpc.t0+n.mpc.tex);
  if n.mpc.t0_param!=Any
    setvalue(n.mpc.t0_param,copy(n.mpc.t0));
  end
  return nothing
end

"""
status=autonomousControl!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/1/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function autonomousControl!(n)
  mpcUpdate!(n)
 # update states based off of n.mpc.X0p
  # updateStates!(n)  TODO update this function so it does the following code
  for st in 1:n.numStates   # update states based off of n.mpc.X0p
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(n.r.x0_con[st,1], (n.mpc.X0p[st]+n.X0_tol[st]));
      JuMP.setRHS(n.r.x0_con[st,2],-(n.mpc.X0p[st]-n.X0_tol[st]));
    else
      JuMP.setRHS(n.r.x0_con[st],n.mpc.X0p[st]);
    end
  end
  optimize!(n)
  return n.r.status
end


end # module
