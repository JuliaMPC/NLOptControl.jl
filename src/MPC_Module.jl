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
      MPC

type MPC
  # constants
  tp::Float64          # predication time (if finalTimeDV == true -> this is not known before optimization)
  tex::Float64         # execution horizon time
  max_iter::Int64      # maximum number of iterations

  # variables
  goal_reached::Bool           # flag  to indicate that goal has been reached
  t0_actual                    # actual intial time
  t0::Float64                  # mpc inital time
  tf::Float64                  # mpc final time
  t0_param::Any                # parameter for mpc t0
  X0p::Array{Float64,1}        # predicted initial states
  X0                          # array of all actual intial states (n.X0 leads because we are in simulation)

  # options
  PredictX0::Bool
  FixedTp::Bool
end

function MPC()
  MPC(0.0,
      0.0,
      0,
      false,
      0.0,
      0.0,
      0.0,
      Any,
      Float64[],          # predicted intial state conditions
      [Float64[]],
      true,
      true);
end

########################################################################################
# mpc functions
########################################################################################
"""
initializeMPC(n,r)

initializeMPC!(n,r;FixedTp=c.m.FixedTp,PredictX0=c.m.PredictX0,tp=c.m.tp,tex=copy(c.m.tex),max_iter=c.m.mpc_max_iter);

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function initializeMPC!(n,r;FixedTp::Bool=true,PredictX0::Bool=true,tp::Float64=5.0,tex::Float64=0.5,max_iter::Int64=10)
  n.mpc::MPC=MPC();
  n.mpc.FixedTp=FixedTp;
  n.mpc.PredictX0=PredictX0;
  n.mpc.tp=tp;
  n.mpc.tex=tex;
  n.mpc.max_iter=max_iter;
  n.mpc.t0=0.0;
  n.mpc.t0_actual=0.0;
  n.mpc.tf=tex;
  nothing
end
#initializeMPC!(n,r,c)=initializeMPC!(n,r;FixedTp=c.m.FixedTp,PredictX0=c.m.PredictX0,tp=c.m.tp,tex=c.m.tex,max_iter=c.m.mpc_max_iter);

"""
# for simulating the model of the plant given control commands
simModel(pa,X0,t,U,t0,tf)

# TODO eventually the "plant" will be different from the "model"
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simModel(n,pa,X0,t,U,t0,tf)
  n.stateEquations(pa,X0,t,U,t0,tf);
end

"""
# NOTE: this function is called inside simPlant()
updateX0(n,r,X0;(:userUpdate=>true))    # user defined update of X0
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 4/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateX0!(n,r,args...;kwargs...)
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
      n.X0[st]=r.dfs_plant[end][n.state.name[st]][end];
    end
  end
  updateStates!(n,r)
  append!(n.mpc.X0,[copy(n.X0)])
  nothing
end

"""
# for simulating the plant model given control commands
simPlant(pa,X0,t,U,t0,tf)

# TODO eventually the "plant" will be different from the "model"
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simPlant!(n,r,s,pa,X0,t,U,t0,tf)
  sol=n.stateEquations(pa,X0,t,U,t0,tf);
  plant2dfs!(n,r,s,U,sol);
  t0p=updateX0!(n,r);
  nothing
end

"""
# TODO eventually the "plant" will be different from the "model"
t0p=predictX0(n,pa,r)
# also shifts the intial time, but this is accounted for in updateMPC()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function predictX0!(n,pa,r)

  if n.mpc.tf==n.mpc.t0
    error("n.mpc.tf==n.mpc.t0")
  end

  if n.mpc.FixedTp || r.eval_num==1
    t0p=n.mpc.tex;
  else
    t0p=r.dfs_opt[end][:t_solve][1]
  end
  # based off of "current X0". Even though we may have the next X0 we should not (i.e.look at @show length(n.mpc.X0)). It is because it is a simulation (in reality they would be running in parallel)
  sol=simModel(n,pa,n.mpc.X0[r.eval_num],r.t_ctr+n.mpc.t0,r.U,n.mpc.t0,n.mpc.t0+t0p)
  n.mpc.X0p=sol(n.mpc.t0+t0p)[:];

  return t0p
end
"""
# this function currently only works when the "actual model"==simPlant()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/9/2017, Last Modified: 4/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function driveStraight!(n,pa,r,s;t0::Float64=n.mpc.t0,tf::Float64=n.mpc.tf)

  # add these signals to r so that they can be used for predictions during optimization
  r.U=0*Matrix{Float64}(n.numControlPoints,n.numControls);
  r.t_ctr=Vector(Ranges.linspace(t0,tf,n.numControlPoints)); # gives a bunch of points
  postProcess!(n,r,s;(:Init=>true)); r.eval_num=1;          # to make solutions line up

  # simulate the "actual vehicle" response
  simPlant!(n,r,s,pa,n.X0,r.t_ctr+n.mpc.t0,r.U,t0,tf)
  nothing
end

"""
updateStates(n,r)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/17/2017, Last Modified: 2/21/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateStates!(n,r)
  for st in 1:n.numStates   # update states
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(r.x0_con[st,1], (n.X0[st]+n.X0_tol[st]));
      JuMP.setRHS(r.x0_con[st,2],-(n.X0[st]-n.X0_tol[st]));
   else
     JuMP.setRHS(r.x0_con[st],n.X0[st]);
   end
  end
  nothing
end

"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 4/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function mpcUpdate!(n,pa,r)
  if n.mpc.PredictX0        # predict where X0 will be when optimized signal is actually sent to the vehicle
    t0p=predictX0!(n,pa,r); # predicted start time -> important for time varying constraints
  else
    t0p=0; n.mpc.X0p=n.mpc.X0[r.eval_num];  # current "known" plant states
  end
  n.mpc.t0=copy(n.mpc.t0_actual+t0p);     # there are two different time scales-> the plant leads by t0p
  n.mpc.tf=copy(n.mpc.t0+n.mpc.tex);
  if n.mpc.t0_param!=Any
    setvalue(n.mpc.t0_param,copy(n.mpc.t0));
  end
  nothing
end

"""
status=autonomousControl(mdl,n,r,s,pa)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/1/2017, Last Modified: 4/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function autonomousControl!(mdl,n,r,s,pa)
  mpcUpdate!(n,pa,r)

  for st in 1:n.numStates   # update states based off of n.mpc.X0p
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(r.x0_con[st,1], (n.mpc.X0p[st]+n.X0_tol[st]));
      JuMP.setRHS(r.x0_con[st,2],-(n.mpc.X0p[st]-n.X0_tol[st]));
    else
      JuMP.setRHS(r.x0_con[st],n.mpc.X0p[st]);
    end
  end

  return optimize!(mdl,n,r,s)
 end

end # module
