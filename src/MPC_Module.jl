module MPC_Module

using JuMP

include("Base.jl")
using .Base

export
      initializeMPC,
      updateStates,
      updateX0,
      predictX0,
      mpcParams,
      mpcUpdate,
      simPlant,
      simModel,
      autonomousControl,
      MPC

type MPC
  # constants
  tp::Float64          # predication time (if finalTimeDV == true -> this is not known before optimization)
  tex::Float64         # execution horizon time
  max_iter::Int64      # maximum number of iterations

  # variables
  goal_reached::Bool           # flag  to indicate that goal has been reached
  t0::Float64                  # inital time
  tf::Float64                  # final time
  t0_param::Any                # parameter for t0
  X0p::Array{Float64,1}        # predicted initial states
  u::Array{Array{Float64,1},1} # all control signals

  # options
  PredictX0::Bool
end

function MPC()
  MPC(0.0,
      0.0,
      0,
      false,
      0.0,
      0.0,
      Any,
      Float64[],          # predicted intial state conditions
      Vector{Float64}[],
      true);
end

"""
initializeMPC(n,r,X0;(PredictX0=>true))
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function initializeMPC(n,r,X0;kwargs...)
  kw = Dict(kwargs);
  # check to see if the user would like to initialize the optimization X0 at a predicted future state
  if !haskey(kw,:Init); kw_=Dict(:PredictX0 => false); n.mpc.PredictX0=get(kw_,:PredictX0,0);
  else; n.mpc.PredictX0=get(kw,:PredictX0,0);
  end

  r.eval_num=1;                          # first evaluation number
  updateX0(n,r,X0;(:userUpdate=>true));
  mpcUpdate(n,r,r.eval_num);                        # initialize mpc parameters before running optimization
end

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
# for simulating the plant model given control commands
simPlant(pa,X0,t,U,t0,tf)

# TODO eventually the "plant" will be different from the "model"
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simPlant(n,r,s,pa,X0,t,U,t0,tf)
  sol=n.stateEquations(pa,X0,t,U,t0,tf);
  plant2dfs(n,r,s,U,sol);
  updateX0(n,r);
end

"""
# TODO eventually the "plant" will be different from the "model"
# TODO eventually have fixedTp be a model parameter if it is useful
predictX0(n,pa,r;(fixedTp=>true))
predictX0(n,pa,r;(fixedTp=>false))
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 4/7/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function predictX0(n,pa,r;kwargs...)
  kw = Dict(kwargs);

  # check to see if the user would like to make the prediction based off of previous solve time
  # this in case the user would like to update as quicly as posible instead of a fixed time
    # will make predictions worse, but may be more important to run quickly in particular applications
  if !haskey(kw,:Init); kw_=Dict(:fixedTp=> true); fixedTp=get(kw_,:fixedTp,0);
  else; fixedTp=get(kw,:fixedTp,0);
  end

 sol=simModel(n,pa,n.X0,r.t_ctr+n.mpc.t0,r.U,n.mpc.t0,n.mpc.tf)

  if fixedTp || r.eval_num==1
    n.mpc.X0p=sol(n.mpc.tf)[:];
  else
    if n.mpc.tf<n.mpc.t0+r.dfs_opt[r.eval_num][:t_solve][end] #TODO this should not be a vector!
      error("\nThe solution to the differential equation was not calculated to this late of a time. Consider decreasing `s.max_time` to a value less than n.mpc.tf. Or increase the amount of time that the model is simulated for\n")
    end
    n.mpc.X0p=sol(n.mpc.t0+r.dfs_opt[r.eval_num][:t_solve][end])[:]; # assume that it will take as long as the previos one could do the previous 5 or something...
  end
end


"""
updateX0(n,r)                           # uses solution from current plant data
updateX0(n,r,X0;(:userUpdate=>true))    # user defined update of X0
# updates intial states

# make sure that you call this before updateStates
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 3/10/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateX0(n,r,args...;kwargs...)
  kw = Dict(kwargs);
  # check to see how the intial states are being updated
  if !haskey(kw,:userUpdate); kw_ = Dict(:userUpdate => false); userUpdate = get(kw_,:userUpdate,0);
  else; userUpdate = get(kw,:userUpdate,0);
  end
  if userUpdate
    X0=args[1];
    if length(X0) != n.numStates
      error(string("\n Length of X0 must match number of states \n"));
    end
    n.X0=X0;
  else # update using current location of plant
    for st in 1:n.numStates
      n.X0[st]=r.dfs_plant[end][n.state.name[st]][end];
    end
  end
end

########################################################################################
# mpc functions
########################################################################################
"""
#TODO move this back to OCP? or make it more general for updating parameters -> move parameters out of this?
updateStates(n,r)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/17/2017, Last Modified: 2/21/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateStates(n,r)
  for st in 1:n.numStates   # update states
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(r.x0_con[st,1], (n.X0[st]+n.X0_tol[st]));
      JuMP.setRHS(r.x0_con[st,2],-(n.X0[st]-n.X0_tol[st]));
   else
     JuMP.setRHS(r.x0_con[st],n.X0[st]);
   end
  end
end

"""
# initialize mpc parameter settings
# this function is designed to be called once
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/6/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function mpcParams(n;tp::Float64=0.0,tex::Float64=0.0,max_iter::Int64=10)
  n.mpc::MPC = MPC(); # reset
  n.mpc.tp = tp;
  n.mpc.tex = tex;
  n.mpc.max_iter = max_iter;
end
mpcParams(n,c)=mpcParams(n;tp=c.m.tp,tex=c.m.tex,max_iter=c.m.mpc_max_iter);   # mpc setup

function mpcUpdate(n,r,ii;goal_reached::Bool=false)
  r.eval_num=ii;                         # update the evaluation number
  n.mpc.goal_reached = goal_reached;
  if !n.mpc.goal_reached
    n.mpc.t0 = n.mpc.tex*(r.eval_num-1);
    n.mpc.tf = n.mpc.tex*r.eval_num;
    setvalue(n.mpc.t0_param,n.mpc.t0);
  end
end

"""
status=autonomousControl(mdl,n,r,s,params);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/1/2017, Last Modified: 4/8/2017 \n
--------------------------------------------------------------------------------------\n
"""
function autonomousControl(mdl,n,r,s,params)
 for st in 1:n.numStates   # update states based off of n.mpc.X0p
  if any(!isnan(n.X0_tol[st]))
    JuMP.setRHS(r.x0_con[st,1], (n.mpc.X0p[st]+n.X0_tol[st]));
    JuMP.setRHS(r.x0_con[st,2],-(n.mpc.X0p[st]-n.X0_tol[st]));
  else
    JuMP.setRHS(r.x0_con[st],n.mpc.X0p[st]);
  end
  end
  optimize(mdl,n,r,s)
 end

end # module
