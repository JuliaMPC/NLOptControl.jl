module MPC_Module

using JuMP

export
      updateStates,
      updateX0,
      mpcParams,
      mpcUpdate,
      simPlant,
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
      Vector{Float64}[]);
end

########################################################################################
# mpc functions
########################################################################################
"""
#TODO move this back to OCP? or make it more general for updating parameters -> move parameters out of this?
updateStates(n,r,params,X0,SA,UX)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/17/2017, Last Modified: 2/21/2017 \n
--------------------------------------------------------------------------------------\n
"""
function updateStates(n,r,params,X0,SA,UX)
  setvalue(params[3], SA)   # update desired steering angle
  setvalue(params[2], UX)   # update speed
  for st in 1:n.numStates   # update states
    if any(!isnan(n.X0_tol[st]))
      JuMP.setRHS(r.x0_con[st,1], (X0[st]+n.X0_tol[st]));
      JuMP.setRHS(r.x0_con[st,2],-(X0[st]-n.X0_tol[st]));
   else
     JuMP.setRHS(r.x0_con[st],X0[st]);
    end
  end
end

"""
updateX0(n,r)                           # uses solution from current plant data
updateX0(n,r,X0;(:userUpdate=>true))    # user defined update of X0
# updates intial states
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
mpcParams(n,c)=mpcParams(n;tp=c.m.tp,tex=c.m.tex,max_iter=c.m.max_iter);   # mpc setup

function mpcUpdate(n,r;goal_reached::Bool=false)
  n.mpc.goal_reached = goal_reached;
  if !n.mpc.goal_reached
    n.mpc.t0 = n.mpc.tex*(r.eval_num-1);
    n.mpc.tf = n.mpc.tex*r.eval_num;
    setvalue(n.mpc.t0_param,n.mpc.t0);
  end
end

"""
r=simPlant(n,r)
# for simulating the plant model given commands
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 2/14/2017 \n
--------------------------------------------------------------------------------------\n
"""
function simPlant(n,r;plantModel::Function=[])
  #push!(r.dfs,dvs2dfs(n,r));
#TODO finish
end

end # module
