"""
r=postProcess(n,r)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/27/2017, Last Modified: 2/10/2017 \n
--------------------------------------------------------------------------------------\n
"""
function postProcess(n::NLOpt,r::Result)
  if n.integrationMethod==:ps
    if n.finalTimeDV
      t = [scale_tau(n.ts[int],n.t0,getvalue(n.tf)) for int in 1:n.Ni];     # scale time from [-1,1] to [t0,tf]
    else
      t = [scale_tau(n.ts[int],n.t0,n.tf) for int in 1:n.Ni];
    end
    r.t_ctr= [idx for tempM in t for idx = tempM[1:end-1]];
    r.t_st = [r.t_ctr;t[end][end]];
  elseif n.integrationMethod==:tm
    if n.finalTimeDV
      r.t_ctr =  append!([0.0],cumsum(getvalue(n.dt)));
    else
      r.t_ctr =  append!([0.0],cumsum(n.dt));
    end
    r.t_st = r.t_ctr;
  end
  r.X =zeros(Float64,n.numStatePoints,n.numStates);
  r.U =zeros(Float64,n.numControlPoints,n.numControls);
  for st in 1:n.numStates
    r.X[:,st] = getvalue(r.x[:,st]);
  end
  for ctr in 1:n.numControls
    r.U[:,ctr] = getvalue(r.u[:,ctr]);
  end

  # TODO make an MPC option -> i.e. there are certain things that do not need to be done
  evalConstraints(n,r);
  push!(r.dfs_opt,opt2dfs(r));
  push!(r.dfs,dvs2dfs(n,r));
  return r
end
