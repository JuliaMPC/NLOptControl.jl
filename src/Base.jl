module Base

using JuMP
using Ranges
using DataFrames

# These functions are required for MPC_Module.jl
export
  evalConstraints!,
  dvs2dfs,
  plant2dfs!,
  opt2dfs,
  postProcess!,
  optimize!,
  scale_tau

"""
scale_tau(τ,x₀,xₙ)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 12/25/2016 \n
--------------------------------------------------------------------------------------\n
"""
function scale_tau(τ,x₀,xₙ)
  (xₙ - x₀)/2*τ + (xₙ + x₀)/2;
end

"""
plant2dfs(n,r,s,u,sol)
# TODO: sometimes plant control models have different states and controls - > take this into account
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function plant2dfs!(n,r,s,u,sol)
  row, column=size(u)
  t_sample=Ranges.linspace(sol.t[1],sol.t[end],row);
  dfs_plant=DataFrame();
  dfs_plant[:t]=t_sample;

  for st in 1:n.numStates
    dfs_plant[n.state.name[st]]=[sol(t)[st] for t in t_sample];
  end

  for ctr in 1:n.numControls
    dfs_plant[n.control.name[ctr]]= u[ctr];
  end

  if s.reset
    r.dfs_plant=[dfs_plant];
  else
    push!(r.dfs_plant,dfs_plant);
  end
  nothing
end


"""
dvs2dfs(n,r)

# funtionality to save state and control data from optimization
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 2/10/2017 \n
--------------------------------------------------------------------------------------\n
"""
function dvs2dfs(n,r)
  dfs=DataFrame()
  dfs[:t]=r.t_st + n.mpc.t0;
  for i in 1:n.numStates
    dfs[n.state.name[i]]=r.X[:,i];
  end
  for i in 1:n.numControls
    if n.integrationMethod==:ts
      dfs[n.control.name[i]]=r.U[:,i];
    else
      dfs[n.control.name[i]]=[r.U[:,i];0];
    end
  end
  return dfs
end

"""
opt2dfs(r)

# funtionality to save optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 2/20/2017 \n
--------------------------------------------------------------------------------------\n
"""
function opt2dfs(r;kwargs...)

  kw = Dict(kwargs);
  # check to see if the user is initializing while compensating for control delay
  if !haskey(kw,:Init); kw_ = Dict(:Init => false); Init = get(kw_,:Init,0);
  else; Init = get(kw,:Init,0);
  end
  dfs_opt=DataFrame()

  if !Init
    dfs_opt[:t_solve]=r.t_solve
    dfs_opt[:obj_val]=r.obj_val
    dfs_opt[:status]=r.status
    dfs_opt[:iter_num]=r.iter_nums
  else
    dfs_opt[:t_solve]=0.0
    dfs_opt[:obj_val]=0.0
    dfs_opt[:status]=:Init
    dfs_opt[:iter_num]=0
  end

  return dfs_opt
end

"""
con2dfs(r)

# funtionality to save constraint data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/20/2017, Last Modified: 2/20/2017 \n
--------------------------------------------------------------------------------------\n
"""
function con2dfs(r)
  dfs_con=DataFrame()
  dfs_con[:con_val]=r.constraint.value;
  return dfs_con
end

"""
postProcess!(n,r,s)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/27/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function postProcess!(n,r,s;kwargs...)

  kw = Dict(kwargs);
  # check to see if the user is initializing while compensating for control delay
  if !haskey(kw,:Init); kw_ = Dict(:Init => false); Init = get(kw_,:Init,0);
  else; Init = get(kw,:Init,0);
  end

  if !Init
    if n.integrationMethod==:ps
      if n.finalTimeDV
        t = [scale_tau(n.ts[int],0.0,getvalue(n.tf)) for int in 1:n.Ni];     # scale time from [-1,1] to [t0,tf]
      else
        t = [scale_tau(n.ts[int],0.0,n.tf) for int in 1:n.Ni];
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

    if s.save
      if s.evalConstraints
        evalConstraints(n,r);
      end
      push!(r.dfs,dvs2dfs(n,r));
      push!(r.dfs_con,con2dfs(r));
      push!(r.dfs_opt,opt2dfs(r));
    end
  else  # no optimization run -> somtimes you drive straight to get started
    push!(r.dfs,nothing);
    push!(r.dfs_con,nothing);
    push!(r.dfs_opt,opt2dfs(r,;(:Init=>true)));
  end
  nothing
end


"""
status=optimize(mdl,n,r,s)

# solves JuMP model and saves optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/6/2017, Last Modified: 2/20/2017 \n
--------------------------------------------------------------------------------------\n
"""
function optimize!(mdl,n,r,s;Iter::Int64=0)
  t1 = time(); status=JuMP.solve(mdl); t2 = time();
  if s.save
    r.status=status;
    r.t_solve=t2-t1;
    r.obj_val=getobjectivevalue(mdl);
    r.iter_nums=Iter; # iteration number for a higher level algorithm
    r.eval_num=r.eval_num+1;  
    postProcess!(n,r,s);
  end
  #postProcess!(n,r,s);
  return status
end

"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 4/7/2017 \n
--------------------------------------------------------------------------------------\n
"""

function evalConstraints!(n,r)
  r.constraint.value=[];   # reset values
  r.constraint.nums=[]; s=1;
  for i = 1:length(r.constraint.handle)
    if r.constraint.name[i]==:dyn_con  # state constraits
      dfs=Vector{DataFrame}(n.numStates);
      con=DataFrame(step=1);
      l=0;
      for st in 1:n.numStates
        if n.integrationMethod==:ps
          temp=[getdual(r.constraint.handle[i][int][:,st]) for int in 1:n.Ni];
          vals=[idx for tempM in temp for idx=tempM];
          dfs[st]=DataFrame(step=1:sum(n.Nck);Dict(n.state.name[st] => vals)...);
          l=l+length(vals);
        else
          dfs[st]=DataFrame(step=1:n.N;Dict(n.state.name[st] => getdual(r.constraint.handle[i][:,st]))...);
          l=l+length(r.constraint.handle[i][:,st]);
        end
        if st==1; con=dfs[st]; else; con=join(con,dfs[st],on=:step); end
      end
    else
      S=JuMP.size(r.constraint.handle[i])
      if length(S)==1
        con = DataFrame(step=1:length(r.constraint.handle[i]);Dict(r.constraint.name[i] => getdual(r.constraint.handle[i][:]))...);
        l=S[1];
      elseif length(S)==2
        dfs=Vector{DataFrame}(S[1]);
        con=DataFrame(step=1);
        for idx in 1:S[1]
          dfs[idx] = DataFrame(step=1:S[2];Dict(r.constraint.name[i] => getdual(r.constraint.handle[i][idx,:]))...);
          if idx==1; con=dfs[idx]; else; con=join(con,dfs[idx],on=:step); end
        end
        l=S[1]*S[2];
      end
    end
    f=s+l-1;
    num=(i,r.constraint.name[i],@sprintf("length = %0.0f",l),string("indecies in g(x) = "),(s,f));
    push!(r.constraint.nums,num);
    push!(r.constraint.value,con);
    s=f+1;
  end
  nothing
end


end # module
