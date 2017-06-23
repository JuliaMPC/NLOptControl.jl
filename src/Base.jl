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
  scale_tau,
  try_import

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
plant2dfs!(n,sol)
# TODO: sometimes plant control models have different states and controls - > take this into account
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 6/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
function plant2dfs!(n,sol)
  row, column=size(n.r.u)
  t_sample=Ranges.linspace(sol.t[1],sol.t[end],row);
  dfs_plant=DataFrame();
  dfs_plant[:t]=t_sample;

  for st in 1:n.numStates
    dfs_plant[n.state.name[st]]=[sol(t)[st] for t in t_sample];
  end

  for ctr in 1:n.numControls
    dfs_plant[n.control.name[ctr]]= n.r.U[ctr];
  end

  if n.s.reset
    n.r.dfs_plant=[dfs_plant];
  else
    push!(n.r.dfs_plant,dfs_plant);
  end
  return nothing
end

"""
dvs2dfs(n)

# funtionality to save state and control data from optimization
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 5/28/2017 \n
--------------------------------------------------------------------------------------\n
"""
function dvs2dfs(n)
  dfs=DataFrame()
  dfs[:t]=n.r.t_st + n.mpc.t0;
  for i in 1:n.numStates
    dfs[n.state.name[i]]=n.r.X[:,i];
  end
  for i in 1:n.numControls
    if n.s.integrationMethod==:tm
      dfs[n.control.name[i]]=n.r.U[:,i];
    else
      dfs[n.control.name[i]]=[n.r.U[:,i];NaN];
    end
  end
  return dfs
end

"""
opt2dfs(n)
# funtionality to save optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 5/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function opt2dfs(n;kwargs...)
  kw = Dict(kwargs);
  # check to see if the user is initializing while compensating for control delay
  if !haskey(kw,:Init); kw_ = Dict(:Init => false); Init = get(kw_,:Init,0);
  else; Init = get(kw,:Init,0);
  end
  dfs_opt=DataFrame()

  if !Init
    dfs_opt[:t_solve]=n.r.t_solve
    dfs_opt[:obj_val]=n.r.obj_val
    dfs_opt[:status]=n.r.status
    dfs_opt[:iter_num]=n.r.iter_nums
  else
    dfs_opt[:t_solve]=0.0
    dfs_opt[:obj_val]=0.0
    dfs_opt[:status]=:Init
    dfs_opt[:iter_num]=0
  end
  return dfs_opt
end

"""
con2dfs(n)

# funtionality to save constraint data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/20/2017, Last Modified: 5/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function con2dfs(n)
  dfs_con=DataFrame()
  dfs_con[:con_val]=n.r.constraint.value;
  return dfs_con
end

"""
postProcess!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/27/2017, Last Modified: 5/28/2017 \n
--------------------------------------------------------------------------------------\n
"""
function postProcess!(n;kwargs...)
  kw = Dict(kwargs);
  # check to see if the user is initializing while compensating for control delay
  if !haskey(kw,:Init);Init=false;
  else; Init=get(kw,:Init,0);
  end

  if !Init
    if n.s.integrationMethod==:ps
      if n.s.finalTimeDV
        t=[scale_tau(n.ts[int],0.0,getvalue(n.tf)) for int in 1:n.Ni];     # scale time from [-1,1] to [t0,tf]
      else
        t=[scale_tau(n.ts[int],0.0,n.tf) for int in 1:n.Ni];
      end
      n.r.t_ctr=[idx for tempM in t for idx = tempM[1:end-1]];
      n.r.t_st=[n.r.t_ctr;t[end][end]];

    elseif n.s.integrationMethod==:tm
      if n.s.finalTimeDV
        n.r.t_ctr=append!([n.t0],cumsum(getvalue(n.dt)));
      else
        n.r.t_ctr=append!([n.t0],cumsum(n.dt));
      end
      n.r.t_st = n.r.t_ctr;
    end
    n.r.X=zeros(Float64,n.numStatePoints,n.numStates);
    n.r.U=zeros(Float64,n.numControlPoints,n.numControls);
    for st in 1:n.numStates
      n.r.X[:,st] = getvalue(n.r.x[:,st]);
    end
    for ctr in 1:n.numControls
      n.r.U[:,ctr] = getvalue(n.r.u[:,ctr]);
    end

    if n.s.evalConstraints
      evalConstraints!(n);
    end

    if n.s.save
      push!(n.r.dfs,dvs2dfs(n));
      push!(n.r.dfs_con,con2dfs(n));
      push!(n.r.dfs_opt,opt2dfs(n));
    end
  elseif n.s.save  # no optimization run -> somtimes you drive straight to get started
    push!(n.r.dfs,nothing);
    push!(n.r.dfs_con,nothing);
    push!(n.r.dfs_opt,opt2dfs(n,;(:Init=>true)));
  else
    warn("postProcess.jl did not do anything")
  end
  return nothing
end


"""
optimize!(n)

# solves JuMP model and saves optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/6/2017, Last Modified: 5/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function optimize!(n;Iter::Int64=0)
  t1=time(); status=JuMP.solve(n.mdl); t2=time();
  if n.s.save
    n.r.status=status;
    n.r.t_solve=t2-t1;
    n.r.obj_val=getobjectivevalue(n.mdl);
    n.r.iter_nums=Iter; # possible iteration number for a higher level algorithm
    n.r.eval_num=n.r.eval_num+1;
  end
  postProcess!(n);      # temporarily save data
  return nothing
end

"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 6/19/2017 \n
--------------------------------------------------------------------------------------\n
"""

function evalConstraints!(n)
  n.r.constraint.value=[];   # reset values
  n.r.constraint.nums=[]; s=1;
  for i = 1:length(n.r.constraint.handle)
    if n.r.constraint.name[i]==:dyn_con  # state constraits
      dfs=Vector{DataFrame}(n.numStates);
      con=DataFrame(step=1);
      l=0;
      for st in 1:n.numStates
        if n.s.integrationMethod==:ps
          temp=[getdual(n.r.constraint.handle[i][int][:,st]) for int in 1:n.Ni];
          vals=[idx for tempM in temp for idx=tempM];
          dfs[st]=DataFrame(step=1:sum(n.Nck);Dict(n.state.name[st] => vals)...);
          l=l+length(vals);
        else
          dfs[st]=DataFrame(step=1:n.N;Dict(n.state.name[st] => getdual(n.r.constraint.handle[i][:,st]))...);
          l=l+length(n.r.constraint.handle[i][:,st]);
        end
        if st==1; con=dfs[st]; else; con=join(con,dfs[st],on=:step); end
      end
    else
      S=0;
      try
        S=JuMP.size(n.r.constraint.handle[i])
      catch
        error("\n For now, the constraints cannot be in this form: \n
        con=@NLconstraint(mdl,n.r.u[1,1]==param); \n
        Write it in array form: \n
          con=@NLconstraint(mdl,[i=1],n.r.u[i,1]==param); \n")
      end
      if length(S)==1
        con = DataFrame(step=1:length(n.r.constraint.handle[i]);Dict(n.r.constraint.name[i] => getdual(n.r.constraint.handle[i][:]))...);
        l=S[1];
      elseif length(S)==2
        dfs=Vector{DataFrame}(S[1]);
        con=DataFrame(step=1);
        for idx in 1:S[1]
          dfs[idx] = DataFrame(step=1:S[2];Dict(n.r.constraint.name[i] => getdual(n.r.constraint.handle[i][idx,:]))...);
          if idx==1; con=dfs[idx]; else; con=join(con,dfs[idx],on=:step); end
        end
        l=S[1]*S[2];
      end
    end
    f=s+l-1;
    num=(i,n.r.constraint.name[i],@sprintf("length = %0.0f",l),string("indecies in g(x) = "),(s,f));
    push!(n.r.constraint.nums,num);
    push!(n.r.constraint.value,con);
    s=f+1;
  end
  return nothing
end

end # module
