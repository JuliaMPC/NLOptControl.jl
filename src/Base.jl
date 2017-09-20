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
  try_import,
  intervals

"""
L = lagrange_basis_poly(x,x_data,Nc,j)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 1/25/2016
Citations: This function was influenced by the lagrange() function [located here](https://github.com/pjabardo/Jacobi.jl/blob/master/src/gauss_quad.jl)
--------------------------------------------------------------------------------------\n
# Input Arguments
* `x`: point to approximate function value at
* `x_data`: x data to used calculate basis polynomials
* `Nc`: order of Lagrange interpolating polynomial
* `j`: index of interest

# Output Arguments
* `L`: Lagrange basis polynomials

A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

"""
function lagrange_basis_poly(x,x_data,Nc,j)
    L = 1;
    for idx in 1:Nc+1 # use all of the data
      if idx!=j
        L = L*(x - x_data[idx])/(x_data[j]-x_data[idx]);
      end
    end
  return L
end

"""
y=interpolate_lagrange(ts[int],ts[int],stateMatrix[int][:,st])
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Created: 1/2/2017, Last Modified: 9/19/2017 \n
--------------------------------------------------------------------------------------\n
"""
function interpolate_lagrange{T<:Number}(x::AbstractArray{T},x_data,y_data)
  Nc = length(x_data) - 1

  if length(x_data)!=length(y_data)
      error(string("\n",
                    "-------------------------------------------------------", "\n",
                    "There is an error with the data vector lengths!!", "\n",
                    "-------------------------------------------------------", "\n",
                    "The following variables should be equal:", "\n",
                    "length(x_data) = ",length(x_data),"\n",
                    "length(y_data) = ",length(y_data),"\n"
                    )
            )
    end
    ns = length(x);
    L = zeros(Float64,Nc+1,ns);
    x = x[:]; x_data = x_data[:]; y_data = y_data[:]; # make sure data is in a column
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
      end
    end
    y = y_data'*L;
    return y
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified: 9/20/2017 \n
--------------------------------------------------------------------------------------\n
"""
function intervals(n,int,x,u)

  if typeof(x[1,1]) == JuMP.Variable
    # states
    x_int=Array{JuMP.Variable}(length(n.Nck_full[int]+1:n.Nck_full[int+1]),n.numStates);
    for st in 1:n.numStates # +1 adds the DV in the next interval
      x_int[:,st] = x[n.Nck_cum[int]+1:n.Nck_cum[int+1]+1,st];
    end

    # controls
    if int!=n.Ni
      u_int=Array{JuMP.Variable}(length(n.Nck_full[int]+1:n.Nck_full[int+1]),n.numControls);
    else                    # -1 -> removing control in last mesh interval
      u_int=Array{JuMP.Variable}(length(n.Nck_full[int]+1:n.Nck_full[int+1]-1),n.numControls);
    end
    for ctr in 1:n.numControls
      if int!=n.Ni          # +1 adds the DV in the next interval
        u_int[:,ctr] = u[n.Nck_cum[int]+1:n.Nck_cum[int+1]+1,ctr];
      else
        u_int[:,ctr] = u[n.Nck_cum[int]+1:n.Nck_cum[int+1],ctr];
      end
    end
  else
    # states
    x_int=Array{Any}(length(n.Nck_full[int]+1:n.Nck_full[int+1]),n.numStates);
    for st in 1:n.numStates # +1 adds the DV in the next interval
      x_int[:,st] = x[n.Nck_cum[int]+1:n.Nck_cum[int+1]+1,st];
    end

    # controls
    if int!=n.Ni
      u_int=Array{Any}(length(n.Nck_full[int]+1:n.Nck_full[int+1]),n.numControls);
    else                    # -1 -> removing control in last mesh interval
      u_int=Array{Any}(length(n.Nck_full[int]+1:n.Nck_full[int+1]-1),n.numControls);
    end
    for ctr in 1:n.numControls
      if int!=n.Ni          # +1 adds the DV in the next interval
        u_int[:,ctr] = u[n.Nck_cum[int]+1:n.Nck_cum[int+1]+1,ctr];
      else
        u_int[:,ctr] = u[n.Nck_cum[int]+1:n.Nck_cum[int+1],ctr];
      end
    end
  end

  return x_int,u_int
end

"""
# part of postProcess!()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Created: 9/19/2017, Last Modified: 9/19/2017 \n
--------------------------------------------------------------------------------------\n
"""
function interpolateLagrange(n; numPts::Int64=100)
  tf = getvalue(n.tf)
  numPts = 100
  # sample points
  n.r.t_polyPts = [linspace(tf/n.Ni*(int-1),tf/n.Ni*int,numPts) for int in 1:n.Ni]
  n.r.X_polyPts = [[zeros(numPts) for int in 1:n.Ni] for st in 1:n.numStates]
  n.r.U_polyPts = [[zeros(numPts) for int in 1:n.Ni] for ctr in 1:n.numControls]

  # time data points
  t_st_int = [n.r.t_st[n.Nck_cum[int]+1:n.Nck_cum[int+1]+1] for int in 1:n.Ni]
  t_ctr_int = [n.r.t_ctr[n.Nck_cum[int]+1:n.Nck_cum[int+1]+1] for int in 1:n.Ni-1]
  t_ctr_int = push!(t_ctr_int, n.r.t_st[n.Nck_cum[n.Ni]+1:n.Nck_cum[n.Ni+1]]) # -1 -> removing control in last mesh interval

  for int in 1:n.Ni

    # controls and states for this interval
    x_int, u_int = intervals(n, int, copy(n.r.X), n.r.U)

    # sample polynomial in interval at n.r.t_polyPts
    for st in 1:n.numStates
      n.r.X_polyPts[st][int] = interpolate_lagrange(n.r.t_polyPts[int], t_st_int[int], x_int[:,st])'
    end

    for ctr in 1:n.numControls
      n.r.U_polyPts[ctr][int] = interpolate_lagrange(n.r.t_polyPts[int], t_ctr_int[int], u_int[:,ctr])'
    end

  end

  nothing
end

"""
scale_tau(tau,ta,tb)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 9/18/2017 \n
--------------------------------------------------------------------------------------\n
"""
function scale_tau(tau,ta,tb)
  (tb - ta)/2*tau + (ta + tb)/2;
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
  if !haskey(kw,:Init); Init=false;
  else; Init=get(kw,:Init,0);
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

    if n.s.evalConstraints#&& n.r.status==:Optimal  # note may want to remove this &&
      evalConstraints!(n);
    end

    if n.s.save
      push!(n.r.dfs,dvs2dfs(n));
      push!(n.r.dfs_con,con2dfs(n));
      push!(n.r.dfs_opt,opt2dfs(n));
      interpolateLagrange(n)
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
    n.r.iter_nums=Iter;    # possible iteration number for a higher level algorithm
    n.r.eval_num=n.r.eval_num+1;
    postProcess!(n);      # temporarily save data
  end
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
          try
            dfs[idx] = DataFrame(step=1:S[2];Dict(n.r.constraint.name[i] => getdual(n.r.constraint.handle[i][idx,:]))...);
          catch
            dfs[idx] = DataFrame(step=1:S[2];Dict(n.r.constraint.name[i] => NaN)...); # fix for case where all of the states are not being constrainted, but some are within some XF_tol
          end
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
