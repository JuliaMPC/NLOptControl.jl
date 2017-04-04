
"""
defineSolver(n; solver=:KNITRO)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 2/9/2017 \n
-------------------------------------------------------------------------------------\n
"""
function defineSolver(n::NLOpt;
                      solver::Symbol=:IPOPT)
  n.solver=solver;
end

"""
linearStateTolerances(n)
# the purpose of this function is to taper the tolerances on the constant state constraints
# the idea is that when doing MPC, the final states are well within the bounds so that the next optimization is not initalized at an infeasible point
# if you want a constant bond, set the slope to zero
# default is a positive slope on the lower bound and a negative slope on the upper bound
# this functionality in not needed for states like position, so you do not need to add a linearStateTol for all states

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/23/2017, Last Modified: 3/25/2017 \n
-------------------------------------------------------------------------------------\n
"""
function linearStateTolerances(n::NLOpt;
                        mXL::Array{Any,1}=falses(n.numStates),
                        mXU::Array{Any,1}=falses(n.numStates))

n.mXL=mXL;  n.mXU=mXU;
for st in 1:n.numStates
  # lower state constraint
  if n.mXL[st]!=false
    if !isnan(n.XL[st])
      for j in 1:n.numStatePoints
      n.XL_var[st,j]=n.XL[st] + n.mXL[st]*(j/n.numStatePoints); # lower
      end
    end
  end

  # upper state constraint
  if n.XU[st]!=false
    if !isnan(n.XU[st])
      for j in 1:n.numStatePoints
      n.XU_var[st,j]=n.XU[st] + n.mXU[st]*(j/n.numStatePoints); # upper
      end
    end
  end
end

end

"""
defineTolerances(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/8/2017, Last Modified: 2/8/2017 \n
-------------------------------------------------------------------------------------\n
"""
function defineTolerances(n::NLOpt;
                          X0_tol::Array{Float64,1}=0.05*ones(Float64,n.numStates,),
                          XF_tol::Array{Float64,1}=0.05*ones(Float64,n.numStates,))
  n.X0_tol=X0_tol; n.XF_tol=XF_tol;
end

"""
n = create_tV(mdl,n)
# define a time vector (n.tV) for use with time varying constraints when (finalTimeDV=>true)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/17/2017, Last Modified: 3/17/2017 \n
--------------------------------------------------------------------------------------\n
"""
function create_tV(mdl::JuMP.Model,n::NLOpt)

  if n.integrationMethod==:ps
    if n.finalTimeDV
      # create mesh points, interval size = tf_var/Ni
      tm = @NLexpression(mdl, [idx=1:n.Ni+1], (idx-1)*n.tf/n.Ni);
      # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
      ts = [Array(Any,n.Nck[int]+1,) for int in 1:n.Ni];
      for int in 1:n.Ni
        ts[int][1:end-1]=@NLexpression(mdl,[j=1:n.Nck[int]], (tm[int+1]-tm[int])/2*n.τ[int][j] +  (tm[int+1]+tm[int])/2);
        ts[int][end]=@NLexpression(mdl, n.tf/n.Ni*int) # append +1 at end of each interval
      end

      tt1 = [idx for tempM in ts for idx = tempM[1:end-1]];
      tmp = [tt1;ts[end][end]];
      n.tV = @NLexpression(mdl,[j=1:n.numStatePoints], n.t0 + tmp[j]);
    else
      error("finish this")
    end
  else
    error("finish this")

    if n.finalTimeDV
      # vector with the design variable in it
      t = Array(Any,n.N+1,1);
      tmp = [0;cumsum(n.dt)];
      n.tV = @NLexpression(mdl,[j=1:n.numStatePoints], n.t0 + tmp[j]);
    else

    end
  end

  return n
end

"""
# integrating JuMP variables
Expr = integrate(mdl,n,u;(:mode=>:control))
Expr = integrate(mdl,n,u,idx=1;C=0.5,(:variable=>:control),(:integrand=>:squared))
Expr = integrate(mdl,n,r.u[:,1];D=1.2,(:variable=>:control),(:integrand=>:squared),(:integrandAlgebra=>:subtract))

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/2/2017, Last Modified: 1/27/2017 \n
--------------------------------------------------------------------------------------\n
"""
function integrate(mdl::JuMP.Model,n::NLOpt,V::Array{JuMP.Variable,1}, args...; C::Float64=1.0,D=0.0,kwargs...)
  kw = Dict(kwargs);
  if !haskey(kw,:integrand); kw_ = Dict(:integrand => :default); integrand = get(kw_,:integrand,0);
  else; integrand = get(kw,:integrand,0);
  end
  if !haskey(kw,:integrandAlgebra); kw_ = Dict(:integrandAlgebra => :default); integrandAlgebra = get(kw_,:integrandAlgebra,0);
  else; integrandAlgebra = get(kw,:integrandAlgebra,0);
  end

  if integrandAlgebra ==:subtract
    if n.integrationMethod==:tm
      if integrand == :default      # integrate V
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl,C*sum((V[j+1]-D)*n.tf/(n.N) for j = 1:n.N));  #TODO fix this.. there is an extra dv here for control, but it does not effect solution
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl,C*sum(0.5*((V[j]-D+V[j+1]-D)*n.tf/(n.N) for j = 1:n.N)));
        end
      elseif integrand == :squared # integrate V^2
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl, C*sum(((V[j+1]-D)^2)*n.tf/(n.N) for j = 1:n.N));
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl, C*sum(0.5*((V[j]-D)^2+(V[j+1]-D)^2)*n.tf/(n.N) for j = 1:n.N));
        end
      else
        error("\n Check :integrand \n")
      end
    elseif n.integrationMethod==:ps
      if !haskey(kw,:mode); kw_ = Dict(:mode => :quadrature); mode = get(kw_,:mode,0);
      else; mode  = get(kw,:mode,0);
      end
      variable = get(kw,:variable,0);
      if variable == :state; Nck_cum  = [0;cumsum(n.Nck+1)];
      elseif variable == :control; Nck_cum = [0;cumsum(n.Nck)];
      else; error("\n Set the variable to either (:variable => :state) or (:variable => :control). \n")
      end

      if mode == :quadrature  #TODO recalculate ws based off of time
        if integrand == :default      # integrate V
          @NLexpression(mdl, temp[int=1:n.Ni], ((n.tf-n.t0)/2)*sum((n.ωₛ[int])[j]*((V[Nck_cum[int]+1:Nck_cum[int+1]])[j]-D) for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, C*sum(temp[int] for int = 1:n.Ni));
        elseif integrand == :squared # integrate V^2
          @NLexpression(mdl, temp[int=1:n.Ni],((n.tf-n.t0)/2)*C*sum((n.ωₛ[int])[j]*((V[Nck_cum[int]+1:Nck_cum[int+1]])[j]-D)*((V[Nck_cum[int] + 1:Nck_cum[int+1]])[j]-D) for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, sum(temp[int] for int = 1:n.Ni));
        else
          error("\n Check :integrand \n")
        end
      elseif mode == :LGRIM# TODO add in option to allow for integration using IMatrix
          error("\n Not implemented yet!! \n")
      end
    end
  elseif integrandAlgebra==:default
    if n.integrationMethod==:tm
      if integrand == :default      # integrate V
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl,C*sum(V[j+1]*n.tf/(n.N) for j = 1:n.N));  #TODO fix this.. there is an extra dv here for control, but it does not effect solution
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl,C*sum(0.5*(V[j]+V[j+1])*n.tf/(n.N) for j = 1:n.N));
        end
      elseif integrand == :squared # integrate V^2
        if n.integrationScheme==:bkwEuler
          Expr =  @NLexpression(mdl, C*sum((V[j+1]^2)*n.tf/(n.N) for j = 1:n.N));
        elseif n.integrationScheme==:trapezoidal
          Expr =  @NLexpression(mdl, C*sum(0.5*(V[j]^2+V[j+1]^2)*n.tf/(n.N) for j = 1:n.N));
        end
      else
        error("\n Check :integrand \n")
      end
    elseif n.integrationMethod==:ps
      if !haskey(kw,:mode); kw_ = Dict(:mode => :quadrature); mode = get(kw_,:mode,0);
      else; mode  = get(kw,:mode,0);
      end
      variable = get(kw,:variable,0);
      if variable == :state; Nck_cum  = [0;cumsum(n.Nck+1)];
      elseif variable == :control; Nck_cum = [0;cumsum(n.Nck)];
      else; error("\n Set the variable to either (:variable => :state) or (:variable => :control). \n")
      end

      if mode == :quadrature  #TODO recalculate ws based off of time
        if integrand == :default      # integrate V
          @NLexpression(mdl, temp[int=1:n.Ni], ((n.tf-n.t0)/2)*sum((n.ωₛ[int])[j]*(V[Nck_cum[int]+1:Nck_cum[int+1]])[j] for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, C*sum(temp[int] for int = 1:n.Ni));
        elseif integrand == :squared # integrate V^2
          @NLexpression(mdl, temp[int=1:n.Ni],((n.tf-n.t0)/2)*C*sum((n.ωₛ[int])[j]*(V[Nck_cum[int]+1:Nck_cum[int+1]])[j]*(V[Nck_cum[int] + 1:Nck_cum[int+1]])[j] for j = 1:n.Nck[int]));
          Expr =  @NLexpression(mdl, sum(temp[int] for int = 1:n.Ni));
        else
          error("\n Check :integrand \n")
        end
      elseif mode == :LGRIM# TODO add in option to allow for integration using IMatrix
          error("\n Not implemented yet!! \n")
      end
    end
  end
  return Expr
end

"""
mdl=build(n)
# build JuMP model
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 2/22/2017 \n
--------------------------------------------------------------------------------------\n
"""
#TODO allow user to pass solver options
function build(n::NLOpt)
  if n.solver==:IPOPT
    mdl=Model(solver=IpoptSolver());
  else n.solver==:KNITRO
    mdl=Model(solver=KnitroSolver());
  end
  return mdl
end

"""
optimize(mdl,n,r,s)

# solves JuMP model and saves optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/6/2017, Last Modified: 2/20/2017 \n
--------------------------------------------------------------------------------------\n
"""
function optimize(mdl::JuMP.Model,n::NLOpt,r::Result,s::Settings;Iter::Int64=0)
  t1 = time(); status=JuMP.solve(mdl); t2 = time();
  if s.save
    push!(r.status,status);
    push!(r.t_solve,(t2 - t1));
    push!(r.obj_val, getobjectivevalue(mdl));
    push!(r.iter_nums,Iter); # iteration number for a higher level algorithm
    r.eval_num=length(r.status);
  end
  postProcess(n,r,s);
  return status
end

########################################################################################
# constraint data functions
########################################################################################
"""
# funtionality to save constraint data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function newConstraint(r::Result,handle,name::Symbol)
  initConstraint(r)
  r.constraint::Constraint = r.constraint
  push!(r.constraint.handle,handle)
  push!(r.constraint.name,name)
end

function initConstraint(r::Result)
  if r.constraint == nothing
    r.constraint = Constraint()
  end
end

function evalConstraints(n::NLOpt, r::Result)
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
end

"""
maxDualInf = evalMaxDualInf(n,r)
# funtionality to evaluate maximum dual infeasibility of problem
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/13/2017, Last Modified: 2/13/2017 \n
--------------------------------------------------------------------------------------\n
"""
function evalMaxDualInf(n::NLOpt, r::Result)
  num=length(r.constraint.handle); dual_con_temp=zeros(num);
  for i = 1:length(r.constraint.handle)
    if !isempty(r.constraint.handle[i])
      if r.constraint.name[i]==:dyn_con  # state constraits
        temp1=zeros(n.numStates);
        for st in 1:n.numStates
          if n.integrationMethod==:ps
            temp=[getdual(r.constraint.handle[i][int][:,st]) for int in 1:n.Ni];
            vals=[idx for tempM in temp for idx=tempM];
            temp1[st]=maximum(vals);
          else
            temp1[st] = maximum(getdual(r.constraint.handle[i][:,st]));
          end
        end
        dual_con_temp[i]=maximum(temp1);
      else
        S=JuMP.size(r.constraint.handle[i])
        if length(S)==1
          dual_con_temp[i]=maximum(getdual(r.constraint.handle[i][:]));
        elseif length(S)==2
          temp1=zeros(S[1]);
          for idx in 1:S[1]
            temp1[idx] = maximum(getdual(r.constraint.handle[i][idx,:]));
          end
          dual_con_temp[i]=maximum(temp1);
        end
      end
    end
  end
  return maximum(maximum(maximum(dual_con_temp)))
end

########################################################################################
# state data functions
########################################################################################
"""
# initialize state names
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""

function initStateNames(n::NLOpt)
  State([Symbol("x$i") for i in 1:n.numStates],
        [String("x$i") for i in 1:n.numStates],
       );
end

function stateNames(n::NLOpt,names,descriptions)
  if !n.define
    error("\n call define() before calling controlNames() \n")
  end
  if length(names)!=n.numStates || length(descriptions)!=n.numStates
    error("\n Check sizes of names and descriptions \n")
  end
  n.state::State = State() # reset
  for i in 1:n.numStates
    push!(n.state.name,names[i])
    push!(n.state.description,descriptions[i])
  end
end

########################################################################################
# control data functions
########################################################################################
"""
# initialize control names
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function initControlNames(n::NLOpt)
  Control([Symbol("u$i") for i in 1:n.numControls],
          [String("u$i") for i in 1:n.numControls]);
end

function controlNames(n::NLOpt,names,descriptions)
  if !n.define
    error("\n call define() before calling controlNames() \n")
  end
  if length(names)!=n.numControls || length(descriptions)!=n.numControls
    error("\n Check sizes of names and descriptions \n")
  end
  n.control::Control = Control() # reset
  for i in 1:n.numControls
    push!(n.control.name,names[i])
    push!(n.control.description,descriptions[i])
  end
end

########################################################################################
# data related functions
########################################################################################

"""
resultsDir(results_dir)
# removes results folder and creates a new one
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 2/17/2017 \n
--------------------------------------------------------------------------------------\n
"""
function resultsDir(results_dir::String)
	if isdir(results_dir)
		rm(results_dir; recursive=true)
		print("\n The old results have all been deleted! \n \n")
	end
	mkdir(results_dir)
end


"""
dvs2dfs(n,r)

# funtionality to save state and control data from optimization
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 2/10/2017 \n
--------------------------------------------------------------------------------------\n
"""
function dvs2dfs(n::NLOpt,r::Result)
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
plant2dfs(n,r,s,u,sol)

# funtionality to save state and control data from plant model

# TODO: sometimes plant control models have different states and controls - > take this into account
# TODO: allow for more general input than output from DifferentialEquations.jl
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/14/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function plant2dfs(n::NLOpt,r::Result,s::Settings,u,sol)

  t_sample=linspace(sol.t[1],sol.t[end],length(u[1]));
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
end

"""
opt2dfs(r)

# funtionality to save optimization data
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/10/2017, Last Modified: 2/20/2017 \n
--------------------------------------------------------------------------------------\n
"""
function opt2dfs(r::Result;kwargs...)

  kw = Dict(kwargs);
  # check to see if the user is initializing while compensating for control delay
  if !haskey(kw,:Init); kw_ = Dict(:Init => false); Init = get(kw_,:Init,0);
  else; Init = get(kw,:Init,0);
  end

  dfs_opt=DataFrame()

  if !Init
    id = find(r.t_solve)
    idx = id[end]
    dfs_opt[:t_solve]=r.t_solve[1:idx]
    dfs_opt[:obj_val]=r.obj_val[1:idx]
    dfs_opt[:status]=r.status[1:idx]
    dfs_opt[:iter_num]=r.iter_nums[idx];
  else
    # find a better spot for this TODO consider eliminating altogether
    push!(r.t_solve,0.0);
    push!(r.obj_val,0.0);
    push!(r.status,:Init);
    push!(r.iter_nums,0);

    dfs_opt[:t_solve]=r.t_solve[1]
    dfs_opt[:obj_val]=r.obj_val[1]
    dfs_opt[:status]=r.status[1]
    dfs_opt[:iter_num]=r.iter_nums[1];
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
function con2dfs(r::Result)
  dfs_con=DataFrame()
  dfs_con[:con_val]=r.constraint.value;
  return dfs_con
end


"""
# maximum(x->maximum(x[:A]), dfs) -> consider
# maximum(x->maximum(filter(y->y !== nothing, x[:A])), dfs)
# minimum(x->maximum(filter(y->y !== nothing, x[:t])), dfs)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/8/2016, Last Modified: 3/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function minDF(dfs,varb)
  k=length(dfs); tmp=1e10*ones(k,1);
  for i in 1:k
    if dfs[i]!=nothing
      tmp[i]=minimum(dfs[i][varb])
    end
  end
  minimum(tmp)  # find the minimum
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/8/2016, Last Modified: 3/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function maxDF(dfs,varb)
  k=length(dfs); tmp=1e-10*ones(k,1);
  for i in 1:k
    if dfs[i]!=nothing
      tmp[i]=maximum(dfs[i][varb])
    end
  end
  maximum(tmp)  # find the minimum
end
