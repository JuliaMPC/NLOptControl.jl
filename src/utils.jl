"""
linearStateTolerances!(n)
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
function linearStateTolerances!(n::NLOpt;
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
  return nothing
end

"""
defineTolerances!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/8/2017, Last Modified: 2/8/2017 \n
-------------------------------------------------------------------------------------\n
"""
function defineTolerances!(n::NLOpt;
                          X0_tol::Array{Float64,1}=0.05*ones(Float64,n.numStates,),
                          XF_tol::Array{Float64,1}=0.05*ones(Float64,n.numStates,))
  n.X0_tol=X0_tol; n.XF_tol=XF_tol;
  return nothing
end

"""
create_tV!(n)
# define a time vector (n.tV) for use with time varying constraints when (finalTimeDV=>true)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/17/2017, Last Modified: 3/17/2017 \n
--------------------------------------------------------------------------------------\n
"""
function create_tV!(n::NLOpt)

  if n.s.integrationMethod==:ps
    if n.s.finalTimeDV
      # create mesh points, interval size = tf_var/Ni
      tm = @NLexpression(n.mdl, [idx=1:n.Ni+1], (idx-1)*n.tf/n.Ni);
      # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
      ts = [Array{Any}(n.Nck[int]+1,) for int in 1:n.Ni];
      for int in 1:n.Ni
        ts[int][1:end-1]=@NLexpression(n.mdl,[j=1:n.Nck[int]], (tm[int+1]-tm[int])/2*n.tau[int][j] +  (tm[int+1]+tm[int])/2);
        ts[int][end]=@NLexpression(n.mdl, n.tf/n.Ni*int) # append +1 at end of each interval
      end
      tt1 = [idx for tempM in ts for idx = tempM[1:end-1]];
      tmp = [tt1;ts[end][end]];
      n.tV = @NLexpression(n.mdl,[j=1:n.numStatePoints], n.mpc.t0_param + tmp[j]);
    else
      error("finish this")
    end
  else
    error("finish this")

    if n.s.finalTimeDV
      # vector with the design variable in it
      t = Array{Any}(n.N+1,1);
      tmp = [0;cumsum(n.dt)];
      n.tV = @NLexpression(n.mdl,[j=1:n.numStatePoints], n.t0 + tmp[j]);
    else

    end
  end
  return nothing
end

"""
obj=integrate!(n,:(u1))
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/2/2017, Last Modified: 9/18/2017 \n
--------------------------------------------------------------------------------------\n
"""
function integrate!(n::NLOpt,V::Expr)
  if n.s.integrationMethod==:ps
    integral_expr=[Array{Any}(n.Nck[int]) for int in 1:n.Ni];
    for int in 1:n.Ni
      x_int,u_int=intervals(n,int,n.r.x,n.r.u);
      L=size(x_int)[1]-1;
      integral_expr[int][:]=NLExpr(n,V,x_int,u_int,L)
    end
    @NLexpression(n.mdl, temp[int=1:n.Ni], (n.tf-n.t0)/2*sum(n.ws[int][j]*integral_expr[int][j] for j = 1:n.Nck[int]) );
    expression=@NLexpression(n.mdl, sum(temp[int] for int = 1:n.Ni));
  elseif n.s.integrationMethod==:tm
    L=size(n.r.x)[1];
    temp=NLExpr(n,V,n.r.x,n.r.u,L);
    if n.s.integrationScheme==:bkwEuler
      expression=@NLexpression(n.mdl,sum(temp[j+1]*n.tf/n.N for j = 1:n.N) );
    elseif n.s.integrationScheme==:trapezoidal
      expression=@NLexpression(n.mdl,sum(0.5*(temp[j]+temp[j+1])*n.tf/n.N for j = 1:n.N) );
    else
      error("\n $(n.s.integrationScheme) not defined in integrationSchemes\n")
    end
  else
    error("\n $(n.s.integrationMethod) not defined in integrationMethods \n")
  end
  return expression
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
function initConstraint!(n)
  if n.r.constraint==nothing
    n.r.constraint=Constraint()
  end
  return nothing
end

function newConstraint!(n,handle,name::Symbol)
  initConstraint!(n)
  n.r.constraint::Constraint=n.r.constraint
  push!(n.r.constraint.handle,handle)
  push!(n.r.constraint.name,name)
  return nothing
end


"""
maxDualInf = evalMaxDualInf(n)
# funtionality to evaluate maximum dual infeasibility of problem
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/13/2017, Last Modified: 2/13/2017 \n
--------------------------------------------------------------------------------------\n
"""
function evalMaxDualInf(n::NLOpt)
  num=length(n.r.constraint.handle); dual_con_temp=zeros(num);
  for i = 1:length(n.r.constraint.handle)
    if !isempty(n.r.constraint.handle[i])
      if n.r.constraint.name[i]==:dyn_con  # state constraits
        temp1=zeros(n.numStates);
        for st in 1:n.numStates
          if n.s.integrationMethod==:ps
            temp=[getdual(n.r.constraint.handle[i][int][:,st]) for int in 1:n.Ni];
            vals=[idx for tempM in temp for idx=tempM];
            temp1[st]=maximum(vals);
          else
            temp1[st] = maximum(getdual(n.r.constraint.handle[i][:,st]));
          end
        end
        dual_con_temp[i]=maximum(temp1);
      else
        S=JuMP.size(n.r.constraint.handle[i])
        if length(S)==1
          dual_con_temp[i]=maximum(getdual(n.r.constraint.handle[i][:]));
        elseif length(S)==2
          temp1=zeros(S[1]);
          for idx in 1:S[1]
            temp1[idx] = maximum(getdual(n.r.constraint.handle[i][idx,:]));
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
        [String("x$i") for i in 1:n.numStates]);
end

function states!(n::NLOpt,names;descriptions=[])
  if !n.define
    error("\n call define() before calling stateNames() \n")
  end
  if length(names)!=n.numStates
    error("\n Check size of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.numStates
    error("\n Check size of descriptions \n")
  end
  n.state::State = State() # reset
  for i in 1:n.numStates
    if names[i]==:xxx
      error("xxx is OFF limits for a state name; please choose something else. \n")
    end
    push!(n.state.name,names[i])
    if !isempty(descriptions)
        push!(n.state.description,descriptions[i])
    end
  end
  return nothing
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

function controls!(n::NLOpt,names;descriptions=[])
  if !n.define
    error("\n call define() before calling controlNames() \n")
  end
  if length(names)!=n.numControls
    error("\n Check sizes of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.numControls
    error("\n Check size of descriptions \n")
  end

  n.control::Control = Control() # reset
  for i in 1:n.numControls
    if names[i]==:uuu
      error("uuu is OFF limits for a control name; please choose something else. \n")
    end
    push!(n.control.name,names[i])
    if !isempty(descriptions)
      push!(n.control.description,descriptions[i])
    end
  end
  return nothing
end

########################################################################################
# data related functions
########################################################################################

"""
description = string(
" *  \n ")

Dir!(r;results_name,description=description)
# removes results folder and creates a new one
# TODO consider putting in a warning or some sort of an interaction with user
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/26/2017, Last Modified: 5/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function resultsDir!(n;results_name::String = "",description::DataFrame = DataFrame())
 results_dir=string(n.r.main_dir,"/results/",results_name)  # define directories
 n.r.results_dir=results_dir;

 if isdir(n.r.results_dir)
   rm(n.r.results_dir; recursive=true)
   print("\n The old results have all been deleted! \n \n")
 end
 mkdir(n.r.results_dir)# create directory

 cd(n.r.results_dir)
   CSV.write("description.csv", description)
 cd(n.r.main_dir)
 return nothing
end


"""
------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/26/2017, Last Modified: 2/6/2018 \n
--------------------------------------------------------------------------------------\n
"""
function savePlantData!(n)
  cd(n.r.results_dir)
    CSV.write("plant_data.csv", n.r.dfs_plantPts);
  cd(n.r.main_dir)
  return nothing
end

"""
------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 10/3/2017, Last Modified: 11/11/2017 \n
--------------------------------------------------------------------------------------\n
"""

function saveData(n)
  # all polynomial data
  dfs=DataFrame();

  dfs[:t] = n.r.t_pts
  for st in 1:n.numStates # state
    dfs[n.state.name[st]]=n.r.X_pts[:,st];
  end

  for ctr in 1:n.numControls # control
    dfs[n.control.name[ctr]]=n.r.U_pts[:,ctr];
  end

  if n.s.evalCostates && n.s.integrationMethod == :ps && n.s.evalConstraints
    for st in 1:n.numStates # state
      dfs[Symbol(n.state.name[st],:cs)]=n.r.CS_pts[:,st];
    end
  end

  cd(n.r.results_dir)
    CSV.write("st_ctr.csv",n.r.dfs[end]); # assuming only want the last one is needed
    CSV.write("st_ctr_poly.csv",dfs);
  cd(n.r.main_dir)
  return nothing
end


"""
------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 9/14/2017, Last Modified: 9/19/2017 \n
--------------------------------------------------------------------------------------\n
"""
function saveBenchMarkData!(n)
  first=2
  dfs=DataFrame();
  temp = [n.r.dfs[jj][:t][1:end-1,:] for jj in first:length(n.r.dfs)]; # time
  U=[idx for tempM in temp for idx=tempM]; dfs[:t]=U;

  for st in 1:n.numStates # state
    temp = [n.r.dfs[jj][n.state.name[st]][1:end-1,:] for jj in first:length(n.r.dfs)];
    U=[idx for tempM in temp for idx=tempM];
    dfs[n.state.name[st]]=U;
  end

  for ctr in 1:n.numControls # control
    temp = [n.r.dfs[jj][n.control.name[ctr]][1:end-1,:] for jj in first:length(n.r.dfs)];
    U=[idx for tempM in temp for idx=tempM];
    dfs[n.control.name[ctr]]=U;
  end

  # save optimization times
  temp = [n.r.dfs_opt[jj][n.control.name[ctr]][1:end-1,:] for jj in first:length(n.r.dfs)];


  cd(n.r.results_dir)
    CSV.write("bench_data.csv",dfs);
  cd(n.r.main_dir)
  return nothing
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


function try_import(name::Symbol)
 try
      @eval using $name
      return true
  catch e
      return false
  end
end
