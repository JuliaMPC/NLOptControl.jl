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
                                mXL::Array{Any,1}=falses(n.ocp.state.num),
                                mXU::Array{Any,1}=falses(n.ocp.state.num))
  n.ocp.mXL=mXL;  n.ocp.mXU=mXU;
  for st in 1:n.ocp.state.num
    # lower state constraint
    if n.ocp.mXL[st]!=false
      if !isnan(n.ocp.XL[st])
        for j in 1:n.ocp.state.pts
        n.ocp.XL_var[st,j]=n.ocp.XL[st] + n.ocp.mXL[st]*(j/n.ocp.state.pts); # lower
        end
      end
    end

    # upper state constraint
    if n.ocp.XU[st]!=false
      if !isnan(n.ocp.XU[st])
        for j in 1:n.ocp.state.pts
        n.ocp.XU_var[st,j]=n.ocp.XU[st] + n.ocp.mXU[st]*(j/n.ocp.state.pts); # upper
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
                          X0_tol::Array{Float64,1}=0.05*ones(Float64,n.ocp.state.num,),
                          XF_tol::Array{Float64,1}=0.05*ones(Float64,n.ocp.state.num,))
  n.ocp.X0_tol=X0_tol; n.ocp.XF_tol=XF_tol;
  return nothing
end

"""
create_tV!(n)
# define a time vector (n.ocp.tV) for use with time varying constraints when (finalTimeDV=>true)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/17/2017, Last Modified: 3/17/2017 \n
--------------------------------------------------------------------------------------\n
"""
function create_tV!(n::NLOpt)

  if n.s.ocp.integrationMethod==:ps
    # create mesh points, interval size = tf_var/Ni
    tm = @NLexpression(n.ocp.mdl, [idx=1:n.ocp.Ni+1], (idx-1)*n.ocp.tf/n.ocp.Ni);
    # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
    ts = [Array{Any}(n.ocp.Nck[int]+1,) for int in 1:n.ocp.Ni];
    for int in 1:n.ocp.Ni
      ts[int][1:end-1]=@NLexpression(n.ocp.mdl,[j=1:n.ocp.Nck[int]], (tm[int+1]-tm[int])/2*n.ocp.tau[int][j] +  (tm[int+1]+tm[int])/2);
      ts[int][end]=@NLexpression(n.ocp.mdl, n.ocp.tf/n.ocp.Ni*int) # append +1 at end of each interval
    end
    tt1 = [idx for tempM in ts for idx = tempM[1:end-1]];
    tmp = [tt1;ts[end][end]];
    n.ocp.tV = @NLexpression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.mpc.v.t0Param + tmp[j]);
  else
    # create vector with the design variable in it
    t = Array{Any}(n.ocp.N+1,1);
    tm = @NLexpression(n.ocp.mdl, [idx=1:n.ocp.N], n.ocp.tf/n.ocp.N*idx);
    tmp = [0;tm];
    n.ocp.tV = @NLexpression(n.ocp.mdl,[j=1:n.ocp.state.pts], n.mpc.v.t0Param + tmp[j]);
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
  if n.s.ocp.integrationMethod==:ps
    integral_expr = [Array{Any}(n.ocp.Nck[int]) for int in 1:n.ocp.Ni]
    for int in 1:n.ocp.Ni
      x_int,u_int = intervals(n,int,n.r.ocp.x,n.r.ocp.u)
      L = size(x_int)[1]-1
      integral_expr[int][:] = NLExpr(n,V,x_int,u_int,L)
    end
    @NLexpression(n.ocp.mdl, temp[int=1:n.ocp.Ni], (n.ocp.tf-n.ocp.t0)/2*sum(n.ocp.ws[int][j]*integral_expr[int][j] for j = 1:n.ocp.Nck[int]) )
    expression = @NLexpression(n.ocp.mdl, sum(temp[int] for int = 1:n.ocp.Ni))
  elseif n.s.ocp.integrationMethod==:tm
    L = size(n.r.ocp.x)[1];
    temp = NLExpr(n,V,n.r.ocp.x,n.r.ocp.u,L);
    if n.s.ocp.integrationScheme==:bkwEuler
      expression = @NLexpression(n.ocp.mdl, sum(temp[j+1]*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    elseif n.s.ocp.integrationScheme==:trapezoidal
      expression = @NLexpression(n.ocp.mdl, sum(0.5*(temp[j]+temp[j+1])*n.ocp.tf/n.ocp.N for j = 1:n.ocp.N) )
    else
      error("\n $(n.s.ocp.integrationScheme) not defined in integrationSchemes\n")
    end
  else
    error("\n $(n.s.ocp.integrationMethod) not defined in integrationMethods \n")
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
  if n.r.ocp.constraint==nothing
    n.r.ocp.constraint=Constraint()
  end
  return nothing
end

function newConstraint!(n,handle,name::Symbol)
  initConstraint!(n)
  n.r.ocp.constraint::Constraint=n.r.ocp.constraint
  push!(n.r.ocp.constraint.handle,handle)
  push!(n.r.ocp.constraint.name,name)
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
  num=length(n.r.ocp.constraint.handle); dual_con_temp=zeros(num);
  for i = 1:length(n.r.ocp.constraint.handle)
    if !isempty(n.r.ocp.constraint.handle[i])
      if n.r.ocp.constraint.name[i]==:dyn_con  # state constraits
        temp1=zeros(n.ocp.state.num);
        for st in 1:n.ocp.state.num
          if n.s.ocp.integrationMethod==:ps
            temp=[getdual(n.r.ocp.constraint.handle[i][int][:,st]) for int in 1:n.ocp.Ni];
            vals=[idx for tempM in temp for idx=tempM];
            temp1[st]=maximum(vals);
          else
            temp1[st] = maximum(getdual(n.r.ocp.constraint.handle[i][:,st]));
          end
        end
        dual_con_temp[i]=maximum(temp1);
      else
        S=JuMP.size(n.r.ocp.constraint.handle[i])
        if length(S)==1
          dual_con_temp[i]=maximum(getdual(n.r.ocp.constraint.handle[i][:]));
        elseif length(S)==2
          temp1=zeros(S[1]);
          for idx in 1:S[1]
            temp1[idx] = maximum(getdual(n.r.ocp.constraint.handle[i][idx,:]));
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
# initialize states
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function initState(numStates)
  s = State()
  s.num = numStates
  s.name = [Symbol("x$i") for i in 1:numStates]
  s.description = [String("x$i") for i in 1:numStates]
  return s
end

function states!(n::NLOpt,names;descriptions=[])
  if !n.f.ocp.defined
    error("\n call define() before calling states!() \n")
  end
  if length(names)!=n.ocp.state.num
    error("\n Check size of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.ocp.state.num
    error("\n Check size of descriptions \n")
  end

  for i in 1:n.ocp.state.num
    if names[i]==:xxx
      error("xxx is OFF limits for a state name; please choose something else. \n")
    end
  end

   n.ocp.state.name = names
   if !isempty(descriptions)
     n.ocp.state.description = descriptions
   end
  return nothing
end

########################################################################################
# control data functions
########################################################################################
"""
# initialize control
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/7/2017, Last Modified: 3/6/2017 \n
--------------------------------------------------------------------------------------\n
"""
function initControl(numControls)
  c = Control()
  c.num = numControls
  c.name = [Symbol("u$i") for i in 1:numControls]
  c.description = [String("u$i") for i in 1:numControls]
  return c
end

function controls!(n::NLOpt,names;descriptions=[])
  if !n.f.ocp.defined
    error("\n call define() before calling controls!() \n")
  end
  if length(names)!=n.ocp.control.num
    error("\n Check sizes of names \n")
  end
  if !isempty(descriptions) && length(descriptions)!=n.ocp.control.num
    error("\n Check size of descriptions \n")
  end

  for i in 1:n.ocp.control.num
    if names[i]==:uuu
      error("uuu is OFF limits for a control name; please choose something else. \n")
    end
   end

  n.ocp.control.name = names
  if !isempty(descriptions)
    n.ocp.control.description = descriptions
  end
  return nothing
end

########################################################################################
# data related functions
########################################################################################

"""
description = string(
" *  \n ")

resultsDir!(r;resultsName=resultsName,description=description)
# removes results folder and creates a new one
# TODO consider putting in a warning or some sort of an interaction with user
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/26/2017, Last Modified: 5/29/2017 \n
--------------------------------------------------------------------------------------\n
"""
function resultsDir!(n;resultsName::String = "",description::DataFrame = DataFrame())
 results_dir=string(n.r.mainDir,"/results/",resultsName)  # define directories
 n.r.resultsDir=results_dir;

 if isdir(n.r.resultsDir)
   rm(n.r.resultsDir; recursive=true)
   print("\n The old results have all been deleted! \n \n")
 end
 mkdir(n.r.resultsDir)# create directory

 cd(n.r.resultsDir)
   CSV.write("description.csv", description; quotechar = ' ')
 cd(n.r.mainDir)
 return nothing
end


"""
------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 3/26/2017, Last Modified: 2/6/2018 \n
--------------------------------------------------------------------------------------\n
"""
function savePlantData!(n)
  cd(n.r.resultsDir)
    CSV.write("plant_data.csv", n.r.ip.dfsplantPts; quotechar = ' ');
  cd(n.r.mainDir)
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

  dfs[:t] = n.r.ocp.tpts
  for st in 1:n.ocp.state.num # state
    dfs[n.ocp.state.name[st]]=n.r.ocp.Xpts[:,st];
  end

  for ctr in 1:n.ocp.control.num # control
    dfs[n.ocp.control.name[ctr]]=n.r.ocp.Upts[:,ctr];
  end

  if n.s.ocp.evalCostates && n.s.ocp.integrationMethod == :ps && n.s.ocp.evalConstraints
    for st in 1:n.ocp.state.num # state
      dfs[Symbol(n.ocp.state.name[st],:cs)]=n.r.ocp.CSpts[:,st];
    end
  end

  cd(n.r.resultsDir)
    CSV.write("st_ctr.csv",n.r.ocp.dfs[end]; quotechar = ' '); # assuming only want the last one is needed
    CSV.write("st_ctr_poly.csv",dfs; quotechar = ' ');
  cd(n.r.mainDir)
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
  temp = [n.r.ocp.dfs[jj][:t][1:end-1,:] for jj in first:length(n.r.ocp.dfs)]; # time
  U=[idx for tempM in temp for idx=tempM]; dfs[:t]=U;

  for st in 1:n.ocp.state.num # state
    temp = [n.r.ocp.dfs[jj][n.ocp.state.name[st]][1:end-1,:] for jj in first:length(n.r.ocp.dfs)];
    U=[idx for tempM in temp for idx=tempM];
    dfs[n.ocp.state.name[st]]=U;
  end

  for ctr in 1:n.ocp.control.num # control
    temp = [n.r.ocp.dfs[jj][n.ocp.control.name[ctr]][1:end-1,:] for jj in first:length(n.r.ocp.dfs)];
    U=[idx for tempM in temp for idx=tempM];
    dfs[n.ocp.control.name[ctr]]=U;
  end

  # save optimization times
  temp = [n.r.ocp.dfsOpt[jj][n.ocp.control.name[ctr]][1:end-1,:] for jj in first:length(n.r.ocp.dfs)];

  cd(n.r.resultsDir)
    CSV.write("bench_data.csv",dfs; quotechar = ' ');
  cd(n.r.mainDir)
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
