"""
n=define(;numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,-Inf],XU=[1/9,Inf],CL=[-Inf],CU=[Inf])
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 6/18/2017 \n
Citations: \n
----------\n
Initially Influenced by: S. Hughes.  steven.p.hughes@nasa.gov
Source: DecisionVector.m [located here](https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/@DecisionVector/)
-------------------------------------------------------------------------------------\n
"""
function define(;
                numStates::Int64=0,
                numControls::Int64=0,
                X0=fill(NaN,numStates,),
                XF=fill(NaN,numStates,),
                XL=fill(NaN,numStates,),
                XU=fill(NaN,numStates,),
                CL=fill(NaN,numControls,),
                CU=fill(NaN,numControls,)
                )

  n=NLOpt();

  # validate input
  if numControls <= 0
      error("numControls must be > 0","\n");
  end
  if numStates <= 0
      error("numStates must be > 0","\n");
  end
  if length(X0) != numStates
    error(string("\n Length of X0 must match number of states \n"));
  end
  if length(XF) != numStates
    error(string("\n Length of XF must match number of states \n"));
  end
  if length(XL) != numStates
    error(string("\n Length of XL must match number of states \n"));
  end
  if length(XU) != numStates
    error(string("\n Length of XU must match number of states \n"));
  end
  if length(CL) != numControls
    error(string("\n Length of CL must match number of controls \n"));
  end
  if length(CU) != numControls
    error(string("\n Length of CU must match number of controls \n"));
  end

  n.numStates=numStates;
  n.numControls = numControls;
  n.state=initStateNames(n);
  n.control=initControlNames(n);
  n.X0 = X0;
  n.X0_tol = NaN*X0;
  n.XF = XF;
  n.XF_tol = NaN*XF;
  n.XL = XL;
  n.XU = XU;
  n.CL = CL;
  n.CU = CU;
  n.define=true;
  return n
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 7/04/2017 \n
-------------------------------------------------------------------------------------\n
"""
function defineSolver!(n::NLOpt,kw)

  # get the name of the solver
  if haskey(kw,:name); n.s.solver.name=get(kw,:name,0); end
  if try_import(n.s.solver.name)
  else error(string("\n could not import ",n.s.solver.name,"\n consider adding it with: \n Pkg.add(``",n.s.solver.name,``")") )
  end

  # see if the user would like to use a standard set of solver settings for mpc
  if haskey(kw,:mpc_defaults); mpc_defaults=get(kw,:mpc_defaults,0); else mpc_defaults=false; end
  if mpc_defaults
    if n.s.solver.name==:Ipopt
      n.s.solver.settings=_Ipopt_MPC;
    elseif n.s.solver.name==:KNITRO
      n.s.solver.settings=_KNITRO_MPC;
    else
      error(string("solver ",n.s.sover.name, " not defined"))
    end
  else # default solver settings
    if n.s.solver.name==:Ipopt  # NOTE this should already have been done by default, but could get messed up is user is playing with options
      n.s.solver.settings=_Ipopt_defaults;
    elseif n.s.solver.name==:KNITRO
      n.s.solver.settings=_KNITRO_defaults;
    else
      error(string("solver ", n.s.sover.name, " not defined"))
    end
  end

  # modify additional defaults individually
  for (key,value) in kw
    if haskey(n.s.solver.settings,key)
      n.s.solver.settings[key]=value
    elseif key!=:name && key!=:mpc_defaults # ignore the name and default settings option TODO could remove them from the Dict
      error(string(" \n Unknown key: ", kw, " for ", n.s.solver.name, " used in defineSolver!() \n "))
    end
  end

  if n.s.solver.name==:Ipopt
    setsolver(n.mdl,Ipopt.IpoptSolver(;max_cpu_time=n.s.solver.settings[:max_cpu_time],
                               print_level=n.s.solver.settings[:print_level],
                               warm_start_init_point=n.s.solver.settings[:warm_start_init_point],
                               max_iter=n.s.solver.settings[:max_iter],
                               tol=n.s.solver.settings[:tol],
                               dual_inf_tol=n.s.solver.settings[:dual_inf_tol],
                               constr_viol_tol=n.s.solver.settings[:constr_viol_tol],
                               compl_inf_tol=n.s.solver.settings[:compl_inf_tol],
                               acceptable_tol=n.s.solver.settings[:acceptable_tol],
                               acceptable_constr_viol_tol=n.s.solver.settings[:acceptable_constr_viol_tol],
                               acceptable_dual_inf_tol=n.s.solver.settings[:acceptable_dual_inf_tol],
                               acceptable_compl_inf_tol=n.s.solver.settings[:acceptable_compl_inf_tol],
                               acceptable_obj_change_tol=n.s.solver.settings[:acceptable_obj_change_tol],
                               diverging_iterates_tol=n.s.solver.settings[:diverging_iterates_tol]))
  elseif n.s.solver.name==:KNITRO
    setsolver(n.mdl,KnitroSolver(;outlev=n.s.solver.settings[:outlev],
                                 maxit=n.s.solver.settings[:maxit],
                                 maxtime_real=n.s.solver.settings[:maxtime_real],
                                 feastol=n.s.solver.settings[:feastol],
                                 feastol_abs=n.s.solver.settings[:feastol_abs],
                                 ftol=n.s.solver.settings[:ftol],
                                 ftol_iters=n.s.solver.settings[:ftol_iters],
                                 infeastol=n.s.solver.settings[:infeastol],
                                 maxfevals=n.s.solver.settings[:maxfevals],
                                 maxit=n.s.solver.settings[:maxit],
                                 maxtime_cpu=n.s.solver.settings[:maxtime_cpu],
                                 maxtime_real=n.s.solver.settings[:maxtime_real],
                                 opttol=n.s.solver.settings[:opttol],
                                 opttol_abs=n.s.solver.settings[:opttol_abs],
                                 xtol=n.s.solver.settings[:xtol],
                                 xtol_iters=n.s.solver.settings[:xtol_iters],
                                 algorithm=n.s.solver.settings[:algorithm],
                                 bar_initpt=n.s.solver.settings[:bar_initpt],
                                 bar_murule=n.s.solver.settings[:bar_murule],
                                 bar_penaltycons=n.s.solver.settings[:bar_penaltycons],
                                 bar_penaltyrule=n.s.solver.settings[:bar_penaltyrule],
                                 bar_switchrule=n.s.solver.settings[:bar_switchrule],
                                 linesearch=n.s.solver.settings[:linesearch],
                                 linsolver=n.s.solver.settings[:linsolver]))
  else
    error(string("solver ",n.s.sover.name, " not defined"))
  end
  return nothing
end  # function

"""
OCPdef!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 7/1/2017 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef!(n::NLOpt)
  n.r=Result();

  # state variables
  @variable(n.mdl,x[1:n.numStatePoints,1:n.numStates]); n.r.x=x;
  for st in 1:n.numStates
    # lower state constraint
    if !isnan.(n.XL[st])
      if n.mXL[st]==false
        for j in 1:n.numStatePoints
          setlowerbound(n.r.x[j,st], n.XL[st])
        end
      else
        for j in 1:n.numStatePoints
          setlowerbound(n.r.x[j,st],n.XL_var[st,j])
        end
      end
    end

    # upper state constraint
    if !isnan.(n.XU[st])
      if n.XU[st]!=false
        for j in 1:n.numStatePoints
          setupperbound(n.r.x[j,st], n.XU[st])
        end
      else
        for j in 1:n.numStatePoints
          setlowerbound(n.r.x[j,st],n.XU_var[st,j])
        end
      end
    end
  end

  # control variables
  @variable(n.mdl,u[1:n.numControlPoints,1:n.numControls]);n.r.u=u;
  for ctr in 1:n.numControls
    if !isnan.(n.CL[ctr])
      for j in 1:n.numControlPoints
        setlowerbound(n.r.u[j,ctr], n.CL[ctr])
      end
    end
    if !isnan.(n.CU[ctr])
      for j in 1:n.numControlPoints
        setupperbound(n.r.u[j,ctr], n.CU[ctr])
      end
    end
  end

  # boundary constraints
  if any(.!isnan.(n.X0_tol))              # create handles for constraining the entire initial state
    n.r.x0_con=Array{Any}(n.numStates,2); # this is so they can be easily reference when doing MPC
  else
    n.r.x0_con=[];
  end

  if any(.!isnan.(n.XF_tol))              # create handles for constraining the entire final state
    n.r.xf_con=Array{Any}(n.numStates,2); # this is so they can be easily reference when doing MPC
  else
    n.r.xf_con=[];
  end

  for st in 1:n.numStates
    if !isnan(n.X0[st]) # could have a bool for this
      if any(.!isnan.(n.X0_tol)) #NOTE in JuMP: Modifying range constraints is currently unsupported.
        n.r.x0_con[st,1]=@constraint(n.mdl, n.r.x[1,st] <=  (n.X0[st]+n.X0_tol[st]));
        n.r.x0_con[st,2]=@constraint(n.mdl,-n.r.x[1,st] <= -(n.X0[st]-n.X0_tol[st]));
      else
        n.r.x0_con=[n.r.x0_con; @constraint(n.mdl, n.r.x[1,st]==n.X0[st])]
      end
    end
    if !isnan(n.XF[st])
      if any(.!isnan.(n.XF_tol))
        n.r.xf_con[st,1]=@constraint(n.mdl, n.r.x[end,st] <=  (n.XF[st]+n.XF_tol[st]));
        n.r.xf_con[st,2]=@constraint(n.mdl,-n.r.x[end,st] <= -(n.XF[st]-n.XF_tol[st]));
      else
        n.r.xf_con=[n.r.xf_con; @constraint(n.mdl, n.r.x[end,st]==n.XF[st])];
      end
    end
  end

  @NLparameter(n.mdl,t0_param==0.0);   # for now we just start at zero
  n.mpc.t0_param=t0_param;

  if n.s.integrationMethod==:ps
    n.r.dyn_con=[Array{Any}(n.Nck[int],n.numStates) for int in 1:n.Ni];
    dynamics_expr=[Array{Any}(n.Nck[int],n.numStates) for int in 1:n.Ni];

    if n.s.finalTimeDV
      @variable(n.mdl, 0.001 <= tf <=  n.s.tf_max)
      n.tf=tf;
      create_tV!(n)          # make a time vector
    end

    for int in 1:n.Ni
      x_int,u_int=intervals(n,int,n.r.x,n.r.u);

      # dynamics
      L=size(x_int)[1]-1;
      dx=Array{Any}(L,n.numStates)
      for st in 1:n.numStates
        dx[:,st]=DiffEq(n,x_int,u_int,L,st);
      end

      for st in 1:n.numStates # TODO consider multiplying X*D to reduce computations (i.e. remove this for loop for the states)
        if n.s.integrationScheme==:lgrExplicit
          dynamics_expr[int][:,st]=@NLexpression(n.mdl, [j in 1:n.Nck[int]], sum(n.DMatrix[int][j,i]*x_int[i,st] for i in 1:n.Nck[int]+1) - ((n.tf)/2)*dx[j,st]  )
        elseif n.s.integrationScheme==:lgrImplicit
          dynamics_expr[int][:,st]=@NLexpression(n.mdl, [j in 1:n.Nck[int]], x_int[j+1,st] - x_int[1,st] - ((n.tf)/2)*sum(n.IMatrix[int][j,i]*dx[i,st] for i in 1:n.Nck[int]) )
        end
        for j in 1:n.Nck[int]
          n.r.dyn_con[int][j,st]=@NLconstraint(n.mdl, 0==dynamics_expr[int][j,st])
        end
      end

      # additional constraints
      for num in 1:length(n.NLcon)
        addCon(n,x_int,u_int,L,num);
      end
    end
  elseif n.s.integrationMethod==:tm
    n.r.dyn_con=Array{Any}(n.N,n.numStates);
    if n.s.finalTimeDV
     @variable(n.mdl, 0.001 <= tf <= n.s.tf_max)
     n.tf = tf;
    end
    n.dt = n.tf/n.N*ones(n.N,);

    L=size(n.r.x)[1];
    dx=Array{Any}(L,n.numStates)
    for st in 1:n.numStates
      dx[:,st]=DiffEq(n,n.r.x,n.r.u,L,st);
    end

    if n.s.integrationScheme==:bkwEuler
      for st in 1:n.numStates
        n.r.dyn_con[:,st] = @NLconstraint(n.mdl, [j in 1:n.N], n.r.x[j+1,st] - n.r.x[j,st] ==  dx[j+1,st]*n.tf/(n.N) );
      end
    elseif n.s.integrationScheme==:trapezoidal
      for st in 1:n.numStates
        n.r.dyn_con[:,st] = @NLconstraint(n.mdl, [j in 1:n.N], n.r.x[j+1,st] - n.r.x[j,st] == 0.5*(dx[j,st] + dx[j+1,st])*n.tf/(n.N) )
      end
    end
  end

  # save constraint data
  newConstraint!(n,n.r.x0_con,:x0_con);
  newConstraint!(n,n.r.xf_con,:xf_con);
  newConstraint!(n,n.r.dyn_con,:dyn_con);

  # save the current working directory for navigation purposes
  n.r.main_dir=pwd();
  return nothing
end

"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 11/10/2017 \n
-------------------------------------------------------------------------------------\n
"""
function configure!(n::NLOpt; kwargs... )
  kw = Dict(kwargs);

  # final time
  if !haskey(kw,:finalTimeDV);n.s.finalTimeDV=false;
  else; n.s.finalTimeDV=get(kw,:finalTimeDV,0);
  end

  if !haskey(kw,:tf) && !n.s.finalTimeDV
    error("\n If the final is not a design variable pass it as: (:tf=>Float64(some #)) \n
        If the final time is a design variable, indicate that as: (:finalTimeDV=>true)\n")
  elseif haskey(kw,:tf) && !n.s.finalTimeDV
    n.tf = get(kw,:tf,0);
  elseif n.s.finalTimeDV
    n.tf = Any;
  end

  # integrationScheme
  if !haskey(kw,:integrationScheme); n.s.integrationScheme=:lgrExplicit; # default
  else; n.s.integrationScheme=get(kw,:integrationScheme,0);
  end

  if n.s.integrationScheme==:lgrExplicit ||  n.s.integrationScheme==:lgrImplicit
    n.s.integrationMethod=:ps;
  elseif n.s.integrationScheme==:trapezoidal || n.s.integrationScheme==:bkwEuler
    n.s.integrationMethod=:tm;
  else
    error("the :integrationScheme that you specified is not currently implemeted \n")
  end

  if n.s.integrationMethod==:ps
    if haskey(kw,:N)
      error(" \n N is not an appropriate kwargs for :ps methods \n")
    end
    if !haskey(kw,:Nck);n.Nck=[10,10,10,10]; # default
    else; n.Nck = get(kw,:Nck,0);
    end
    n.Ni=length(n.Nck);

    for int in 1:n.Ni
        if (n.Nck[int]<0)
            error("\n Nck must be > 0");
        end
    end
    n.numStatePoints=sum(n.Nck)+1;
    n.numControlPoints=sum(n.Nck);
    n.Nck_full=[0;cumsum(n.Nck+1)];
    n.Nck_cum=[0;cumsum(n.Nck)];

    # initialize node data
    if n.s.integrationScheme==:lgrExplicit ||  n.s.integrationScheme==:lgrImplicit
      taus_and_weights = [gaussradau(n.Nck[int]) for int in 1:n.Ni];
    end
    n.tau=[taus_and_weights[int][1] for int in 1:n.Ni];
    n.w=[taus_and_weights[int][2] for int in 1:n.Ni];
    createIntervals!(n);
    DMatrix!(n);

  elseif n.s.integrationMethod==:tm
    if haskey(kw,:Nck)
      error(" \n Nck is not appropriate kwargs for :tm methods \n")
    end
    if !haskey(kw,:N);n.N=100; # default
    else; n.N = get(kw,:N,0);
    end
    n.numStatePoints=n.N+1;
    n.numControlPoints=n.N+1;
  end
  n.mXL=falses(n.numStates);
  n.mXU=falses(n.numStates);
  n.XL_var=Matrix{Float64}(n.numStates,n.numStatePoints);
  n.XU_var=Matrix{Float64}(n.numStates,n.numStatePoints);

  # solver settings
  if !haskey(kw,:solverSettings);SS=Dict((:name=>:Ipopt)); # default
  else; SS=get(kw,:solverSettings,0); SS=Dict(SS);
  end
  defineSolver!(n,SS)

  # optimal control problem
  OCPdef!(n);

  if n.s.evalCostates
     if n.s.integrationMethod != :ps
       error("costates are only implmented for :ps methods")
     end
     n.s.evalConstraints = true
  end
  return nothing
end
