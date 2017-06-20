"""
n=define(;numStates=2,numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,-Inf],XU=[1/9,Inf],CL=[-Inf],CU=[Inf])
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 6/18/2017 \n
Citations: \n
----------\n
Initially Influenced by: S. Hughes.  steven.p.hughes@nasa.gov
Source: DecisionVector.m [located here](https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/@DecisionVector/)
-------------------------------------------------------------------------------------\n
"""
function define(de;
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
  if isa(de,Array)
    n.DXexpr=de;
    n.stateEquations=DiffEq;
  else
    n.stateEquations=de;
  end

  # validate input
  if  numStates <= 0
      error("\n numStates must be > 0","\n",
              "default value = 0","\n",
            );
  end
  if  numControls <= 0
      error("eventually numControls must be > 0","\n",
            "default value = 0","\n",
            );
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

  n.numStates = numStates;
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
defineSolver!(n;(:name=>:Ipopt))
# To debug KNITRO turn up the output level
# Try to tune KNITRO
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 5/29/2017 \n
-------------------------------------------------------------------------------------\n
"""
function defineSolver!(n::NLOpt;kwargs...)
  kw=Dict(kwargs);

  # get the name of the solver
  if haskey(kw,:name); n.s.solver.name=get(kw,:name,0); end

  # modify defaults
  for (key,value) in kw
    if haskey(n.s.solver.settings,key)
      n.s.solver.settings[key]=value
    elseif key!=:name # ignore the name
      error("Unknown key: ", kw)
    end
  end

  if try_import(n.s.solver.name)
  else error(string("\n could not import ",n.s.solver.name,"\n consider adding it with: \n Pkg.add(``",n.solver.name,``")") )
  end
  if n.s.solver.name==:Ipopt
    NLPsolver=Ipopt.IpoptSolver(print_level=0)
      #=
      NLPsolver=Ipopt.IpoptSolver(max_cpu_time=max_cpu_time,
                                 print_level=0,
                                 warm_start_init_point="yes",
                                 max_iter=max_iter,
                                 tol=infeastol,
                                 dual_inf_tol=1.,
                                 constr_viol_tol=0.0001,
                                 compl_inf_tol=1e-6,
                                 acceptable_tol=1e-6,
                                 acceptable_constr_viol_tol=0.01,
                                 acceptable_dual_inf_tol=1e10,
                                 acceptable_compl_inf_tol=0.01,
                                 acceptable_obj_change_tol=1e20,
                                 diverging_iterates_tol=1e20)
                                 =#
  elseif n.s.solver.name==:KNITRO
      NLPsolver=KNITRO.KnitroSolver(outlev=0)
  else
    error("the :name key needs to be set to either :KNITRO or :Ipopt in defineSolver!()\n ")
  end
  n.mdl=JuMP.Model(solver=NLPsolver)

  # optimal control problem
  OCPdef!(n);
  return nothing
end  # function

"""
OCPdef!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 6/14/2017 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef!(n::NLOpt)
  n.r=Result();

  # state variables
  @variable(n.mdl,x[1:n.numStatePoints,1:n.numStates]); n.r.x=x;
  for st in 1:n.numStates
    # lower state constraint
    if !isnan(n.XL[st])
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
    if !isnan(n.XU[st])
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
    if !isnan(n.CL[ctr])
      for j in 1:n.numControlPoints
        setlowerbound(n.r.u[j,ctr], n.CL[ctr])
      end
    end
    if !isnan(n.CU[ctr])
      for j in 1:n.numControlPoints
        setupperbound(n.r.u[j,ctr], n.CU[ctr])
      end
    end
  end

  # boundary constraints
  n.r.xf_con=[]; # currently modifying the final state constraint (with tolerance) is not needed, can easily ad this functionlity though
  if any(!isnan(n.X0_tol))           # create handles for constraining the enire initial state
    n.r.x0_con=Array(Any,n.numStates,2); # this is so they can be easily reference when doing MPC
  else
    n.r.x0_con=[];
  end

  for st in 1:n.numStates
    if !isnan(n.X0[st]) # could have a bool for this
      if any(!isnan(n.X0_tol)) #NOTE in JuMP: Modifying range constraints is currently unsupported.
        n.r.x0_con[st,1]=@constraint(n.mdl, n.r.x[1,st] <=  (n.X0[st]+n.X0_tol[st]));
        n.r.x0_con[st,2]=@constraint(n.mdl,-n.r.x[1,st] <= -(n.X0[st]-n.X0_tol[st]));
      else
        n.r.x0_con=[n.r.x0_con; @constraint(n.mdl, n.r.x[1,st]==n.X0[st])]
      end
    end
    if !isnan(n.XF[st])
      if isnan(n.XF_tol[st])
        n.r.xf_con=[n.r.xf_con; @constraint(n.mdl, n.r.x[end,st]==n.XF[st])];
      else #TODO fix this as well
        n.r.xf_con=[n.r.xf_con; @constraint(n.mdl, n.XF[st]-n.XF_tol[st] <= n.r.x[end,st] <= n.XF[st]+n.XF_tol[st])];
      end
    end
  end

  @NLparameter(n.mdl,t0_param==0.0);   # for now we just start at zero
  n.mpc.t0_param=t0_param;

  if n.s.integrationMethod==:ps
    n.r.dyn_con=[Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];
    Nck_full=[0;cumsum(n.Nck+1)]; Nck_vars=[0;cumsum(n.Nck)];
    dynamics_expr=[Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];

    if n.s.finalTimeDV
      @variable(n.mdl, 0.001 <= tf <=  n.s.tf_max)
      n.tf=tf;
      create_tV!(n) # make a time vector
    end

    for int in 1:n.Ni
      # states
      x_int=Array(JuMP.Variable,length(Nck_full[int]+1:Nck_full[int+1]),n.numStates);
      for st in 1:n.numStates # +1 adds the DV in the next interval
        x_int[:,st]=n.r.x[Nck_vars[int]+1:Nck_vars[int+1]+1,st];
      end

      # controls
      if int!=n.Ni
        u_int=Array(JuMP.Variable,length(Nck_full[int]+1:Nck_full[int+1]),n.numControls);
      else # -1 -> removing control in last mesh interval
        u_int=Array(JuMP.Variable,length(Nck_full[int]+1:Nck_full[int+1]-1),n.numControls);
      end
      for ctr in 1:n.numControls
        if int!=n.Ni # +1 adds the DV in the next interval
          u_int[:,ctr]=n.r.u[Nck_vars[int]+1:Nck_vars[int+1]+1,ctr];
        else
          u_int[:,ctr]=n.r.u[Nck_vars[int]+1:Nck_vars[int+1],ctr];
        end
      end

      # dynamics
      if isempty(n.DXexpr)
        dx=n.stateEquations(n,x_int,u_int);
      else
        L=size(x_int)[1]-1;
        dx=Array(Any,L,n.numStates)
        for st in 1:n.numStates
          dx[:,st]=n.stateEquations(n,x_int,u_int,L,st);
        end
      end

      for st in 1:n.numStates
        if n.s.integrationScheme==:lgrExplicit
          dynamics_expr[int][:,st]=@NLexpression(n.mdl, [j in 1:n.Nck[int]], sum(n.DMatrix[int][j,i]*x_int[i,st] for i in 1:n.Nck[int]+1) - ((n.tf)/2)*dx[j,st]  )
        end
        for j in 1:n.Nck[int]
          n.r.dyn_con[int][j,st]=@NLconstraint(n.mdl, 0==dynamics_expr[int][j,st])
        end
      end

    end
  elseif n.s.integrationMethod==:tm
    n.r.dyn_con=Array(Any,n.N,n.numStates);
    if n.s.finalTimeDV
     #@variable( mdl, 0.00001 <= dt[1:n.N] <= 0.2) #TODO allow for an varaible array of dts
     @variable(n.mdl, 0.001 <= tf <= n.s.tf_max)
     n.tf = tf;
    end
    n.dt = n.tf/n.N*ones(n.N,);

    if isempty(n.DXexpr)
      dx = n.stateEquations(n,n.r.x,n.r.u);
    else
      L=size(n.r.x)[1];
      dx=Array(Any,L,n.numStates)
      for st in 1:n.numStates
        dx[:,st]=n.stateEquations(n,n.r.x,n.r.u,L,st);
      end
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
Date Create: 1/1/2017, Last Modified: 6/14/2017 \n
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

  if n.s.integrationScheme==:lgrExplicit
    n.s.integrationMethod=:ps;
  elseif (n.s.integrationScheme==:trapezoidal) || (n.s.integrationScheme==:bkwEuler)
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
     n.numPoints = [n.Nck[int] for int in 1:n.Ni];  # number of design variables per interval
     n.numStatePoints=sum(n.Nck)+1;
     n.numControlPoints=sum(n.Nck);

    # initialize node data
    if n.s.integrationScheme==:lgrExplicit
      taus_and_weights = [gaussradau(n.Nck[int]) for int in 1:n.Ni];
    end
    n.τ=[taus_and_weights[int][1] for int in 1:n.Ni];
    n.ω=[taus_and_weights[int][2] for int in 1:n.Ni];
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

  # default solver
  defineSolver!(n;name=:Ipopt);

  return nothing
end
