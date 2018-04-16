"""
n=define(;numControls=1,X0=[0.,1],XF=[0.,-1.],XL=[0.,-Inf],XU=[1/9,Inf],CL=[-Inf],CU=[Inf])
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 4/13/2018 \n
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

  n = NLOpt()

  # validate input
  if numControls <= 0
      error("numControls must be > 0","\n")
  end
  if numStates <= 0
      error("numStates must be > 0","\n")
  end
  if length(X0) != numStates
    error(string("\n Length of X0 must match number of states \n"))
  end
  if length(XF) != numStates
    error(string("\n Length of XF must match number of states \n"))
  end
  if length(XL) != numStates
    error(string("\n Length of XL must match number of states \n"))
  end
  if length(XU) != numStates
    error(string("\n Length of XU must match number of states \n"))
  end
  if length(CL) != numControls
    error(string("\n Length of CL must match number of controls \n"))
  end
  if length(CU) != numControls
    error(string("\n Length of CU must match number of controls \n"))
  end

  n.ocp.state = initState(numStates)
  n.ocp.control = initControl(numControls)
  n.ocp.X0 = X0
  n.ocp.X0_tol = NaN*X0
  n.ocp.XF = XF
  n.ocp.XF_tol = NaN*XF
  n.ocp.XL = XL
  n.ocp.XU = XU
  n.ocp.CL = CL
  n.ocp.CU = CU
  n.f.ocp.defined = true
  return n
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 4/13/2018 \n
-------------------------------------------------------------------------------------\n
"""
function defineSolver!(n::NLOpt,kw)

 if typeof(kw)!=Dict{Symbol,Symbol}
   kw = Dict(kw)
 end

  # get the name of the solver
  if haskey(kw,:name); n.s.ocp.solver.name=get(kw,:name,0); end
  if try_import(n.s.ocp.solver.name)
  else error(string("could not import ",n.s.ocp.solver.name) )
  end

  # see if the user would like to use a standard set of solver settings for mpc
  if haskey(kw,:mpc_defaults); mpc_defaults=get(kw,:mpc_defaults,0); else mpc_defaults=false; end
  if mpc_defaults
    if n.s.ocp.solver.name==:Ipopt
      n.s.ocp.solver.settings=_Ipopt_MPC;
    elseif n.s.ocp.solver.name==:KNITRO
      n.s.ocp.solver.settings=_KNITRO_MPC;
    else
      error(string("solver ",n.s.sover.name, " not defined"))
    end
  else # default solver settings
    if n.s.ocp.solver.name==:Ipopt  # NOTE this should already have been done by default, but could get messed up is user is playing with options
      n.s.ocp.solver.settings=_Ipopt_defaults;
    elseif n.s.ocp.solver.name==:KNITRO
      n.s.ocp.solver.settings=_KNITRO_defaults;
    else
      error(string("solver ", n.s.sover.name, " not defined"))
    end
  end

  # modify additional defaults individually
  for (key,value) in kw
    if haskey(n.s.ocp.solver.settings,key)
      n.s.ocp.solver.settings[key]=value
    elseif key!=:name && key!=:mpc_defaults # ignore the name and default settings option TODO could remove them from the Dict
      error(string(" \n Unknown key: ", kw, " for ", n.s.ocp.solver.name, " used in defineSolver!() \n "))
    end
  end

  if n.s.ocp.solver.name==:Ipopt
    setsolver(n.ocp.mdl,Ipopt.IpoptSolver(;max_cpu_time=n.s.ocp.solver.settings[:max_cpu_time],
                               print_level=n.s.ocp.solver.settings[:print_level],
                               warm_start_init_point=n.s.ocp.solver.settings[:warm_start_init_point],
                               max_iter=n.s.ocp.solver.settings[:max_iter],
                               tol=n.s.ocp.solver.settings[:tol],
                               dual_inf_tol=n.s.ocp.solver.settings[:dual_inf_tol],
                               constr_viol_tol=n.s.ocp.solver.settings[:constr_viol_tol],
                               compl_inf_tol=n.s.ocp.solver.settings[:compl_inf_tol],
                               acceptable_tol=n.s.ocp.solver.settings[:acceptable_tol],
                               acceptable_constr_viol_tol=n.s.ocp.solver.settings[:acceptable_constr_viol_tol],
                               acceptable_dual_inf_tol=n.s.ocp.solver.settings[:acceptable_dual_inf_tol],
                               acceptable_compl_inf_tol=n.s.ocp.solver.settings[:acceptable_compl_inf_tol],
                               acceptable_obj_change_tol=n.s.ocp.solver.settings[:acceptable_obj_change_tol],
                               diverging_iterates_tol=n.s.ocp.solver.settings[:diverging_iterates_tol]))
  elseif n.s.ocp.solver.name==:KNITRO
    setsolver(n.ocp.mdl,KnitroSolver(;outlev=n.s.ocp.solver.settings[:outlev],
                                 maxit=n.s.ocp.solver.settings[:maxit],
                                 maxtime_real=n.s.ocp.solver.settings[:maxtime_real],
                                 feastol=n.s.ocp.solver.settings[:feastol],
                                 feastol_abs=n.s.ocp.solver.settings[:feastol_abs],
                                 ftol=n.s.ocp.solver.settings[:ftol],
                                 ftol_iters=n.s.ocp.solver.settings[:ftol_iters],
                                 infeastol=n.s.ocp.solver.settings[:infeastol],
                                 maxfevals=n.s.ocp.solver.settings[:maxfevals],
                                 maxit=n.s.ocp.solver.settings[:maxit],
                                 maxtime_cpu=n.s.ocp.solver.settings[:maxtime_cpu],
                                 maxtime_real=n.s.ocp.solver.settings[:maxtime_real],
                                 opttol=n.s.ocp.solver.settings[:opttol],
                                 opttol_abs=n.s.ocp.solver.settings[:opttol_abs],
                                 xtol=n.s.ocp.solver.settings[:xtol],
                                 xtol_iters=n.s.ocp.solver.settings[:xtol_iters],
                                 algorithm=n.s.ocp.solver.settings[:algorithm],
                                 bar_initpt=n.s.ocp.solver.settings[:bar_initpt],
                                 bar_murule=n.s.ocp.solver.settings[:bar_murule],
                                 bar_penaltycons=n.s.ocp.solver.settings[:bar_penaltycons],
                                 bar_penaltyrule=n.s.ocp.solver.settings[:bar_penaltyrule],
                                 bar_switchrule=n.s.ocp.solver.settings[:bar_switchrule],
                                 linesearch=n.s.ocp.solver.settings[:linesearch],
                                 linsolver=n.s.ocp.solver.settings[:linsolver]))
  else
    error(string("solver ",n.s.sover.name, " not defined"))
  end
  return nothing
end  # function

"""
OCPdef!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 4/13/2018 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef!(n::NLOpt)
  # state variables
  @variable(n.ocp.mdl,x[1:n.ocp.state.pts,1:n.ocp.state.num]); n.r.ocp.x=x;
  for st in 1:n.ocp.state.num
    # lower state constraint
    if !isnan.(n.ocp.XL[st])
      if n.ocp.mXL[st]==false
        for j in 1:n.ocp.state.pts
          setlowerbound(n.r.ocp.x[j,st], n.ocp.XL[st])
        end
      else
        for j in 1:n.ocp.state.pts
          setlowerbound(n.r.ocp.x[j,st],n.ocp.XL_var[st,j])
        end
      end
    end

    # upper state constraint
    if !isnan.(n.ocp.XU[st])
      if n.ocp.XU[st]!=false
        for j in 1:n.ocp.state.pts
          setupperbound(n.r.ocp.x[j,st], n.ocp.XU[st])
        end
      else
        for j in 1:n.ocp.state.pts
          setlowerbound(n.r.ocp.x[j,st],n.ocp.XU_var[st,j])
        end
      end
    end
  end

  # control variables
  @variable(n.ocp.mdl,u[1:n.ocp.control.pts,1:n.ocp.control.num]);n.r.ocp.u=u;
  for ctr in 1:n.ocp.control.num
    if !isnan.(n.ocp.CL[ctr])
      for j in 1:n.ocp.control.pts
        setlowerbound(n.r.ocp.u[j,ctr], n.ocp.CL[ctr])
      end
    end
    if !isnan.(n.ocp.CU[ctr])
      for j in 1:n.ocp.control.pts
        setupperbound(n.r.ocp.u[j,ctr], n.ocp.CU[ctr])
      end
    end
  end

  # boundary constraints
  if any(.!isnan.(n.ocp.X0_tol))              # create handles for constraining the entire initial state
    n.r.ocp.x0Con = Array{Any}(n.ocp.state.num,2) # this is so they can be easily reference when doing MPC
  else
    n.r.ocp.x0Con = []
  end

  if any(.!isnan.(n.ocp.XF_tol))              # create handles for constraining the entire final state
    n.r.ocp.xfCon = Array{Any}(n.ocp.state.num,2) # this is so they can be easily reference when doing MPC
  else
    n.r.ocp.xfCon = []
  end

  for st in 1:n.ocp.state.num
    if !isnan(n.ocp.X0[st]) # could have a bool for this
      if !isnan(n.ocp.X0_tol[st]) #NOTE in JuMP: Modifying range constraints is currently unsupported.
        n.r.ocp.x0Con[st,1] = @constraint(n.ocp.mdl, n.r.ocp.x[1,st] <=  (n.ocp.X0[st]+n.ocp.X0_tol[st]))
        n.r.ocp.x0Con[st,2] = @constraint(n.ocp.mdl,-n.r.ocp.x[1,st] <= -(n.ocp.X0[st]-n.ocp.X0_tol[st]))
      else
        n.r.ocp.x0Con = [n.r.ocp.x0Con; @constraint(n.ocp.mdl, n.r.ocp.x[1,st]==n.ocp.X0[st])]
      end
    end
    if !isnan(n.ocp.XF[st])
      if any(.!isnan.(n.ocp.XF_tol))
        n.r.ocp.xfCon[st,1] = @constraint(n.ocp.mdl, n.r.ocp.x[end,st] <=  (n.ocp.XF[st]+n.ocp.XF_tol[st]));
        n.r.ocp.xfCon[st,2] = @constraint(n.ocp.mdl,-n.r.ocp.x[end,st] <= -(n.ocp.XF[st]-n.ocp.XF_tol[st]));
      else
        n.r.ocp.xfCon = [n.r.ocp.xfCon; @constraint(n.ocp.mdl, n.r.ocp.x[end,st]==n.ocp.XF[st])];
      end
    end
  end

  @NLparameter(n.ocp.mdl,t0_param==0.0);   # for now we just start at zero
  n.ocp.t0 = t0_param  # NOTE consider making a constraint that t0 < tf

  if n.s.ocp.integrationMethod==:ps
    n.r.ocp.dynCon = [Array{Any}(n.ocp.Nck[int],n.ocp.state.num) for int in 1:n.ocp.Ni]
    dynamics_expr = [Array{Any}(n.ocp.Nck[int],n.ocp.state.num) for int in 1:n.ocp.Ni]

    if n.s.ocp.finalTimeDV
      @variable(n.ocp.mdl, 0.001 <= tf <=  n.s.ocp.tfMax)
      n.ocp.tf = tf
      create_tV!(n)          # make a time vector
    end

    for int in 1:n.ocp.Ni
      x_int,u_int = intervals(n,int,n.r.ocp.x,n.r.ocp.u)

      # dynamics
      L = size(x_int)[1]-1;
      dx = Array{Any}(L,n.ocp.state.num)
      for st in 1:n.ocp.state.num
        dx[:,st] = DiffEq(n,x_int,u_int,L,st)
      end

      for st in 1:n.ocp.state.num # TODO consider multiplying X*D to reduce computations (i.e. remove this for loop for the states)
        if n.s.ocp.integrationScheme==:lgrExplicit
          dynamics_expr[int][:,st] = @NLexpression(n.ocp.mdl, [j in 1:n.ocp.Nck[int]], sum(n.ocp.DMatrix[int][j,i]*x_int[i,st] for i in 1:n.ocp.Nck[int]+1) - ((n.ocp.tf)/2)*dx[j,st]  )
        elseif n.s.ocp.integrationScheme==:lgrImplicit
          dynamics_expr[int][:,st] = @NLexpression(n.ocp.mdl, [j in 1:n.ocp.Nck[int]], x_int[j+1,st] - x_int[1,st] - ((n.ocp.tf)/2)*sum(n.ocp.IMatrix[int][j,i]*dx[i,st] for i in 1:n.ocp.Nck[int]) )
        end
        for j in 1:n.ocp.Nck[int]
          n.r.ocp.dynCon[int][j,st] = @NLconstraint(n.ocp.mdl, 0==dynamics_expr[int][j,st])
        end
      end

      # additional constraints
      for num in 1:length(n.ocp.NLcon)
         ch = addCon(n,x_int,u_int,L,num)
         newConstraint!(n,ch,Symbol(string("ch",num))) # TODO could let the user name these
      end
    end
  elseif n.s.ocp.integrationMethod==:tm
    n.r.ocp.dynCon = Array{Any}(n.ocp.N,n.ocp.state.num)
    if n.s.ocp.finalTimeDV
     @variable(n.ocp.mdl, 0.001 <= tf <= n.s.ocp.tfMax)
     n.ocp.tf = tf
     create_tV!(n)          # make a time vector
    end
    n.ocp.dt = n.ocp.tf/n.ocp.N*ones(n.ocp.N,)

    L=size(n.r.ocp.x)[1];
    dx=Array{Any}(L,n.ocp.state.num)
    for st in 1:n.ocp.state.num
      dx[:,st]=DiffEq(n,n.r.ocp.x,n.r.ocp.u,L,st)
    end

    if n.s.ocp.integrationScheme==:bkwEuler
      for st in 1:n.ocp.state.num
        n.r.ocp.dynCon[:,st] = @NLconstraint(n.ocp.mdl, [j in 1:n.ocp.N], n.r.ocp.x[j+1,st] - n.r.ocp.x[j,st] ==  dx[j+1,st]*n.ocp.tf/(n.ocp.N) );
      end
    elseif n.s.ocp.integrationScheme==:trapezoidal
      for st in 1:n.ocp.state.num
        n.r.ocp.dynCon[:,st] = @NLconstraint(n.ocp.mdl, [j in 1:n.ocp.N], n.r.ocp.x[j+1,st] - n.r.ocp.x[j,st] == 0.5*(dx[j,st] + dx[j+1,st])*n.ocp.tf/(n.ocp.N) )
      end
    end

    # additional constraints
    for num in 1:length(n.ocp.NLcon)
       ch = addCon(n,n.r.ocp.x,n.r.ocp.u,L,num)
       newConstraint!(n,ch,Symbol(string("ch",num))) # TODO could let the user name these
    end
  end

  # save constraint data
  newConstraint!(n,n.r.ocp.x0Con,:x0_con)
  newConstraint!(n,n.r.ocp.xfCon,:xf_con)
  newConstraint!(n,n.r.ocp.dynCon,:dyn_con)

  # save the current working directory for navigation purposes
  n.r.mainDir = pwd()
  return nothing
end

"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 11/10/2017 \n
-------------------------------------------------------------------------------------\n
"""
function configure!(n::NLOpt; kwargs... )
  kw = Dict(kwargs)

  # final time
  if !haskey(kw,:finalTimeDV);n.s.ocp.finalTimeDV=false;
  else; n.s.ocp.finalTimeDV=get(kw,:finalTimeDV,0);
  end

  if !haskey(kw,:tf) && !n.s.ocp.finalTimeDV
    error("\n If the final is not a design variable pass it as: (:tf=>Float64(some #)) \n
        If the final time is a design variable, indicate that as: (:finalTimeDV=>true)\n")
  elseif haskey(kw,:tf) && !n.s.ocp.finalTimeDV
    n.ocp.tf = get(kw,:tf,0)
  elseif n.s.ocp.finalTimeDV
    n.ocp.tf = Any
  end

  # integrationScheme
  if !haskey(kw,:integrationScheme); n.s.ocp.integrationScheme=:lgrExplicit; # default
  else; n.s.ocp.integrationScheme=get(kw,:integrationScheme,0);
  end

  if n.s.ocp.integrationScheme==:lgrExplicit ||  n.s.ocp.integrationScheme==:lgrImplicit
    n.s.ocp.integrationMethod = :ps
  elseif n.s.ocp.integrationScheme==:trapezoidal || n.s.ocp.integrationScheme==:bkwEuler
    n.s.ocp.integrationMethod = :tm
  else
    error("the :integrationScheme that you specified is not currently implemeted \n")
  end

  if n.s.ocp.integrationMethod==:ps
    if haskey(kw,:N)
      error(" \n N is not an appropriate kwargs for :ps methods \n")
    end
    if !haskey(kw,:Nck);n.ocp.Nck=[10,10,10,10]; # default
    else; n.ocp.Nck = get(kw,:Nck,0);
    end
    n.ocp.Ni = length(n.ocp.Nck)

    for int in 1:n.ocp.Ni
        if (n.ocp.Nck[int]<0)
            error("\n Nck must be > 0");
        end
    end
    n.ocp.state.pts = sum(n.ocp.Nck) + 1
    n.ocp.control.pts = sum(n.ocp.Nck)
    n.ocp.Nck_full = [0;cumsum(n.ocp.Nck+1)]
    n.ocp.Nck_cum = [0;cumsum(n.ocp.Nck)]

    # initialize node data
    if n.s.ocp.integrationScheme==:lgrExplicit ||  n.s.ocp.integrationScheme==:lgrImplicit
      taus_and_weights = [gaussradau(n.ocp.Nck[int]) for int in 1:n.ocp.Ni];
    end
    n.ocp.tau = [taus_and_weights[int][1] for int in 1:n.ocp.Ni]
    n.ocp.w = [taus_and_weights[int][2] for int in 1:n.ocp.Ni]
    createIntervals!(n)
    DMatrix!(n)

  elseif n.s.ocp.integrationMethod==:tm
    if haskey(kw,:Nck)
      error(" \n Nck is not appropriate kwargs for :tm methods \n")
    end
    if !haskey(kw,:N);n.ocp.N=100; # default
    else; n.ocp.N = get(kw,:N,0);
    end
    n.ocp.state.pts = n.ocp.N + 1
    n.ocp.control.pts = n.ocp.N + 1
  end
  n.ocp.mXL = falses(n.ocp.state.num)
  n.ocp.mXU = falses(n.ocp.state.num)
  n.ocp.XL_var = Matrix{Float64}(n.ocp.state.num,n.ocp.state.pts)
  n.ocp.XU_var = Matrix{Float64}(n.ocp.state.num,n.ocp.state.pts)

  # solver settings
  if !haskey(kw,:solverSettings);SS = Dict((:name=>:Ipopt)); # default
  else; SS = get(kw,:solverSettings,0); SS = Dict(SS);
  end
  defineSolver!(n,SS)

  # optimal control problem
  OCPdef!(n)

  if n.s.ocp.evalCostates
     if n.s.ocp.integrationMethod != :ps
       error("costates are only implmented for :ps methods")
     end
     n.s.ocp.evalConstraints = true
  end
  return nothing
end
