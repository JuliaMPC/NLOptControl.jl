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
    numStates::Int = 0,
    numControls::Int = 0,
    X0 = fill(NaN,numStates,),
    XF = fill(NaN,numStates,),
    XL = fill(NaN,numStates,),
    XU = fill(NaN,numStates,),
    CL = fill(NaN,numControls,),
    CU = fill(NaN,numControls,),
    XS = fill(1.,numStates,),
    CS = fill(1.,numControls,)
)::NLOpt

    n = NLOpt()

    # validate input
    if numControls <= 0
        @error "numControls ($numControls) must be > 0","\n"
    end
    if numStates <= 0
        @error "numStates ($numStates) must be > 0","\n"
    end
    if length(X0) != numStates
        @error "\n Length of X0 ($(length(X0))) must match number of states ($numStates)\n"
    end
    if length(XF) != numStates
        @error "\n Length of XF ($(length(XF))) must match number of states ($numStates)\n"
    end
    if length(XL) != numStates
        @error "\n Length of XL ($(length(XL))) must match number of states ($numStates)\n"
    end
    if length(XU) != numStates
        @error "\n Length of XU ($(length(XU))) must match number of states ($numStates)\n"
    end
    if length(XS) != numStates
        @error "\n Length of XS ($(length(XS))) must match number of states ($numStates)\n"
    end
    if length(CL) != numControls
        @error "\n Length of CL ($(length(CL))) must match number of controls ($numControls)\n"
    end
    if length(CU) != numControls
        @error "\n Length of CU ($(length(CU))) must match number of controls ($numControls)\n"
    end
    if length(CS) != numControls
        @error "\n Length of CS ($(length(CS))) must match number of controls ($numControls)\n"
    end

    n.ocp.state = initState(numStates)
    n.ocp.control = initControl(numControls)
    n.ocp.X0 = X0
    n.ocp.X0_tol = fill(NaN, size(X0))
    n.ocp.XF = XF
    n.ocp.XF_tol = fill(NaN, size(XF))
    n.ocp.XL = XL
    n.ocp.XU = XU
    n.ocp.CL = CL
    n.ocp.CU = CU
    n.ocp.XS = XS
    n.ocp.CS = CS
    n.f.ocp.defined = true

    return n

end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/9/2017, Last Modified: 12/02/2019 \n
Last Modified by: Alexander Buck
-------------------------------------------------------------------------------------\n
"""
function defineSolver!(n::NLOpt, kw::Dict{Symbol,Symbol})

    # Get the name of the solver
    if haskey(kw,:name)
        n.s.ocp.solver.name=get(kw,:name,0)
    end

    # Import solver into Julia environment
    try_import(n.s.ocp.solver.name) ? nothing : @error "Could not import $(n.s.ocp.solver.name) into Julia environment"

    # Check if user provided MPC defaults
    if haskey(kw, :mpc_defaults)
        mpc_defaults = get(kw,:mpc_defaults,0)
    else
        mpc_defaults=false
    end

    # Set solver defaults
    if mpc_defaults
        if n.s.ocp.solver.name == :Ipopt
            n.s.ocp.solver.settings = NLOptControl.NLOptBase._Ipopt_MPC
        else
            @error "Solver $(n.s.ocp.solver.name) not defined"
        end
    else # default solver settings
        if n.s.ocp.solver.name == :Ipopt  # NOTE this should already have been done by default, but could get messed up is user is playing with options
            n.s.ocp.solver.settings = NLOptControl.NLOptBase._Ipopt_defaults;
        else
            @error "Solver $(n.s.ocp.solver.name) not defined"
        end
    end

    # Modify additional defaults individually
    for (key,value) in kw
        if haskey(n.s.ocp.solver.settings, key)
            n.s.ocp.solver.settings[key] = value
        elseif key != :name && key != :mpc_defaults # ignore the name and default settings option # TODO could remove them from the Dict
            @error " \n Unknown key: $kw for $(n.s.ocp.solver.name) used in defineSolver!() \n "
        end
    end

    # Set solver
    if n.s.ocp.solver.name == :Ipopt
        setsolver(
            n.ocp.mdl,
            Ipopt.IpoptSolver(;
                max_cpu_time               = n.s.ocp.solver.settings[:max_cpu_time],
                print_level                = n.s.ocp.solver.settings[:print_level],
                warm_start_init_point      = n.s.ocp.solver.settings[:warm_start_init_point],
                max_iter                   = n.s.ocp.solver.settings[:max_iter],
                tol                        = n.s.ocp.solver.settings[:tol],
                dual_inf_tol               = n.s.ocp.solver.settings[:dual_inf_tol],
                constr_viol_tol            = n.s.ocp.solver.settings[:constr_viol_tol],
                compl_inf_tol              = n.s.ocp.solver.settings[:compl_inf_tol],
                acceptable_tol             = n.s.ocp.solver.settings[:acceptable_tol],
                acceptable_constr_viol_tol = n.s.ocp.solver.settings[:acceptable_constr_viol_tol],
                acceptable_dual_inf_tol    = n.s.ocp.solver.settings[:acceptable_dual_inf_tol],
                acceptable_compl_inf_tol   = n.s.ocp.solver.settings[:acceptable_compl_inf_tol],
                acceptable_obj_change_tol  = n.s.ocp.solver.settings[:acceptable_obj_change_tol],
                diverging_iterates_tol     = n.s.ocp.solver.settings[:diverging_iterates_tol]
            )
        )
    else
        @error "Solver $(n.s.ocp.solver.name) not defined"
    end

    return nothing

end  # function

"""
OCPdef!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 12/06/2019 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef!(n::NLOpt)

    # State variables
    @variable(n.ocp.mdl, x[1:n.ocp.state.pts, 1:n.ocp.state.num])
    n.r.ocp.x = Matrix{Any}(undef, n.ocp.state.pts, n.ocp.state.num)
    for st in 1:n.ocp.state.num
        n.r.ocp.x[:,st] = @NLexpression(n.ocp.mdl, [j=1:n.ocp.state.pts], n.ocp.XS[st]*x[j,st])
    end
    #  n.r.ocp.x=x;

    for st in 1:n.ocp.state.num
        # lower state constraint
        if !isnan.(n.ocp.XL[st])
        if n.ocp.mXL[st] == false
            for j in 1:n.ocp.state.pts
            setlowerbound(x[j,st], n.ocp.XL[st]/n.ocp.XS[st])
            end
        else
            for j in 1:n.ocp.state.pts
            setlowerbound(x[j,st], n.ocp.XL_var[st,j]/n.ocp.XS[st])
            end
        end
        end

        # upper state constraint
        if !isnan.(n.ocp.XU[st])
        if n.ocp.XU[st] != false
            for j in 1:n.ocp.state.pts
            setupperbound(x[j,st], n.ocp.XU[st]/n.ocp.XS[st])
            end
        else
            for j in 1:n.ocp.state.pts
            setlowerbound(x[j,st], n.ocp.XU_var[st,j]/n.ocp.XS[st])
            end
        end
        end
    end

    # control variables
    @variable(n.ocp.mdl,u[1:n.ocp.control.pts,1:n.ocp.control.num])
    n.r.ocp.u = Matrix{Any}(undef, n.ocp.control.pts,n.ocp.control.num)
    for ctr in 1:n.ocp.control.num
        n.r.ocp.u[:,ctr] = @NLexpression(n.ocp.mdl, [j=1:n.ocp.control.pts], n.ocp.CS[ctr]*u[j,ctr])
    end
    #n.r.ocp.u=u;

    for ctr in 1:n.ocp.control.num
        if !isnan.(n.ocp.CL[ctr])
        for j in 1:n.ocp.control.pts
            setlowerbound(u[j,ctr], n.ocp.CL[ctr]/n.ocp.CS[ctr])
        end
        end
        if !isnan.(n.ocp.CU[ctr])
        for j in 1:n.ocp.control.pts
            setupperbound(u[j,ctr], n.ocp.CU[ctr]/n.ocp.CS[ctr])
        end
        end
    end

    # boundary constraints
    if n.s.ocp.x0slackVariables   # create handles for constraining the entire initial state
        n.r.ocp.x0Con = Array{Any}(undef, n.ocp.state.num,2) # this is so they can be easily reference when doing MPC
    else
        n.r.ocp.x0Con = []
    end

    if n.s.ocp.xFslackVariables # create handles for constraining the entire final state
        n.r.ocp.xfCon = Array{Any}(undef, n.ocp.state.num,2) # this is so they can be easily reference when doing MPC
    else
        n.r.ocp.xfCon = []
    end

    if n.s.ocp.x0slackVariables # currently adding slack variables for all initial states even if there are no constraints on them
        @variable(n.ocp.mdl,x0s[1:n.ocp.state.num])
        n.ocp.x0s = x0s
        n.r.ocp.x0sCon = Array{Any}(undef, n.ocp.state.num)
    end
    if n.s.ocp.xFslackVariables # currently adding slack variables for all final states even if there are no constraints on them
        @variable(n.ocp.mdl,xFs[1:n.ocp.state.num])
        n.ocp.xFs = xFs
        n.r.ocp.xfsCon = Array{Any}(undef, n.ocp.state.num)
    end

    for st in 1:n.ocp.state.num
        if !isnan(n.ocp.X0[st]) # could have a bool for this
        if n.s.ocp.x0slackVariables
            n.r.ocp.x0Con[st,1] = @constraint(n.ocp.mdl, x[1,st]*n.ocp.XS[st] -n.ocp.x0s[st] <=  n.ocp.X0[st])
            n.r.ocp.x0Con[st,2] = @constraint(n.ocp.mdl,-x[1,st]*n.ocp.XS[st] -n.ocp.x0s[st] <= -n.ocp.X0[st])
            if !isnan(n.ocp.X0_tol[st])
            n.r.ocp.x0sCon[st] = @constraint(n.ocp.mdl, n.ocp.x0s[st] <= n.ocp.X0_tol[st])
            end
        else
            n.r.ocp.x0Con = [n.r.ocp.x0Con; @constraint(n.ocp.mdl, x[1,st]*n.ocp.XS[st] == n.ocp.X0[st])]
        end
        end
        if !isnan(n.ocp.XF[st])
        if n.s.ocp.xFslackVariables
            n.r.ocp.xfCon[st,1] = @constraint(n.ocp.mdl, x[end,st]*n.ocp.XS[st] -n.ocp.xFs[st] <=  n.ocp.XF[st] )
            n.r.ocp.xfCon[st,2] = @constraint(n.ocp.mdl,-x[end,st]*n.ocp.XS[st] -n.ocp.xFs[st] <= -n.ocp.XF[st] )
            if !isnan(n.ocp.XF_tol[st])
            n.r.ocp.xfsCon[st] = @constraint(n.ocp.mdl, n.ocp.xFs[st] <= n.ocp.XF_tol[st])
            end
        else
            n.r.ocp.xfCon = [n.r.ocp.xfCon; @constraint(n.ocp.mdl, x[end,st]*n.ocp.XS[st] == n.ocp.XF[st])]
        end
        end
    end

    @NLparameter(n.ocp.mdl,t0_param == 0.0);   # for now we just start at zero
    n.ocp.t0 = t0_param  # NOTE consider making a constraint that t0 < tf

    if n.s.ocp.integrationMethod == :ps
        n.r.ocp.dynCon = [Array{Any}(undef, n.ocp.Nck[int],n.ocp.state.num) for int in 1:n.ocp.Ni]
        dynamics_expr = [Array{Any}(undef, n.ocp.Nck[int],n.ocp.state.num) for int in 1:n.ocp.Ni]

        if n.s.ocp.finalTimeDV
        @variable(n.ocp.mdl, n.s.ocp.tfMin <= tf <=  n.s.ocp.tfMax)
        n.ocp.tf = tf
        create_tV!(n)          # make a time vector
        end

        for int in 1:n.ocp.Ni
        x_int,u_int = intervals(n,int,n.r.ocp.x,n.r.ocp.u)

        # dynamics
        L = size(x_int)[1]-1;
        dx = Matrix{Any}(undef, L,n.ocp.state.num)
        for st in 1:n.ocp.state.num
            dx[:,st] = DiffEq(n,x_int,u_int,L,st)
        end

        for st in 1:n.ocp.state.num # TODO consider multiplying X*D to reduce computations (i.e. remove this for loop for the states)
            if n.s.ocp.integrationScheme == :lgrExplicit
            dynamics_expr[int][:,st] = @NLexpression(n.ocp.mdl, [j in 1:n.ocp.Nck[int]], sum(n.ocp.DMatrix[int][j,i]*x_int[i,st] for i in 1:n.ocp.Nck[int]+1) - ((n.ocp.tf)/2)*dx[j,st]  )
            elseif n.s.ocp.integrationScheme==:lgrImplicit
            dynamics_expr[int][:,st] = @NLexpression(n.ocp.mdl, [j in 1:n.ocp.Nck[int]], x_int[j+1,st] - x_int[1,st] - ((n.ocp.tf)/2)*sum(n.ocp.IMatrix[int][j,i]*dx[i,st] for i in 1:n.ocp.Nck[int]) )
            end
            for j in 1:n.ocp.Nck[int]
            n.r.ocp.dynCon[int][j,st] = @NLconstraint(n.ocp.mdl, 0. == dynamics_expr[int][j,st])
            end
        end

        # additional constraints
        for num in 1:length(n.ocp.NLcon)
            ch = addCon(n,x_int,u_int,L,num)
            newConstraint!(n,ch,Symbol(string("ch",num))) # TODO could let the user name these
        end
        end
    elseif n.s.ocp.integrationMethod == :tm
        n.r.ocp.dynCon = Array{Any}(undef, n.ocp.N,n.ocp.state.num)
        if n.s.ocp.finalTimeDV
        @variable(n.ocp.mdl, n.s.ocp.tfMin <= tf <= n.s.ocp.tfMax)
        n.ocp.tf = tf
        create_tV!(n)          # make a time vector
        end
        n.ocp.dt = n.ocp.tf/n.ocp.N*ones(n.ocp.N,)

        L = size(n.r.ocp.x)[1]
        dx = Array{Any}(undef, L,n.ocp.state.num)
        for st in 1:n.ocp.state.num
        dx[:,st] = DiffEq(n,n.r.ocp.x,n.r.ocp.u,L,st)
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
Date Create: 1/1/2017, Last Modified: 12/06/2019 \n
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

  # x0 slack variables
  if !haskey(kw,:x0slackVariables)
    if !any(isnan.(n.ocp.X0_tol))
      n.s.ocp.x0slackVariables = true
    else
      n.s.ocp.x0slackVariables = false # default
    end
  else; n.s.ocp.x0slackVariables=get(kw,:x0slackVariables,0);  # TODO check to see if user is using tolerances and passed false to configure
  end

  # xF slack variables
  if !haskey(kw,:xFslackVariables)
    if !any(isnan.(n.ocp.XF_tol))
      n.s.ocp.xFslackVariables = true
    else
      n.s.ocp.xFslackVariables = false # default
    end
  else; n.s.ocp.xFslackVariables=get(kw,:xFslackVariables,0);
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
    if !haskey(kw,:Nck);n.ocp.Nck = [10,10,10,10]; # default
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
    n.ocp.Nck_full = [0;cumsum(n.ocp.Nck.+1)]
    n.ocp.Nck_cum = [0;cumsum(n.ocp.Nck)]

    # initialize node data
    if n.s.ocp.integrationScheme==:lgrExplicit ||  n.s.ocp.integrationScheme==:lgrImplicit
      taus_and_weights = [ gaussradau(n.ocp.Nck[int]) for int in 1:n.ocp.Ni ];
    end
    n.ocp.tau = [ taus_and_weights[int][1] for int in 1:n.ocp.Ni ]
    n.ocp.w = [ taus_and_weights[int][2] for int in 1:n.ocp.Ni ]
    createIntervals!(n)
    DMatrix!(n)

  elseif n.s.ocp.integrationMethod==:tm
    if haskey(kw,:Nck)
      error(" \n Nck is not appropriate kwargs for :tm methods \n")
    end
    if !haskey(kw,:N);n.ocp.N = 100; # default
    else; n.ocp.N = get(kw,:N,0);
    end
    n.ocp.state.pts = n.ocp.N + 1
    n.ocp.control.pts = n.ocp.N + 1
  end
  n.ocp.mXL = falses(n.ocp.state.num)
  n.ocp.mXU = falses(n.ocp.state.num)
  n.ocp.XL_var = Matrix{Float64}(undef, n.ocp.state.num, n.ocp.state.pts)
  n.ocp.XU_var = Matrix{Float64}(undef, n.ocp.state.num, n.ocp.state.pts)

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
