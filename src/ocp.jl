"""
n,r=OCPdef(mdl,n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 1/30/2017 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef(mdl::JuMP.Model,n::NLOpt, args...)
  if length(args)==1; params=args[1]; paramsON=true; else paramsON=false; end

  # state and control variables
  @variable(mdl, x[1:n.numStatePoints,1:n.numStates]); # +1 for the last interval
  for st in 1:n.numStates
    if !isnan(n.XL[st])
      for j in 1:n.numStatePoints
        setlowerbound(x[j,st], n.XL[st])
      end
    end
    if !isnan(n.XU[st])
      for j in 1:n.numStatePoints
        setupperbound(x[j,st], n.XU[st])
      end
    end
  end
  @variable(mdl, u[1:n.numControlPoints,1:n.numControls]);  #TODO make sure that there are not too many controls
  for ctr in 1:n.numControls
    if !isnan(n.CL[ctr])
      for j in 1:n.numControlPoints
        setlowerbound(u[j,ctr], n.CL[ctr])
      end
    end
    if !isnan(n.CU[ctr])
      for j in 1:n.numControlPoints
        setupperbound(u[j,ctr], n.CU[ctr])
      end
    end
  end

  # boundary constraints
  x0_con=[];xf_con=[];
  for st in 1:n.numStates
    if !isnan(n.X0[st])
       x0_con = [x0_con; @constraint(mdl, x[1,st]==n.X0[st])]
    end
    if !isnan(n.XF[st])
      xf_con = [xf_con; @constraint(mdl, x[end,st]==n.XF[st])];
    end
  end

  if n.integrationMethod==:ps
    dyn_con  = [Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];
    Nck_st  = [0;cumsum(n.Nck+1)]; Nck_ctr = [0;cumsum(n.Nck)];
    dynamics_expr  = [Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];

    if n.finalTimeDV
      @variable(mdl, 0.1 <= tf <=  n.tf_max)
      n.tf = tf;
    end

    for int in 1:n.Ni
      # states
      x_int = Array(Any,length(Nck_st[int]+1:Nck_st[int+1]),n.numStates);
      for st in 1:n.numStates # NOTE we use Nck_ctr for state indexing
        x_int[:,st] = x[Nck_ctr[int]+1:Nck_ctr[int+1]+1,st];  #+1 adds the DV in the next interval
      end
      # controls
      u_int = Array(Any,length(Nck_ctr[int]+1:Nck_ctr[int+1]),n.numControls);
      for ctr in 1:n.numControls
        u_int[:,ctr] = u[Nck_ctr[int]+1:Nck_ctr[int+1],ctr];
      end
      # dynamics
      if paramsON
        dx = n.stateEquations(mdl,n,x_int,u_int,params);
      else
        dx = n.stateEquations(mdl,n,x_int,u_int);
      end
      for st in 1:n.numStates
        if n.integrationScheme==:lgrExplicit
          dynamics_expr[int][:,st] = @NLexpression(mdl, [j in 1:n.Nck[int]], sum(n.DMatrix[int][j,i]*x_int[i,st] for i in 1:n.Nck[int]+1) - ((n.tf-n.t0)/2)*dx[j,st]  )
        end
        for j in 1:n.Nck[int]
          dyn_con[int][j,st] = @NLconstraint(mdl, 0 == dynamics_expr[int][j,st])
        end
      end

    end
  elseif n.integrationMethod==:tm
    dyn_con  = Array(Any,n.N,n.numStates);
    if n.finalTimeDV
     #@variable( mdl, 0.00001 <= dt[1:n.N] <= 0.2) #TODO allow for an varaible array of dts
     @variable(mdl, 0.01 <= tf <= n.tf_max)
     n.tf = tf;
    end
    n.dt = n.tf/n.N*ones(n.N,);
    if paramsON
      dx = n.stateEquations(mdl,n,x,u,params);
    else
      dx = n.stateEquations(mdl,n,x,u);
    end
    if n.integrationScheme==:bkwEuler
      for st in 1:n.numStates
        dyn_con[:,st] = @NLconstraint(mdl, [j in 1:n.N], 0 == x[j+1,st] - x[j,st] - dx[j+1,st]*n.tf/(n.N) );
      end
    elseif n.integrationScheme==:trapezoidal
      for st in 1:n.numStates
        dyn_con[:,st] = @NLconstraint(mdl, [j in 1:n.N], 0 == x[j+1,st] - x[j,st] - 0.5*(dx[j,st] + dx[j+1,st])*n.tf/(n.N) )
      end
    end
  end
  # store results
  r = Result();
  r.x=x; r.u=u; r.x0_con=x0_con; r.xf_con=xf_con;
  return n, r
end
