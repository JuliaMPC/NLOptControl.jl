"""
OCPdef!(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 5/28/2017 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef!(n::NLOpt,args...)

  if length(args)==1; params=args[1]; paramsON=true; else paramsON=false; end
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
      @variable(n.mdl, 0.001 <= tf <=  n.tf_max)
      n.tf=tf;
      create_tV!(n) # make a time vector
    end

    for int in 1:n.Ni
      # states
      x_int=Array(Any,length(Nck_full[int]+1:Nck_full[int+1]),n.numStates);
      for st in 1:n.numStates # +1 adds the DV in the next interval
        x_int[:,st]=n.r.x[Nck_vars[int]+1:Nck_vars[int+1]+1,st];
      end

      # controls
      if int!=n.Ni
        u_int=Array(Any,length(Nck_full[int]+1:Nck_full[int+1]),n.numControls);
      else # -1 -> removing control in last mesh interval
        u_int=Array(Any,length(Nck_full[int]+1:Nck_full[int+1]-1),n.numControls);
      end

      for ctr in 1:n.numControls
        if int!=n.Ni # +1 adds the DV in the next interval
          u_int[:,ctr]=n.r.u[Nck_vars[int]+1:Nck_vars[int+1]+1,ctr];
        else
          u_int[:,ctr]=n.r.u[Nck_vars[int]+1:Nck_vars[int+1],ctr];
        end
      end
      # dynamics
      if paramsON
        dx=n.stateEquations(n,x_int,u_int,params);
      else
        dx=n.stateEquations(n,x_int,u_int);
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
    if n.finalTimeDV
     #@variable( mdl, 0.00001 <= dt[1:n.N] <= 0.2) #TODO allow for an varaible array of dts
     @variable(n.mdl, 0.001 <= tf <= n.tf_max)
     n.tf = tf;
    end
    n.dt = n.tf/n.N*ones(n.N,);
    if paramsON
      dx = n.stateEquations(n,n.r.x,n.r.u,params);
    else
      dx = n.stateEquations(n,n.r.x,n.r.u); # TODO for now leaving the reduntant terms so that the diff eq functions are consistent
    end
    if n.s.integrationScheme==:bkwEuler
      for st in 1:n.numStates
        n.r.dyn_con[:,st] = @NLconstraint(n.mdl, [j in 1:n.N], 0 == n.r.x[j+1,st]-n.r.x[j,st] - dx[j+1,st]*n.tf/(n.N) );
      end
    elseif n.s.integrationScheme==:trapezoidal
      for st in 1:n.numStates
        n.r.dyn_con[:,st] = @NLconstraint(n.mdl, [j in 1:n.N], 0 == n.r.x[j+1,st]-n.r.x[j,st] - 0.5*(dx[j,st] + dx[j+1,st])*n.tf/(n.N) )
      end
    end
  end

  # save constraint data
  newConstraint!(n.r,n.r.x0_con,:x0_con);
  newConstraint!(n.r,n.r.xf_con,:xf_con);
  newConstraint!(n.r,n.r.dyn_con,:dyn_con);

  # save the current working directory for navigation purposes
  n.r.main_dir=pwd();
  nothing
end
