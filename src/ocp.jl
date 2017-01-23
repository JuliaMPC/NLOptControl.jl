"""
n,x,u,c=OCPdef(mdl,n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 1/23/2017 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef(mdl::JuMP.Model,n::NLOpt)

  if n.integrationMethod==:ps
    @variable(mdl, x[1:sum(n.numPoints)+1,1:n.numStates]); # +1 for the last interval
    for st in 1:n.numStates
      for j in 1:sum(n.numPoints-1)
        setlowerbound(x[j,st], n.XL[st])
        setupperbound(x[j,st], n.XU[st])
      end
    end
    @variable(mdl, u[1:sum(n.numPoints),1:n.numControls]);  #TODO make sure that there are not too many controls
    for ctr in 1:n.numControls
      for j in 1:sum(n.numPoints)
        setlowerbound(u[j,ctr], n.CL[ctr])
        setupperbound(u[j,ctr], n.CU[ctr])
      end
    end
  elseif n.integrationMethod==:tm
    @variable(mdl, x[1:n.N+1,1:n.numStates]); # +1 for the last interval
    for st in 1:n.numStates
      for j in 1:n.N+1
        setlowerbound(x[j,st], n.XL[st])
        setupperbound(x[j,st], n.XU[st])
      end
    end
    @variable(mdl, u[1:n.N+1,1:n.numControls]);  #TODO make sure that there are not too many controls
    for ctr in 1:n.numControls
      for j in 1:n.N+1
        setlowerbound(u[j,ctr], n.CL[ctr])
        setupperbound(u[j,ctr], n.CU[ctr])
      end
    end
  end
  # boundary constraints
  x0_con = @constraint(mdl, [st in 1:n.numStates], x[1,st]  == n.X0[st]);
  xf_con = @constraint(mdl, [st in 1:n.numStates], x[end,st] == n.XF[st]);

  if n.integrationMethod==:ps
    dyn_con  = [Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];
    Nck_st  = [0;cumsum(n.Nck+1)]; Nck_ctr = [0;cumsum(n.Nck)];
    dx = [Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];
    if n.finalTimeDV
      dynamics_expr  = [Array(Any,n.Nck[int],n.numStates) for int in 1:n.Ni];
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
      if n.finalTimeDV
        @variable(mdl, 0.01 <= n.tf <= Inf) #TODO allow user to pass ranges
        n.ts, n.ωₛ = create_intervals(mdl,n);
        n.DMatrix = D_matrix(mdl,n);
        for int in 1:n.Ni
          for st in 1:n.numStates
            dx[int][:,st] = @NLexpression(mdl, [j=1:n.Nck[int]], n.stateEquations([Nck_ctr[int]+1:Nck_ctr[int+1]+1,:][j],u[Nck_ctr[int]+1:Nck_ctr[int+1],:][j],st) );
          end
        end
        for st in 1:n.numStates
          if integrationScheme==:lgrExplicit
            dynamics_expr[int][:,st] = @NLexpression(mdl, [j in 1:n.Nck[int]], sum(DMatrix[int][j,i]*x_int[i,st] -dx[int][i,st]  for i in 1:n.Nck[int]))
          end
          for j in 1:n.Nck[int]
            dyn_con[int][j,st] = @NLconstraint(mdl, 0 == dynamics_expr[int][j,st])
          end
        end
      else
        # redo these every interval
        n=LGR_matrices(n);       # calculate LGR matrices
        dx = [n.stateEquations(x_int,u_int,st) for st in 1:n.numStates];
        dynamics = Array(Any,length(Nck_ctr[int]+1:Nck_ctr[int+1]),n.numStates);
        for st in 1:n.numStates
          dynamics[:,st] = n.DMatrix[int]*x_int[:,st] - dx[st];
          for j in 1:n.Nck[int] # add dynamics constraints
            dyn_con[int][j,st] = @constraint(mdl, 0 == dynamics[j,st])
          end
        end
      end

    end
  elseif n.integrationMethod==:tm
    dyn_con  = Array(Any,n.N,n.numStates);
    dx = Array(Any,n.N+1,n.numStates)

    if n.finalTimeDV
     #@variable( mdl, 0.00001 <= dt[1:n.N] <= 0.2) #TODO allow for an varaible array of dts
     @variable(mdl, 0.01 <= n.tf <= Inf) #TODO allow user to pass ranges
    end
    n.dt = n.tf/n.N*ones(n.N,);

    for st in 1:n.numStates
      dx[:,st] = @NLexpression(mdl, [j=1:n.N+1], n.stateEquations(x,u,st)[j] );
    end

    if n.integrationScheme==:bkwEuler
      for st in 1:n.numStates
        dyn_con[:,st] = @NLconstraint(mdl, [j in 1:n.N], 0 == x[j+1,st] - x[j,st] - dx[j+1,st]*n.dt[j] )
      end
    elseif n.integrationScheme==:trapezoidal
      for st in 1:n.numStates
        dyn_con[:,st] = @NLconstraint(mdl, [j in 1:n.N], 0 == x[j+1,st] - x[j,st] - (dx[j,st] + dx[j+1,st])*n.dt[j]/2 )
      end
    end
  end
  c = [x0_con, xf_con, dyn_con];  # pack up all constraints
  return n,x,u,c
end

#= old
#JuMP.register(mdl, :create_intervals_JuMP, 1, create_intervals_JuMP, autodiff=true)
#JuMP.register(mdl, :poldif_JuMP, 1, poldif_JuMP, autodiff=true)
#JuMP.register(mdl, :D_matrix_JuMP, 1, D_matrix_JuMP, autodiff=true)
=#
