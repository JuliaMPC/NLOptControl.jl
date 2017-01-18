
"""
OCPdef(mdl,nlp,ps)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/14/2017, Last Modified: 1/17/2017 \n
--------------------------------------------------------------------------------------\n
"""
function OCPdef(mdl::JuMP.Model, nlp::NLP_data, ps::PS_data)
  # inequality constraints and design variable definitions
  @unpack numStates, numStatePoints,  XL, XU = nlp  #TODO with current formulation numStatePoints == numControlPoints consider simplifying
  @variable(mdl, x[1:sum(numStatePoints-1)+1,1:numStates]); # +1 for the last interval
  for st in 1:numStates
    for j in 1:sum(numStatePoints-1)
      setlowerbound(x[j,st], XL[st])
      setupperbound(x[j,st], XU[st])
    end
  end

  @unpack finalTimeDV = nlp
  if finalTimeDV
    @variable(mdl, tf_var)
    @unpack Ni, Nck, τ, ω = ps
    const Ni_const = Ni  #TODO either way these can all be constaints ! --> simplify
    const Nck_const = Nck
    const τ_const = τ
    const ω_const = ω
    #JuMP.register(mdl, :create_intervals_JuMP, 1, create_intervals_JuMP, autodiff=true)
    #JuMP.register(mdl, :poldif_JuMP, 1, poldif_JuMP, autodiff=true)
    #JuMP.register(mdl, :D_matrix_JuMP, 1, D_matrix_JuMP, autodiff=true)
    DMatrix_JuMP = D_matrix_JuMP(mdl,tf_var,Nck_const,Ni_const,τ_const,ω_const);
    state_eqs = [Array(Any,Nck[int],numStates) for int in 1:Ni];
    dynamics_expr  = [Array(Any,Nck[int],numStates) for int in 1:Ni];

  end

  @unpack numControls, numControlPoints, CL, CU = nlp
  @variable(mdl, u[1:sum(numControlPoints),1:numControls]);
  for ctr in 1:numControls
    for j in 1:sum(numControlPoints)
      setlowerbound(u[j,ctr], CL[ctr])
      setupperbound(u[j,ctr], CU[ctr])
    end
  end

  # boundary constraints
  @unpack X0, XF = nlp
  for st in 1:numStates
    @constraint(mdl, x[1,st]   == X0[st]);
    @constraint(mdl, x[end,st] == XF[st]);
  end

  # state equation constraints
  @unpack Ni, Nck = ps;
  Nck_st  = [0;cumsum(Nck+1)];
  Nck_ctr = [0;cumsum(Nck)];

  ## TEMP
  const g = 1.62519; # m/s^2
  state_eqs = [Array(Any,Nck[int],numStates) for int in 1:Ni];
  for int in 1:Ni
    state_eqs[int][:,1] = @NLexpression(mdl, [j=1:Nck[int]], x[Nck_ctr[int]+1:Nck_ctr[int+1]+1,2][j]);
    state_eqs[int][:,2] = @NLexpression(mdl, [j=1:Nck[int]], u[Nck_ctr[int]+1:Nck_ctr[int+1],1][j]-g);
  end
  ## TEMP - eventually push back up or initialize based off of user
  for int in 1:Ni
    # states
    x_int = Array(Any,length(Nck_st[int]+1:Nck_st[int+1]),numStates);
    for st in 1:numStates # NOTE we use Nck_ctr for state indexing
      x_int[:,st] = x[Nck_ctr[int]+1:Nck_ctr[int+1]+1,st];  #+1 adds the DV in the next interval
    end
    # controls
    u_int = Array(Any,length(Nck_ctr[int]+1:Nck_ctr[int+1]),numControls);
    for ctr in 1:numControls
      u_int[:,ctr] = u[Nck_ctr[int]+1:Nck_ctr[int+1],ctr];
    end
    # dynamics
    if finalTimeDV
      for st in 1:numStates
      #  SE = stateEquations(x_int,u_int,st);  # need to generate @NLexpression() to work in the following
        dynamics_expr[int][:,st] = @NLexpression(mdl, [j in 1:Nck[int]], sum(DMatrix_JuMP[int][j,i]*x_int[i,st] - state_eqs[int][i,st]  for i in 1:Nck[int]))
        for j in 1:Nck[int]
          @NLconstraint(mdl, 0 == dynamics_expr[int][j,st])
        end
      end
    else
      # calculate LGR matrices - > IMatrix and DMatrix TODO make option to use IMatrix or DMatrix
      @unpack stateEquations = nlp
      LGR_matrices(ps,nlp);
      @unpack DMatrix = ps;
      dynamics = Array(Any,length(Nck_ctr[int]+1:Nck_ctr[int+1]),numStates);
      for st in 1:numStates
        dynamics[:,st] = DMatrix[int]*x_int[:,st] - stateEquations(x_int,u_int,st)
      end
      # add dynamics constraints
      for st in 1:numStates
        for j in 1:Nck[int]
          @constraint(mdl, 0 == dynamics[j,st])
        end
      end

    end

  end

 return x,u
end
