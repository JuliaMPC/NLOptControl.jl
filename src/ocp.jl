
function OCPdef(mdl::JuMP.Model, nlp::NLP_data, ps::PS_data)
  # inequality constraints and design variable definitions
  @unpack numStates, numStatePoints,  XL, XU = nlp
  @variable(mdl, x[1:sum(numStatePoints),1:numStates]);
  for st in 1:numStates
    for j in 1:sum(numStatePoints)
      setlowerbound(x[j,st], XL[st])
      setupperbound(x[j,st], XU[st])
    end
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

  # state continuity constraints
  @unpack Ni, Nck = ps;
  Nck_st = [1;cumsum(Nck+1)];
  Nck_ctr = [1;cumsum(Nck)];
  for int in 2:Ni
    for st in 1:numStates
      @constraint(mdl, x[Nck_st[int],st] == x[Nck_st[int]+1,st])
    end
    for ctr in 1:numControls #TODO why do we need this?
      @constraint(mdl, u[Nck_ctr[int],ctr] == u[Nck_ctr[int]+1,ctr])
    end
  end

  # calculate LGR matrices - > IMatrix and DMatrix
  LGR_matrices(ps,nlp); # TODO if the final time is changing, DMatrix will change!--> needs to be recalcualted during optimization
  @unpack DMatrix, t0, tf = ps; # TODO look into regestering the DMatrix

  # state equation constraints
  Nck_st  = [0;cumsum(Nck+1)];
  Nck_ctr = [0;cumsum(Nck)];
  for int in 1:Ni
    # states
    x_int = Array(Any,length(Nck_st[int]+1:Nck_st[int+1]),numStates);
    for st in 1:numStates
      x_int[:,st] = x[Nck_st[int]+1:Nck_st[int+1],st];
    end
    # controls
    u_int = Array(Any,length(Nck_ctr[int]+1:Nck_ctr[int+1]),numControls);
    for ctr in 1:numControls
      u_int[:,ctr] = u[Nck_ctr[int]+1:Nck_ctr[int+1],ctr];
    end
    # dynamics
    @unpack stateEquations = nlp
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

 return x,u
end
