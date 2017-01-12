"""
c, ceq = constraints(nlp,ps)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/8/2017, Last Modified: 1/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function constraints(nlp::NLP_data,ps::PS_data)

  # calculate equality constraints
  ode_constraint(nlp,ps); # TODO consider running these in parallel
  continuity_constraint(nlp,ps);
  boundary_constraint(nlp,ps);
  @unpack odeConstraint, continuityConstraint, boundaryConstraint = ps
  ceq1 = [idx for tempM in odeConstraint for idx = tempM];
  ceq2 = [idx for tempM in continuityConstraint for idx = tempM];
  ceq3 = [idx for tempM in boundaryConstraint for idx = tempM];
  ceq = [ceq1; ceq2; ceq3];

  # calculate inequality constraints
  inequality_constraint(nlp,ps);
  @unpack stateConstraint, controlConstraint = ps
  c1 = [idx for tempM in stateConstraint for idx = tempM];
  c2 = [idx for tempM in controlConstraint for idx = tempM];
  c = [c1; c2];
  return c, ceq
end

"""
ode_constraint(nlp,ps)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/8/2017, Last Modified: 1/12/2017 \n
--------------------------------------------------------------------------------------\n
"""
function ode_constraint(nlp::NLP_data,ps::PS_data)
    @unpack FMatrix, DMatrix, stateMatrix, Ni, odeConstraint = ps
    FMatrix = F_matrix(nlp,ps);
    for int in 1:Ni
        odeConstraint[int] = DMatrix[int]*stateMatrix[int]-FMatrix[int];
    end
    @pack ps = odeConstraint, FMatrix
end

"""
continuity_constraint(nlp,ps)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/7/2017, Last Modified: 1/7/2017 \n
--------------------------------------------------------------------------------------\n
"""
function continuity_constraint(nlp::NLP_data,ps::PS_data)
    @unpack stateMatrix, controlMatrix, Ni, continuityConstraint = ps
    @unpack numStates, numControls = nlp
    for int in 1:Ni-1
        continuityConstraint[int] = [stateMatrix[int][end,:] - stateMatrix[int+1][1,:];
                                   controlMatrix[int][end,:] - controlMatrix[int+1][1,:]];
    end
    @pack ps = continuityConstraint
end

"""
boundary_constraint(nlp,ps)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/7/2017, Last Modified: 1/8/2017 \n
--------------------------------------------------------------------------------------\n
"""
function boundary_constraint(nlp::NLP_data,ps::PS_data)
    @unpack stateMatrix, Ni, boundaryConstraint = ps
    @unpack X0, XF = nlp
    boundaryConstraint = [X0 - stateMatrix[1][1,:];
                          XF - stateMatrix[Ni][end,:]];
    @pack ps = boundaryConstraint
end


"""
inequality_constraint(nlp,ps)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/8/2017, Last Modified: 1/8/2017 \n
--------------------------------------------------------------------------------------\n
"""
function inequality_constraint(nlp::NLP_data,ps::PS_data)
    @unpack Ni, stateMatrix, stateConstraint = ps
    @unpack numStates, XL, XU = nlp
    for st in 1:numStates
        for int in 1:Ni;
            L = [idx for tempM in stateMatrix for idx = XL[st] - tempM[:,st]];
            U = [idx for tempM in stateMatrix for idx = tempM[:,st] - XU[st]];
            stateConstraint[st]= [L U];
        end
    end

    @unpack controlMatrix, controlConstraint = ps
    @unpack numControls, CL, CU = nlp
    for ctr in 1:numControls
        for int in 1:Ni;
            L = [idx for tempM in controlMatrix for idx = CL[ctr] - tempM[:,ctr]];
            U = [idx for tempM in controlMatrix for idx = tempM[:,ctr] - CU[ctr]];
            controlConstraint[ctr]= [L U];
        end
    end
    @pack ps = stateConstraint, controlConstraint
end
