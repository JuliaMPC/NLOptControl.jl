"""
LGR_matrices(ps,nlp)
n=LGR_matrices(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/3/2017, Last Modified: 1/23/2017 \n
--------------------------------------------------------------------------------------\n
"""
function LGR_matrices(n::NLOpt)

    # check input
    check_ts = zeros(Float64,n.Ni);
    for int in 1:n.Ni
        check_ts[int]=maximum(n.ts[int]);
    end
    if maximum(check_ts) < 10*eps()
      error("\n ts is full of zeros: make sure that you call create_intervals() first to calculate ts! \n
                NOTE: This may have occured because  (:finalTimeDV => true) and the final time dv is not working properly!! \n")
    end

    D = [zeros((n.Nck[int]+1),(n.Nck[int]+1)) for int in 1:n.Ni];
    for int in 1:n.Ni
        D[int] = poldif(n.ts[int], 1) # +1 is already appended onto ts
    end

    n.DMatrix = [zeros((n.Nck[int]),(n.Nck[int]+1)) for int in 1:n.Ni];
    DM = [zeros((n.Nck[int]),(n.Nck[int])) for int in 1:n.Ni];
    for int in 1:n.Ni
        n.DMatrix[int] = D[int][1:end-1,:]; # [Nck]X[Nck+1]
        DM[int] = D[int][1:end-1,2:end];    # [Nck]X[Nck]
      #  IMatrix[int] = inv(DM[int]);        # I = inv[D_{2:N_k+1}]
    end
    return n
end
"""
n = D_matrix(mdl,n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/15/2017, Last Modified: 1/23/2017 \n
--------------------------------------------------------------------------------------\n
"""
function D_matrix(mdl::JuMP.Model,n::NLOpt)
  #D = [Matrix((Nck[int]+1),(Nck[int]+1)) for int in 1:Ni];
  D = polyDiff(mdl,n) # +1 is already appended onto ts

  n.DMatrix = [Matrix((n.Nck[int]),(n.Nck[int]+1)) for int in 1:n.Ni];
#  DM = [Matrix((Nck[int]),(Nck[int])) for int in 1:Ni];
  for int in 1:n.Ni
      n.DMatrix[int] = D[int][1:end-1,:];   # [Nck]X[Nck+1]
      #TODO make IMatrix and option
      #DM[int] = D[int][1:end-1,2:end];         # [Nck]X[Nck]
    #  IMatrix[int] = inv(DM[int]);        # I = inv[D_{2:N_k+1}]
  end
  return n
end
