"""
LGR_matrices(ps,nlp)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/3/2017, Last Modified: 1/4/2017 \n
--------------------------------------------------------------------------------------\n
"""
function LGR_matrices(ps::PS_data,nlp::NLP_data)
    @unpack Nck, Ni, ts = ps
    @unpack DMatrix, IMatrix = ps

    # check input
    check_ts = zeros(Float64,Ni);
    for int in 1:Ni
        check_ts[int]=maximum(ts[int]);
    end
    if maximum(check_ts) < 10*eps()
      error("\n ts is full of zeros: make sure that you call create_intervals() first to calculate ts! \n
                NOTE: This may have occured because  (:finalTimeDV => true) and the final time dv is not working properly!! \n")
    end

    D = [zeros((Nck[int]+1),(Nck[int]+1)) for int in 1:Ni];
    for int in 1:Ni
        D[int] = poldif(ts[int], 1) # +1 is already appended onto ts
    end

    DM = [zeros((Nck[int]),(Nck[int])) for int in 1:Ni];
    for int in 1:Ni
        DMatrix[int] = D[int][1:end-1,:];   # [Nck]X[Nck+1]
        DM[int] = D[int][1:end-1,2:end];    # [Nck]X[Nck]
        IMatrix[int] = inv(DM[int]);        # I = inv[D_{2:N_k+1}]
    end
    @pack ps = DMatrix, IMatrix
end
"""
DMatrix_JuMP=D_matrix_JuMP(tf_var,Nck_const,Ni_const,τ_const,ω_const)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/15/2017, Last Modified: 1/16/2017 \n
--------------------------------------------------------------------------------------\n
"""
function D_matrix_JuMP(tf_var,Nck_const,Ni_const,τ_const,ω_const)
  # calculate ts and ωₛ based off of tf
  ts_JuMP, ωₛ_JuMP = create_intervals_JuMP(tf_var,Nck_const,Ni_const,τ_const,ω_const);

  Nck=Nck_const; Ni=Ni_const;
  D = [Matrix((Nck[int]+1),(Nck[int]+1)) for int in 1:Ni];
  for int in 1:Ni
      D[int] = poldif_JuMP(ts_JuMP[int]) # +1 is already appended onto ts
  end

  DMatrix_JuMP = [Matrix((Nck[int]),(Nck[int]+1)) for int in 1:Ni];
  DM = [Matrix((Nck[int]),(Nck[int])) for int in 1:Ni];
  for int in 1:Ni
      DMatrix_JuMP[int] = D[int][1:end-1,:];   # [Nck]X[Nck+1]
      DM[int] = D[int][1:end-1,2:end];         # [Nck]X[Nck]
      #TODO make IMatrix and option
    #  IMatrix[int] = inv(DM[int]);        # I = inv[D_{2:N_k+1}]
  end
  return DMatrix_JuMP
end
