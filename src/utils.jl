"""
ps, nlp = initialize_NLP(numStates=1,numControls=3);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 1/2/2017 \n
Citations: \n
----------\n
Influenced by: S. Hughes.  steven.p.hughes@nasa.gov
Source: DecisionVector.m [located here](https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/@DecisionVector/)
--------------------------------------------------------------------------\n
"""

# should only call this once
function initialize_NLP(args...;numStates::Int64=0,numControls::Int64=0,Ni::Int64=10,Nck::Array{Int64,1}=4*ones(Int64,Ni,), kwargs...)
    # validate input
    if length(Nck) != Ni
        error("length(Nck) != Ni");
    end
    for int in 1:Ni
        if (Nck[int]<0)
            error("Nck must be > 0");
        end
    end
    if  Ni <= 0
        error("Ni must be > 0");
    end
    if  numStates <= 0
        error("numStates must be > 0","\n",
                "default value = 0","\n",
              );
    end
    if  numControls <= 0
        error("eventually numControls must be > 0","\n",
              "default value = 0","\n",
              );
    end

    # initialize node data TODO -> eventually make different PS methods available

    taus_and_weights = [gaussradau(Nck[int]) for int in 1:Ni];
    τ = [taus_and_weights[int][1] for int in 1:Ni];
    ω = [taus_and_weights[int][2] for int in 1:Ni];

    # calculate general properties
    numStatePoints = [(Nck[int]+1) for int in 1:Ni]; # number of states dvs (decisionVector) for each interval and each state
    numControlPoints = [Nck[int] for int in 1:Ni];   # number of control dvs for each interval and each control

    # calculate length of vectors
    lengthStateVector = sum([numStatePoints[int]*numStates for int in 1:Ni]);
    lengthControlVector = sum([numControlPoints[int]*numControls for int in 1:Ni]);
    lengthDecVector = lengthStateVector + lengthControlVector + 2; # + 2 is for t0 and tf

    if length(args) == 0
        decisionVector = zeros(lengthDecVector,);
        t0 = 0.0;
        tf = 0.0;
    else
        # validate optional input
        if  length(args[1]) != lengthDecVector
            error(string("length of decisionVector must be = ",lengthDecVector));
        end
        if  length(args[2]) != 1
            error(string("length of t0 must be = 1"));
        end
        if  length(args[3]) != 1
            error(string("length of tf must be = 1"));
        end
        decisionVector = args[1];
        t0 = args[2];
        tf = args[3];
    end
    ####################################################################################
    ############################## Indices #############################################
    ####################################################################################

    # ==================================================================================
    #__________________________ State Indices ___________________________________________
    # ===================================================================================
    # General Properties
    stateStartIdx = 1;
    stateStopIdx = stateStartIdx + lengthStateVector -1; # -1 because we start on 1
    st_sum = [0; cumsum(numStatePoints)];

    # Organize Tuples Each Entire Mesh Grid Of All Consecutive State Vectors
    stateIdx = [((st_sum[int])*numStates + stateStartIdx,
    (st_sum[int] + numStatePoints[int])*numStates + stateStartIdx -1)
    for int in 1:Ni];

    # Break Tuples For Each State Variable within Entire Mesh Grid
    stateIdx_all = [(stateIdx[int][1] - numStatePoints[int]*(st-numStates),
                     stateIdx[int][2] - numStatePoints[int]*(st-numStates +1))
    for int in 1:Ni for st in numStates:-1:1]  # does the outer loop first

    # Organize Tuples by Individual States
    organize_state_array = zeros(Int64,Ni*numStates,); idx=1;
    for int in 1:Ni
      for st in 1:numStates
        organize_state_array[idx] = int + (st-1)*numStates;
        idx=idx+1;
      end
    end
    stateIdx_st = [( stateIdx_all[jj][1], # all states are near each other
                     stateIdx_all[jj][2])
    for jj in organize_state_array]
    # ==================================================================================
    #_________________________ Control Indices _________________________________________
    # ==================================================================================
    # General Properties
    controlStartIdx = stateStopIdx + 1;
    controlStopIdx = controlStartIdx + lengthControlVector -1; # -1 because we start on 1
    ctr_sum = [0; cumsum(numControlPoints)];

    # Organize Tuples Each Entire Mesh Grid Of All Consecutive Control Vectors
    controlIdx = [((ctr_sum[int])*numControls + controlStartIdx,
     (ctr_sum[int] + numControlPoints[int])*numControls + controlStartIdx -1)
    for int in 1:Ni]

    # Break Tuples For Each Control Variable within Entire Mesh Grid
    controlIdx_all = [(controlIdx[int][1] - numControlPoints[int]*(ctr-numControls),
                     controlIdx[int][2] - numControlPoints[int]*(ctr-numControls +1))
    for int in 1:Ni for ctr in numControls:-1:1]  # does the outer loop first

    # Organize Tuples by Individual Controls
    organize_control_array = zeros(Int64,Ni*numControls,); idx=1;
    for int in 1:Ni
      for ctr in 1:numControls
        organize_control_array[idx] = int + (ctr-1)*numControls;
        idx=idx+1;
      end
    end
    controlIdx_ctr = [(controlIdx_all[jj][1], # all controls are near each other
                      controlIdx_all[jj][2])
    for jj in organize_state_array]
    # ==================================================================================
    #___________________________ Time Indies ____________________________________________
    # ===================================================================================
    timeStartIdx = controlStopIdx + 1;
    timeStopIdx = timeStartIdx + 1;
    # ==================================================================================
    #___________________________ Check Indices __________________________________________
    # ===================================================================================
    if timeStopIdx != lengthDecVector
      error(string("\n",
                    "-------------------------------------", "\n",
                    "There is an error with the indecies!!", "\n",
                    "-------------------------------------", "\n",
                    "The following variables should be equal:", "\n",
                    "timeStopIdx = ",timeStopIdx,"\n",
                    "lengthDecVector = ",lengthDecVector,"\n"
                    )
            )
    end
    ####################################################################################
    ############################## Matrices ############################################
    ####################################################################################
    # CAN ONLY DO THIS IF to and tf where provided
    # each row contains Xi in the stateMatrix where the size = Nck[int]XnumStates
    # Xi is a vector of ALL of the states at point i
    stateMatrix=[zeros((Nck[int]+1)*numStates,) for int in 1:Ni];
    controlMatrix=[zeros(Nck[int]*numControls,) for int in 1:Ni];

    # approximate the derivative --> needed in defect constraints
    D = [zeros((Nck[int]+1),(Nck[int]+1)) for int in 1:Ni]
    for int in 1:Ni
      D[int] = poldif([τ[int];1], 1) # append +1 onto τ
    end
    DMatrix = [zeros((Nck[int]),(Nck[int]+1)) for int in 1:Ni];
    DM = [zeros((Nck[int]),(Nck[int])) for int in 1:Ni];
    IMatrix = DM;
    for int in 1:Ni
      DMatrix[int] = D[int][1:end-1,:];   # turn into a [Nck]X[Nck+1] sized matrix
      DM[int] = D[int][1:end-1,:1:end-1]; # turn into a [Nck]X[Nck] sized matrix
      IMatrix[int] = inv(DM[int]);        # integration matrix
    end
    # ==================================================================================
    #___________________________ Debugging _____________________________________________
    # ===================================================================================
    if false #TODO  make print_level an option
      print(string("lengthStateVector = ", lengthStateVector),"\n")
      print(string("lengthControlVector = ", lengthControlVector),"\n")
      print(string("stateStartIdx = ", stateStartIdx),"\n")
      print(string("stateStopIdx = ", stateStopIdx),"\n")
      print(string("controlStartIdx = ", controlStartIdx),"\n")
      print(string("controlStopIdx = ", controlStopIdx),"\n")
      print(string("timeStartIdx = ", timeStartIdx),"\n")
      print(string("timeStopIdx = ", timeStopIdx),"\n")

      print(string("typeof(Nck) = ",typeof(Nck),"\n"))
      print(string("typeof(Ni) = ",typeof(Ni),"\n"))
      print(string("typeof(τ) = ",typeof(τ),"\n"))
      print(string("typeof(ω) = ",typeof(ω),"\n"))
      print(string("typeof(t0) = ",typeof(t0),"\n"))
      print(string("typeof(tf) = ",typeof(tf),"\n"))
      print(string("typeof(stateMatrix) = ",typeof(stateMatrix),"\n"))
      print(string("typeof(controlMatrix) = ",typeof(controlMatrix),"\n"))
      print(string("typeof(DMatrix) = ",typeof(DMatrix),"\n"))
      print(string("typeof(IMatrix) = ",typeof(IMatrix),"\n"))
  end
  # ==================================================================================
  #_________________________ Initialize Problem Data __________________________________
  # ===================================================================================
  ps = PS_data(Nck=Nck,
               Ni=Ni,
                τ=τ,
                ω=ω,
               t0=t0,
               tf=tf,
      stateMatrix=stateMatrix,
    controlMatrix=controlMatrix,
          DMatrix=DMatrix,
          IMatrix=IMatrix
              );
  nlp = NLP_data(     numStates=numStates,
                      numControls=numControls,
                      numStatePoints=numStatePoints,
                      numControlPoints=numControlPoints,
                      lengthControlVector=lengthControlVector,
                      lengthStateVector=lengthStateVector,
                      lengthDecVector=lengthDecVector,
                      stateIdx=stateIdx,
                      controlIdx=controlIdx,
                      stateIdx_all=stateIdx_all,
                      controlIdx_all=controlIdx_all,
                      controlIdx_ctr=controlIdx_ctr,
                      stateIdx_st=stateIdx_st,
                      timeStartIdx=timeStartIdx,
                      timeStopIdx=timeStopIdx,
                      decisionVector=decisionVector
                 );
    return ps, nlp
end

function scale_tau(τ::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*τ + (xₙ + x₀)/2;
end

function scale_w(ω::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*ω;
end

"""
di, tm, ts, ωₛ=create_intervals(t0,tf,Ni,Nck,τ,ω);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 1/2/2017 \n
--------------------------------------------------------------------------\n
"""
function create_intervals(t0::Float64,tf::Float64,Ni::Int64,Nck::Array{Int64,1},τ::Array{Vector{Float64},1},ω::Array{Vector{Float64},1})
  di = (tf - t0)/Ni; # interval size
  # create mesh points
  tm = zeros(Float64,Ni+1); tm[1] = t0;
  for idx in 1:Ni
    tm[idx+1] = tm[idx] + di;
  end
    # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
    ts=[[scale_tau(τ[int],tm[int],tm[int+1]);di*int] for int in 1:Ni];
    ωₛ=[scale_w(ω[int],tm[int],tm[int+1]) for int in 1:Ni];
    return di, tm, ts, ωₛ
end

"""
L = lagrange_basis_poly(x,x_data,Nc,j)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 12/28/2016
Citations: This function was influenced by the lagrange() function [located here](https://github.com/pjabardo/Jacobi.jl/blob/master/src/gauss_quad.jl)
--------------------------------------------------------------------------------------\n
# Input Arguments
* `x`: point to approximate function value at
* `x_data`: x data to used calculate basis polynomials
* `Nc`: order of Lagrange interpolating polynomial
* `j`: index of interest

# Output Arguments
* `L`: Lagrange basis polynomials

A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

"""
function lagrange_basis_poly(x,x_data,Nc,j)
    L = 1;
    for idx in 1:Nc+1 # use all of the data
      if idx!=j
        L = L*(x - x_data[idx])/(x_data[j]-x_data[idx]);
      end
    end
  return L
end

function lagrange_basis_poly!{T<:Number}(x::AbstractArray{T},x_data,Nc,L::AbstractArray{T})
   if Nc > length(x_data) -1
      error("Maximum N value = length(x_data)-1")
    end
    ns = length(x);
    L = zeros(Float64,Nc+1,ns);
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
      end
    end
    return L
end
lagrange_basis_poly{T<:Number}(x::AbstractArray{T},x_data,Nc) = lagrange_basis_poly!(x::AbstractArray{T},x_data,Nc,zeros(x))

function interpolate_lagrange{T<:Number}(x::AbstractArray{T},x_data,y_data,Nc)
    if Nc > length(x_data) -1
      error("Maximum N value = length(x_data)-1")
    end
    ns = length(x);
    L = zeros(Float64,Nc+1,ns);
    x = x[:]; x_data = x_data[:]; y_data = y_data[:]; # make sure data is in a column
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
      end
    end
    y = y_data'*L;
    return y
end

"""
D = poldif(x, malpha, B...);
--------------------------------------------------------------------------\n
Last modifed for julia on December 23, 2016 by Huckleberry Febbo\n
Original Author: JJ.A.C. Weideman, S.C. Reddy 1998\n
Original Function Name: poldif.m  |  Source: [matlabcentral](https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite)\n
--------------------------------------------------------------------------\n
# Input Arguments
* `x::Vector`: Vector of N distinct nodes.
* `malpha::Int65`: M, the number of derivatives required
    * Note:     0 < M < N-1
*  `B`: Currently this functionality is not tested with this input (only works for constant weights)

# Output Arguments
* `D::Array{Float64,2}`: Differention Matrix
    * DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M

The function DM =  poldif(x, maplha, B) computes the differentiation matrices D1, D2, ..., DM on arbitrary nodes. The function is called with either two or three input arguments. If two input arguments are supplied, the weight function is assumed to be constant.   If three arguments are supplied, the weights should be defined as the second and third arguments. (CURRENTLY NOT TESTED IN julia)

"""

function poldif(x, malpha, B...)
  B_given = Bool(length(B));
         N = length(x);
         x = x[:];                     # Make sure x is a column vector

  if !B_given                       # Check if constant weight function
         M = malpha;                   # is to be assumed.
     alpha = ones(N,1);
         B = zeros(M,N);
  elseif B_given
     alpha = malpha(:);                # Make sure alpha is a column vector
         M = length(B[:,1]);           # First dimension of B is the number
  end                                  # of derivative matrices to be computed

        #  I = eye(N);                  # Identity matrix.
          #L = logical(I);              # Logical identity matrix.
          L=eye(Bool, N, N);
      # XX = x(:,ones(1,N));
         XX = repmat(x,1,N);
         DX = XX-XX';                  # DX contains entries x(k)-x(j).

      DX[L] = ones(N,1);               # Put 1's one the main diagonal.

          c = alpha.*prod(DX,2);       # Quantities c(j).

          #C = c[:,ones(1,N)];
          C = repmat(c,1,N);
          C = C./C';                   # Matrix with entries c(k)/c(j).

          Z = 1./DX;                   # Z contains entries 1/(x(k)-x(j))
       Z[L] = zeros(N,1);              # with zeros on the diagonal.

          X = Z';                      # X is same as Z', but with
       #X[L] = [];                      # diagonal entries removed.
       flag = trues(size(X));
    flag[L] = false;
          X = X[flag];
          X = reshape(X,N-1,N);

          Y = ones(N-1,N);             # Initialize Y and D matrices.
          D = eye(N);                  # Y is matrix of cumulative sums,

  DM = zeros(Float64,N,N,M);                                   # D differentiation matrices.
  for ell = 1:M
          temp = reshape(B[ell,:],1,N)
          Y   = cumsum([temp; ell*Y[1:N-1,:].*X]); # Diagonals
          D   = ell*Z.*(C.*repmat(diag(D),1,N) - D);   # Off-diagonals
       D[L]   = Y[N,:];                                # Correct the diagonal
  DM[:,:,ell] = D;                                     # Store the current D

  end

  return squeeze(DM,3)
end

# the following code uses the function lepoly() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
function lepoly(n::Int64,x,derivative_option::Bool)
# lepoly  Legendre polynomial of degree n
    # y=lepoly(n,x) is the Legendre polynomial
    # The degree should be a nonnegative integer
    # The argument x should be on the closed interval [-1,1];
    # [dy,y]=lepoly(n,x) also returns the values of 1st-order
    #  derivative of the Legendre polynomial stored in dy
# Last modified on August 30, 2011
# Verified with the chart in http://keisan.casio.com/has10/SpecExec.cgi
  if !derivative_option
       if n==0 y=ones(size(x));  return y end
       if n==1 y=x; return y end
       polylst=ones(size(x)); poly=x;   # L_0(x)=1, L_1(x)=x
       polyn=0;
       for k=2:n                      # Three-term recurrence relation:
  	     polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k; # kL_k(x)=(2k-1)xL_{k-1}(x)-(k-1)L_{k-2}(x)
         polylst=poly; poly=polyn;
       end
       y=polyn;
  else
     if n==0 y=ones(size(x)); dy=zeros(size(x)); return dy,y end
     if n==1 y=x; dy=ones(size(x)); return dy, y end

     # TODO MAKE SURE THIS IS CALCULATING SIZE correctly!! i.e. N=4 and size(N+1,N+1) = 1
      polylst=ones(size(x)); pderlst=zeros(size(x));poly=x; pder=ones(size(x));
       # L_0=1, L_0'=0, L_1=x, L_1'=1
      polyn=0;
      pdern=0;
      for k=2:n                          # Three-term recurrence relation:
        polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k; # kL_k(x)=(2k-1)xL_{k-1}(x)-(k-1)L_{k-2}(x)
        pdern=pderlst+(2*k-1)*poly; # L_k'(x)=L_{k-2}'(x)+(2k-1)L_{k-1}(x)
     	  polylst=poly; poly=polyn;
    	  pderlst=pder; pder=pdern;
      end
      y=polyn; dy=pdern;
      return dy, y
  end

end
