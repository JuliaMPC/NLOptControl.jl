"""
ps, nlp = initialize_NLP(numStates=2,numControls=2,Ni=4,Nck=[3, 3, 7, 2]);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 1/9/2017 \n
Citations: \n
----------\n
Influenced by: S. Hughes.  steven.p.hughes@nasa.gov
Source: DecisionVector.m [located here](https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/@DecisionVector/)
-------------------------------------------------------------------------------------\n
"""
# should only call this once
function initialize_NLP(
  args...;
  numStates::Int64=0,
  numControls::Int64=0,
  Ni::Int64=10,
  Nck::Array{Int64,1}=4*ones(Int64,Ni,),
  stateEquations::Function=[],
  X0::Array{Float64,1}=zeros(Float64,numStates,1),
  XF::Array{Float64,1}=zeros(Float64,numStates,1),
  XL::Array{Float64,1}=zeros(Float64,numStates,1),
  XU::Array{Float64,1}=zeros(Float64,numStates,1),
  CL::Array{Float64,1}=zeros(Float64,numControls,1),
  CU::Array{Float64,1}=zeros(Float64,numControls,1),
  kwargs...)

    # validate input
    if length(Nck) != Ni
        error("\n length(Nck) != Ni \n");
    end
    for int in 1:Ni
        if (Nck[int]<0)
            error("\n Nck must be > 0");
        end
    end
    if  Ni <= 0
        error("\n Ni must be > 0 \n");
    end
    if  numStates <= 0
        error("\n numStates must be > 0","\n",
                "default value = 0","\n",
              );
    end
    if  numControls <= 0
        error("eventually numControls must be > 0","\n",
              "default value = 0","\n",
              );
    end
    if length(X0) != numStates
      error(string("\n Length of X0 must match number of states \n"));
    end
    if length(XF) != numStates
      error(string("\n Length of XF must match number of states \n"));
    end
    if length(XL) != numStates
      error(string("\n Length of XL must match number of states \n"));
    end
    if length(XU) != numStates
      error(string("\n Length of XU must match number of states \n"));
    end
    if length(CL) != numControls
      error(string("\n Length of CL must match number of controls \n"));
    end
    if length(CU) != numControls
      error(string("\n Length of CU must match number of controls \n"));
    end
    # initialize node data TODO -> eventually make different PS methods available

    taus_and_weights = [gaussradau(Nck[int]) for int in 1:Ni];
    τ = [taus_and_weights[int][1] for int in 1:Ni];
    ω = [taus_and_weights[int][2] for int in 1:Ni];

    # initialize scaled variables as zeros
    ts = 0*τ;
    ωₛ = 0*ω;

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
            error(string("Length of decisionVector must be = ",lengthDecVector),"\n");
        end
        if  length(args[2]) != 1
            error(string("\n Length of t0 must be = 1 \n"));
        end
        if  length(args[3]) != 1
            error(string("\n Length of tf must be = 1 \n"));
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

    if numStates == 1
      stateIdx_all = [(-99,-99) for int in 1:1];
      stateIdx_st = [(-99,-99) for int in 1:1];
    else  # only calculate these if we have more than one state
       # Break Tuples For Each State Variable within Entire Mesh Grid
      stateIdx_all = [(stateIdx[int][1] - numStatePoints[int]*(st-numStates),
                       stateIdx[int][2] - numStatePoints[int]*(st-numStates +1))
                                      for int in 1:Ni for st in numStates:-1:1]  # does the outer loop first

      if Ni == 1
        stateIdx_st = [(-99,-99) for int in 1:1];
      else         # Organize Tuples by Individual States
        organize_state_array = zeros(Int64,Ni*numStates,); idx=1;
      for int in 1:Ni
        for st in 1:numStates
          organize_state_array[idx] = int + (st-1)*numStates;
          idx=idx+1;
        end
      end
      stateIdx_st = [( stateIdx_all[jj][1], # all states are near each other
                     stateIdx_all[jj][2])   # TODO something is wrong with this!!
            for jj in organize_state_array]
      end
    end
 #TODO try to make a single way to reference states
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

    if numControls == 1
      controlIdx_all = [(-99,-99) for int in 1:1];
      controlIdx_ctr = [(-99,-99) for int in 1:1];
    else # only calculate these if we have more than one control
      # Break Tuples For Each Control Variable within Entire Mesh Grid
      controlIdx_all = [(controlIdx[int][1] - numControlPoints[int]*(ctr-numControls),
                       controlIdx[int][2] - numControlPoints[int]*(ctr-numControls +1))
      for int in 1:Ni for ctr in numControls:-1:1]  # does the outer loop first
        if Ni == 1
          controlIdx_ctr = [(-99,-99) for int in 1:1];
        else           # Organize Tuples by Individual Controls
          organize_control_array = zeros(Int64,Ni*numControls,); idx=1;
          for int in 1:Ni
            for ctr in 1:numControls
              organize_control_array[idx] = int + (ctr-1)*numControls;
              idx=idx+1;
            end
          end
          controlIdx_ctr = [(controlIdx_all[jj][1], # all controls are near each other
                            controlIdx_all[jj][2])
                    for jj in organize_control_array]
        end
    end
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
    # each row contains Xi in the stateMatrix where the size = Nck[int]XnumStates
    # Xi is a vector of ALL of the states at point i
    stateMatrix=[zeros((Nck[int]+1),numStates) for int in 1:Ni];
    controlMatrix=[zeros(Nck[int],numControls) for int in 1:Ni];

    DMatrix = [zeros((Nck[int]),(Nck[int]+1)) for int in 1:Ni];
    IMatrix = [zeros((Nck[int]),(Nck[int])) for int in 1:Ni];
    FMatrix = [zeros((Nck[int]),numStates) for int in 1:Ni];

    odeConstraint = [zeros(Float64,Nck[int],numStates) for int in 1:Ni]; # may also want to call this stateEqConstraint
    continuityConstraint = [zeros(Float64,numStates+numControls,) for int in 1:Ni-1];
    boundaryConstraint = zeros(Float64,2*numStates,);
    stateConstraint = [zeros(Float64,sum(numStatePoints),2) for st in 1:numStates];
    controlConstraint = [zeros(Float64,sum(numControlPoints),2) for ctr in 1:numControls];
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
  ps = PS_data(           Nck=Nck,
                           Ni=Ni,
                            τ=τ,
                           ts=ts,
                            ω=ω,
                           ωₛ=ωₛ,
                           t0=t0,
                           tf=tf,
                  stateMatrix=stateMatrix,
                controlMatrix=controlMatrix,
                      DMatrix=DMatrix,
                      IMatrix=IMatrix,
                      FMatrix=FMatrix,
                 odeConstraint=odeConstraint,  # might also be state constraint, but this is the same as the variable name below
         continuityConstraint=continuityConstraint,
           boundaryConstraint=boundaryConstraint,
              stateConstraint=stateConstraint,
            controlConstraint=controlConstraint
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
                      decisionVector=decisionVector,
                      stateEquations=stateEquations,
                      X0=X0,
                      XF=XF,
                      XL=XL,
                      XU=XU,
                      CL=CL,
                      CU=CU
                 );
    return ps, nlp
end
"""
nlp2ocp(nlp,ps);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/2/2017, Last Modified: 1/4/2017 \n
--------------------------------------------------------------------------\n
"""
# to do, pack decisionVector into nlp
function nlp2ocp(nlp::NLP_data,ps::PS_data)
    @unpack t0, tf, stateMatrix, controlMatrix, Ni, Nck = ps
    @unpack stateIdx_all, controlIdx_all, timeStartIdx, timeStopIdx = nlp
    @unpack stateIdx, controlIdx = nlp
    @unpack numStates, numControls, lengthDecVector = nlp
    @unpack decisionVector = nlp

    if length(decisionVector)!=lengthDecVector
      error(string("\n",
                    "-------------------------------------", "\n",
                    "There is an error with the indecies!!", "\n",
                    "-------------------------------------", "\n",
                    "The following variables should be equal:", "\n",
                    "length(decisionVector) = ",length(decisionVector),"\n",
                    "lengthDecVector = ",lengthDecVector,"\n"
                    )
            )
    end
    # update parameters
    t0 = decisionVector[timeStartIdx];
    tf = decisionVector[timeStopIdx];

    # the state matrix is sized according to eq. (40) in the GPOPS II article
    # n is the total number of states -> the individual states are columns
    # V[int]      = [X11               X21      ...      Xn1;
    #                X12               X22      ...      Xn2;
    #                .                  .                 .
    #                .                  .                 .
    #                .                  .                 .
    #         X1_{Nck[int]+1}    X2_{Nck[int]+1}   Xn_{Nck[int]+1}

    stateMatrix = [zeros(Nck[int]+1, numStates) for int in 1:Ni];
    idx = 1;
    for int in 1:Ni
        for st in 1:numStates
            if numStates > 1
              stateMatrix[int][:,st] = decisionVector[stateIdx_all[idx][1]:stateIdx_all[idx][2]]
            else # use indexing for single state variable
              stateMatrix[int][:,st] = decisionVector[stateIdx[idx][1]:stateIdx[idx][2]]
           end
          idx+=1;
        end
    end

    controlMatrix = [zeros(Nck[int], numControls) for int in 1:Ni];
    idx = 1;
    for int in 1:Ni
        for ctr in 1:numControls
            if numControls > 1
              controlMatrix[int][:,ctr] = decisionVector[controlIdx_all[idx][1]:controlIdx_all[idx][2]];
            else
              controlMatrix[int][:,ctr] = decisionVector[controlIdx[idx][1]:controlIdx[idx][2]];
            end
            idx+=1;
        end
    end
    @pack ps = t0, tf, stateMatrix, controlMatrix
end
