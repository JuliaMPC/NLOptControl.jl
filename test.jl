using NLOptControl
using Polynomials
using Plots
using FastGaussQuadrature
using Parameters
pyplot()

"""
ps, nlp = initialize_NLP(Nc,Ni,numStates,numControls);
ps, nlp = initialize_NLP(Nc,Ni,numStates,numControls,stateVector,controlVector,decisionVector, t0, tf);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 1/1/2017 \n
Citations: \n
----------\n
Original Author: S. Hughes.  steven.p.hughes@nasa.gov
Source: DecisionVector.m [located here](https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/@DecisionVector/)
--------------------------------------------------------------------------\n
"""

# should only call this once
function initialize_NLP(Nc::Int64,Ni::Int64,numStates::Int64,numControls::Int64, args...)

    # validate input
    if  Nc <= 0
        error("Nc must be > 0");
    end
    if  Ni <= 0
        error("Ni must be > 0");
    end
    if  numStates <= 0
        error("numStates must be > 0");
    end
    if  numControls <= 0
        error("eventually numControls must be > 0");
    end

    # initialize node data TODO -> eventually make different PS methods available
    τ, ω = gaussradau(Nc);

    # calculate general properties
    numStatePoints = Nc*Ni;
    numControlPoints = Nc*Ni;

    # calculate length of vectors
    lengthStateVector = numStatePoints*numStates;
    lengthControlVector = numControlPoints*numControls;
    lengthDecVector = lengthStateVector + lengthControlVector + 2; # + 2 is for t0 and tf

    if length(args) == 0
        stateVector = zeros(lengthStateVector,);
        controlVector = zeros(lengthControlVector,);
        decisionVector = zeros(lengthDecVector,);
        t0 = 0.0;
        tf = 0.0;
    else
        # validate optional input
        if  length(args[1]) != lengthStateVector
            error(string("length of stateVector must be = ",lengthStateVector));
        end
        if  length(args[2]) != lengthControlVector
            error(string("length of controlVector must be = ",lengthControlVector));
        end
        if  length(args[3]) != lengthDecVector
            error(string("length of decisionVector must be = ",lengthDecVector));
        end
        if  length(args[4]) != 1
            error(string("length of t0 must be = 1"));
        end
        if  length(args[5]) != 1
            error(string("length of tf must be = 1"));
        end
        stateVector = args[1];
        controlVector = args[2];
        decisionVector = args[3];
        t0 = args[4];
        tf = args[5];
    end

    # determine indecies within overall decision vector of all variables (i.e. decisionVector)
    stateStartIdx = 1;
    stateStopIdx = stateStartIdx + lengthStateVector -1; # -1 because we start on 1
    controlStartIdx = stateStopIdx + 1;
    controlStopIdx = controlStartIdx + lengthControlVector -1; # -1 because we start on 1
    timeStartIdx = controlStopIdx + 1;
    timeStopIdx = timeStartIdx + 1;

    if 1==1 #TODO eventually make print_level an option
      print(string("lengthStateVector = ", lengthStateVector),"\n")
      print(string("lengthControlVector = ", lengthControlVector),"\n")

      print(string("stateStartIdx = ", stateStartIdx),"\n")
      print(string("stateStopIdx = ", stateStopIdx),"\n")
      print(string("controlStartIdx = ", controlStartIdx),"\n")
      print(string("controlStopIdx = ", controlStopIdx),"\n")
      print(string("timeStartIdx = ", timeStartIdx),"\n")
      print(string("timeStopIdx = ", timeStopIdx),"\n")
    end

    # check indecies
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

    # initialize problem data
    ps = PS_data(Nc=Nc,
                      Ni=Ni,
                       τ=τ,
                       ω=ω,
                      t0=t0,
                      tf=tf);

    nlp = NLP_data(numStates=numStates,
                        numStatePoints=numStatePoints,
                        numControls=numControls,
                        numControlPoints=numControlPoints,
                        lengthControlVector=lengthControlVector,
                        lengthStateVector=lengthStateVector,
                        lengthDecVector=lengthDecVector,
                        stateStartIdx=stateStartIdx,
                        stateStopIdx=stateStopIdx,
                        controlStartIdx=controlStartIdx,
                        controlStopIdx=controlStopIdx,
                        timeStartIdx=timeStartIdx,
                        timeStopIdx=timeStopIdx,
                        stateVector=stateVector,
                        controlVector=controlVector,
                        decisionVector=decisionVector
                        );

    return ps, nlp
end

Nc = 3;
Ni = 1;
numStates = 1;
numControls = 1;

ps, nlp = initialize_NLP(Nc,Ni,numStates,numControls);
