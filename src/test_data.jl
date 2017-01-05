############
# TEST DATA
############
function generate_Fake_data(nlp::NLP_data,ps::PS_data,γ)
    @unpack decisionVector, lengthControlVector, lengthStateVector, numStates, stateIdx, stateIdx_all = nlp
    fake_control_data = zeros(lengthControlVector,);

    @unpack Ni, t0, tf, ts = ps

    if length(fake_control_data)!=lengthControlVector
      error(string("\n",
                    "-------------------------------------", "\n",
                    "There is an error with the indecies!!", "\n",
                    "-------------------------------------", "\n",
                    "The following variables should be equal:", "\n",
                    "length(fake_control_data) = ",length(fake_control_data),"\n",
                    "lengthControlVector = ",lengthControlVector,"\n"
                    )
            )
    end

    # each row contains all of the data for an interval
    if numStates > 1
      fake_state = [[polyval(γ,ts[int]),zeros(Float64,length(ts[int]),numStates-1)] for int in 1:Ni ];
    else
      fake_state = [polyval(γ,ts[int]) for int in 1:Ni ];
    end

    fake_state_data = zeros(lengthStateVector,);
    idx=1;
    for int in 1:Ni
        for st in 1:numStates # turn into vector
          if numStates > 1
            fake_state_data[stateIdx_all[idx][1]:stateIdx_all[idx][2]] = fake_state[int][st];
          else
            fake_state_data[stateIdx[idx][1]:stateIdx[idx][2]] = fake_state[int];
          end
         idx+=1;
      end
    end
    if length(fake_state_data)!=lengthStateVector
      error(string("\n",
                    "-------------------------------------", "\n",
                    "There is an error with the indecies!!", "\n",
                    "-------------------------------------", "\n",
                    "The following variables should be equal:", "\n",
                    "length(fake_state_data) = ",length(fake_state_data),"\n",
                    "lengthStateVector = ",lengthStateVector,"\n"
                    )
            )
    end
    decisionVector=[fake_state_data[:];fake_control_data;t0;tf]; # for now looking at no controls
    @pack  nlp = decisionVector
end
############
# TEST DATA
############
