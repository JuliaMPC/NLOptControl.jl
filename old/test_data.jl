############
# TEST DATA
############
function generate_Fake_data(nlp::NLP_data,ps::PS_data,γS,γC)
    @unpack decisionVector = nlp
    @unpack Ni, t0, tf, ts = ps

    # state data
    @unpack lengthStateVector, numStates, stateIdx, stateIdx_all = nlp
    if numStates > 1     # each row contains all of the data for an interval
      fake_state = [[polyval(γS[st],ts[int]),numStates-1] for int in 1:Ni for st in 1:numStates];
    else
      fake_state = [polyval(γS[1],ts[int]) for int in 1:Ni ];
    end

    fake_state_data = zeros(lengthStateVector,);
    for idx in 1:Ni*numStates
      if numStates > 1
        fake_state_data[stateIdx_all[idx][1]:stateIdx_all[idx][2]] = fake_state[idx][1];
      else
        fake_state_data[stateIdx[idx][1]:stateIdx[idx][2]] = fake_state[idx]; # was int...a bug
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

    # control data
    @unpack lengthControlVector, numControls, controlIdx, controlIdx_all = nlp
    if numControls > 1     # each row contains all of the data for an interval
      fake_control = [[polyval(γC[ctr],ts[int][1:end-1]),numControls-1] for int in 1:Ni for ctr in 1:numControls];
    else
      fake_control = [polyval(γC[1],ts[int][1:end-1]) for int in 1:Ni ];
    end

    fake_control_data = zeros(lengthControlVector,);
    for idx in 1:Ni*numControls  #TODO these indices do not start at 1!! carful here
      if numControls > 1
        offset = controlIdx_all[1][1]-1; # offset tupples to start them at 1
        fake_control_data[controlIdx_all[idx][1]-offset:controlIdx_all[idx][2]-offset] = fake_control[idx][1];
      else
        offset = controlIdx[1][1]-1; # offset tupples to start them at 1
        fake_control_data[controlIdx[idx][1]-offset:controlIdx[idx][2]-offset] = fake_control[idx];
      end
    end

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

    decisionVector=[fake_state_data[:];fake_control_data[:];t0;tf]; # for now looking at no controls
    @pack  nlp = decisionVector
end
############
# TEST DATA
############
