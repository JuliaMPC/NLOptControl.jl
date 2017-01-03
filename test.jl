using NLOptControl
using Plots
using Polynomials
using FastGaussQuadrature
using Parameters

# initialize basic problem definition
ps, nlp = initialize_NLP(numStates=1,numControls=2,Ni=2,Nck=[3, 5]);

####################################
# perform analytical calcualtions -> for plots
####################################
@unpack t0, tf = ps
t0 = Float64(0); tf = Float64(10);
@pack ps = t0, tf
t = Array(linspace(t0,tf,100));
α₁ =  -0.3; α₂ = 3; α₃ = -8; α₄ =  7;

γ = Poly([α₁,α₂,α₁]); #TODO check on that imported binding warning
y = polyval(γ,t);

# evaluate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,t[end]) - polyval(∫γ,t[1]);
C = Y - polyval(∫γ,t[end]); # constant of integration
∫y = polyval(∫γ,t) + C;

# evaluate the derivative
dγ = polyder(γ);
dy = polyval(dγ,t);

##TEMP## -->to get fake optimization data
@unpack Nck, Ni, t0, tf = ps
taus_and_weights = [gaussradau(Nck[int]) for int in 1:Ni];
τ = [taus_and_weights[int][1] for int in 1:Ni];
ω = [taus_and_weights[int][2] for int in 1:Ni];
di, tm, t_data, ωₛ=create_intervals(t0,tf,Ni,Nck,τ,ω);
@pack ps = τ, ω


@unpack decisionVector, lengthControlVector = nlp
fake_control_data = zeros(lengthControlVector,);
ts = [t_data[1]; t_data[2]]
if Ni > 2
    error("fix ts")
end
fake_state_data = polyval(γ,ts);
decisionVector=[fake_state_data[:];fake_control_data;t0;tf]; # for now looking at no controls
@pack  nlp = decisionVector

##TEMP## -->to get fake optimization data

function nlp2ocp(decisionVector::Array{Float64,1},nlp::NLP_data,ps::PS_data)
    @unpack t0, tf, stateMatrix, controlMatrix, Ni = ps
    @unpack stateIdx = nlp
    @unpack controlIdx = nlp
    @unpack timeStartIdx,timeStopIdx = nlp

    # update parameters
    t0 = decisionVector[timeStartIdx];
    tf = decisionVector[timeStopIdx];

    stateMatrix  = [decisionVector[stateIdx[int][1]:stateIdx[int][2]] for int in 1:Ni];
    controlMatrix = [decisionVector[controlIdx[int][1]:controlIdx[int][2]] for int in 1:Ni];

    @pack ps = t0, tf, stateMatrix, controlMatrix
end

nlp2ocp(decisionVector,nlp,ps);

@unpack_PS_data ps
@unpack_NLP_data nlp
print(nlp,"\n","\n")
print(ps)

#=
using NLOptControl
ps, nlp = initialize_NLP(numStates=1,numControls=3,Ni=2,Nck=[2, 4]);
@unpack_NLP_data nlp
print(nlp)




=#
