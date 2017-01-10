using NLOptControl
using Polynomials
using Plots
using Parameters
pyplot()

############
# TEST DATA
############
t0 = Float64(0); tf = Float64(2);
t = Array(linspace(t0,tf,100));

## STATES ##
# sin and cos taylor series approximation at x = 0
s0 =  0; s1 = 1; s2 = 0;       s3 = -1/(3*2*1); s4 = 0;           s5 = 1/(5*4*3*2*1);
c0 = 1;  c1 = 0; c2 = 1/(2*1); c3 = 0;          c4 = 1/(4*3*2*1); c5 = 0;
γ1S = Poly([s0,s1,s2,s3,s4,s5]);  # state 1
γ2S = Poly([c0,c1,c2,c3,c4,c5]);  # state 2
γS = [γ1S; γ2S]; numStates=length(γS);
yS = [polyval(γS[st],t) for st in 1:numStates];

# evaluate the integrals
∫γS = [polyint(γS[st]) for st in 1:numStates];
YS  = [polyval(∫γS[st],t[end]) - polyval(∫γS[st],t[1]) for st in 1:numStates];
CS  = [YS[st] - polyval(∫γS[st],t[end]) for st in 1:numStates]; # constant of integration
∫yS = [polyval(∫γS[st],t) + CS[st] for st in 1:numStates];

# evaluate the derivatives
dγS = [polyder(γS[st]) for st in 1:numStates];
dyS = [polyval(dγS[st],t) for st in 1:numStates];

## CONTROLS ##
# sin and cos taylor series approximation at x = 0
s0 =  0; s1 = 1; s2 = 0;       s3 = -1/(3*2*1); s4 = 0;           s5 = 1/(5*4*3*2*1);
c0 = 1;  c1 = 0; c2 = 1/(2*1); c3 = 0;          c4 = 1/(4*3*2*1); c5 = 0;
γ1C = Poly([s0,s1,s2,s3,s4,s5]);  # control 1
γ2C = Poly([c0,c1,c2,c3,c4,c5]);  # control 2
γC = [γ1C; γ2C]; numControls=length(γC);
yC = [polyval(γC[ctr],t) for ctr in 1:numControls];

# evaluate the integrals
∫γC = [polyint(γC[ctr]) for ctr in 1:numControls];
YC  = [polyval(∫γC[ctr],t[end]) - polyval(∫γC[ctr],t[1]) for ctr in 1:numControls];
CC  = [YC[ctr] - polyval(∫γC[ctr],t[end]) for ctr in 1:numControls]; # constant of integration
∫yC = [polyval(∫γC[ctr],t) + CC[ctr] for ctr in 1:numControls];

# evaluate the derivatives
dγC = [polyder(γC[ctr]) for ctr in 1:numControls];
dyC = [polyval(dγC[ctr],t) for ctr in 1:numControls];

############
# TEST DATA
############
dx1(t) = 1.0 - 0.5⋅t^2 + 0.0416667⋅t^4;
dx2(t) = 1.0⋅t + 0.166667⋅t^3;
dx = [dx1, dx2];
X0= [yS[1][1];yS[2][1]];
XF= [yS[1][end];yS[2][end]];
XL=[0,-Inf];     XU=[1/9,Inf]; # TODO allow for functions of these so we can calculate them on the fly!
#CL=[-Inf];  CU=[Inf];
CL=[-Inf,-Inf];  CU=[Inf,Inf];

ps, nlp = initialize_NLP(numStates=numStates,numControls=numControls,Ni=5,Nck=[4,4,4,4,4],stateEquations=dx,X0=X0,XF=XF,XL=XL,XU=XU,CL=CL,CU=CU);

@pack ps = t0, tf;  # given in problem def.
@unpack Nck, Ni, t0, tf, τ, ω = ps;
di, tm, ts, ωₛ = create_intervals(t0,tf,Ni,Nck,τ,ω);
@pack ps = τ, ω, ωₛ, ts;
generate_Fake_data(nlp,ps,γS,γC);
nlp2ocp(nlp,ps);

# interpolate function using Lagrange Polynomial | state
@unpack stateMatrix = ps;
@unpack numStates = nlp;
PS = [zeros(Float64,Nck[int]+1,numStates) for int in 1:Ni];
for st in 1:numStates
    for int in 1:Ni
        PS[int][:,st] = interpolate_lagrange(ts[int],ts[int],stateMatrix[int][:,st],Nck[int])
    end
end

# interpolate function using Lagrange Polynomial | control --> had to do Nck[int]-1
@unpack controlMatrix = ps;
@unpack numControls = nlp;
PC = [zeros(Float64,Nck[int]+1,numControls) for int in 1:Ni];
for ctr in 1:numControls
    for int in 1:Ni
        PC[int][:,ctr] = interpolate_lagrange(ts[int],ts[int][1:end-1],controlMatrix[int][:,ctr],Nck[int]-1)
    end
end
# approximate integral using quadrature
stInt1, stIntVal1, ctrInt1, ctrIntVal1 = integrate(ps,nlp)

# calculate LGR matrices - > IMatrix and DMatrix
LGR_matrices(ps,nlp)

# approximate integral using LGRIM
stInt2, stIntVal2, ctrInt2, ctrIntVal2 = integrate(ps,nlp;(:mode=>:LGRIM))

# approximate derivative using LGRDM
dζ = differentiate_state(ps,nlp);
#################

#################
# post processing state
#################
ls = 1.35;
lw = 4;

ip=plot(0,leg=:false)
for st in 1:numStates
  plot!(t,∫yS[st],label=string(@sprintf("act. = %0.3f",∫yS[st][end])," || st. #",st),w=lw);
end
for st in 1:numStates
    legend_bool=true;
    for int in 1:Ni
        if legend_bool
            scatter!(ts[int][1:end-1],stInt1[int][st,1:Nck[int]],marker = (:pentagon, 10, 0.9, :red),leg=:bottomright,label=string(@sprintf("QUAD. = %0.3f", stIntVal1[st])," || st. #",st))
            scatter!(ts[int][1:end-1],stInt2[int][st,1:Nck[int]], marker = (:star5, 10, 0.9, :green),leg=:bottomright,label=string(@sprintf("LGR  = %0.3f", stIntVal2[st])," || st. #",st))
        else # do not show legend a bunch of times
            scatter!(ts[int][1:end-1],stInt1[int][st,1:Nck[int]],marker = (:pentagon, 10, 0.9, :red),leg=:bottomright,label= "",leg=true)
            scatter!(ts[int][1:end-1],stInt2[int][st,1:Nck[int]], marker = (:star5, 10, 0.9, :green),leg=:bottomright,label= "",leg=true)
        end
        legend_bool=false;
    end
end
xlims!(t0,tf*ls)
ylabel!("Integral")
xlabel!(string("Number of Intervals = ", Ni))

dp=plot(0,leg=:false)

for st in 1:numStates
  plot!(t,dyS[st],label=string("act. state #",st),w=lw);
end
for st in 1:numStates
    legend_bool=true;  # new legend for each state
    for int in 1:Ni
        if legend_bool
            scatter!(ts[int][1:end-1],dζ[int][st,1:Nck[int]],marker = (:pentagon, 10, 0.9, :red),label=string(string("LGR")," || st. #",st),leg=:bottomright)
        else # do not show legend a bunch of times
            scatter!(ts[int][1:end-1],dζ[int][st,1:Nck[int]],marker = (:pentagon, 10, 0.9, :red),label= "",leg=true)
        end
        legend_bool=false;
    end
end
xlims!(t0,tf*ls)
ylabel!("Derivative")
#xlabel!(string("State(x) = ",γ))

tF = zeros(Float64,Ni); yF =  zeros(Float64,Ni);
fp=plot(0,leg=:false);
plot!(t,yS,label="act.",w=lw)
for st in 1:numStates
    for int in 1:Ni
        scatter!(ts[int],PS[int][:,st],markersize =10,markershape = :rect,leg=:topright,label=string("# cp. = ",Nck[int]))
        tF[int] = ts[int][end];
        yF[int] = PS[int][end,st];
    end
end
scatter!(tF,yF,markersize = 10,marker = (:star8, 10, 0.9, :black),label=string("end points"))
xlims!(t0,tf*ls)
ylabel!("State")
xlabel!("x --> really time (t)")

plot(ip,dp,fp,layout=(3,1),background_color_subplot=RGB(0.2,0.2,0.2), background_color_legend=RGB(1,1,1))
plot!(foreground_color_grid=RGB(1,1,1))

#################
# post processing control
#################
ls = 1.35;
lw = 4;

ip=plot(0,leg=:false)
for ctr in 1:numControls
  plot!(t,∫yC[ctr],label=string(@sprintf("act. = %0.3f",∫yC[ctr][end])," || ctr. #",ctr),w=lw);
end
for ctr in 1:numControls
    legend_bool=true;
    for int in 1:Ni
        if legend_bool
            scatter!(ts[int][1:end-1],ctrInt1[int][ctr,1:Nck[int]],marker = (:pentagon, 10, 0.9, :red),leg=:bottomright,label=string(@sprintf("QUAD. = %0.3f", ctrIntVal1[ctr])," || ctr. #",ctr))
            scatter!(ts[int][1:end-1],ctrInt2[int][ctr,1:Nck[int]], marker = (:star5, 10, 0.9, :green),leg=:bottomright,label=string(@sprintf("LGR  = %0.3f", ctrIntVal2[ctr])," || ctr. #",ctr))
        else # do not show legend a bunch of times
            scatter!(ts[int][1:end-1],ctrInt1[int][ctr,1:Nck[int]],marker = (:pentagon, 10, 0.9, :red),leg=:bottomright,label= "",leg=true)
            scatter!(ts[int][1:end-1],ctrInt2[int][ctr,1:Nck[int]], marker = (:star5, 10, 0.9, :green),leg=:bottomright,label= "",leg=true)
        end
        legend_bool=false;
    end
end
xlims!(t0,tf*ls)
ylabel!("Integral")
xlabel!(string("Number of Intervals = ", Ni))

tF = zeros(Float64,Ni); yF =  zeros(Float64,Ni);
fp=plot(0,leg=:false);
plot!(t,yS,label="act.",w=lw)
for ctr in 1:numControls
    for int in 1:Ni
        scatter!(ts[int],PC[int][:,ctr],markersize =10,markershape = :rect,leg=:topright,label=string("# cp. = ",Nck[int]))
    end
end
xlims!(t0,tf*ls)
ylabel!("Control")
xlabel!("x --> really time (t)")

plot(ip,fp,layout=(2,1),background_color_subplot=RGB(0.2,0.2,0.2), background_color_legend=RGB(1,1,1))
plot!(foreground_color_grid=RGB(1,1,1))
