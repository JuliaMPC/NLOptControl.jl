using NLOptControl
using Polynomials
using Plots
using FastGaussQuadrature
pyplot()
#default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

# define problem
t0 = Float64(0); tf = Float64(10);  # TODO change and y to x
t = Array(linspace(t0,tf,100));
α₁ =  3; α₂ = -3; α₃ = -8; α₄ =  7;
Nc = Int64(10); # number of collocation points in each interval
Ni = Int64(4); # number of intervals
####################################
# perform analytical calcualtions
####################################
γ = Poly([α₁,α₂,α₁,α₂,α₁,α₂]); #TODO check on that imported binding warning
y = polyval(γ,t)

# evaluate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,t[end]) - polyval(∫γ,t[1]);
C = Y - polyval(∫γ,t[end]); # constant of integration
∫y = polyval(∫γ,t) + C;

# evaluate the derivative
dγ = polyder(γ);
dy = polyval(dγ,t);

####################################
# construct polynomial approximation
####################################
τ, ω = gaussradau(Nc); # number of collocation points per interval

# break the problem up into multiple intervals
di, tm, ts, ωₛ = create_intervals(t0,tf,Ni,Nc,τ,ω)

#TODO check integral matrix

# approximate the integral
ζ = zeros(Float64,Nc,Ni); fτ = zeros(Float64,Nc,Ni);approx_int = Float64(0);
for idx in 1:Ni
  fτ[:,idx] = polyval(γ,ts[:,idx]);
  ζ[:,idx] =  cumsum(ωₛ[:,idx].*fτ[:,idx],1)
  approx_int = approx_int + ζ[end,idx];
end
# approximate the derivative --> needed in defect constraints
dζ = zeros(Float64,Nc,Ni);
for idx in 1:Ni
  D = poldif(ts[:,idx], 1)
  dζ[:,idx] = D*fτ[:,idx]
end

#################
# post processing
#################
fp=plot(0,leg=:false)
plot!(t,y,label=string(γ),w=6)
scatter!(ts,fτ,label=string("ftau with Ni = ",Ni),markershape = :hexagon, markersize=10)

dp=plot(0,leg=:false)
plot!(t,dy,label="derivative",w=6)
plot!(ts,dζ,label="approx. derivative",line=(4,:dash))

ip=plot(0,leg=:false)
plot!(t,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=6)
plot!(ts,ζ,label=@sprintf("approx. integral = %0.3f",approx_int),line=(4,:dash))

plot(fp,ip,dp,layout=(3,1))

#line=(4,:dash,0.6,[:lightgreen :green :darkgreen]),background_color=RGB(0.2,0.2,0.2)
