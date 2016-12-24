using NLOptControl
using Polynomials
using Plots
using FastGaussQuadrature
pyplot()
#default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

# define problem
x₀= Float64(-1); xₙ= Float64(20);  # TODO change x to t and y to x
x = Array(linspace(x₀,xₙ,100));
α₁ =  3; α₂ = -3; α₃ = -8; α₄ =  7;
N = Int64(4); # TODO eventually make it multiple-interval

####################################
# perform analytical calcualtions
####################################
γ = Poly([α₁,α₂,α₁,α₂]); #TODO check on that imported binding warning
# 1/(25x^2+1) TODO check the runge phenomina
y = polyval(γ,x)

# evaluate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,x[end]) - polyval(∫γ,x[1]);
C = Y - polyval(∫γ,x[end]); # constant of integration
∫y = polyval(∫γ,x) + C;

# evaluate the derivative
dγ = polyder(γ);
dy = polyval(dγ,x);

####################################
# construct polynomial approximation
####################################
τ, ω = gaussradau(N);

# scale the problem to the interval
xs = scale_tau(τ,x₀,xₙ) # scale the interval
ωₛ = scale_w(ω,x₀,xₙ)   # scale the weights
#TODO check integral matrix
#TODO figure out how to get the last part of the integral!

# approximate the integral
fτ = polyval(γ,xs);
ζ =  cumsum(ωₛ.*fτ,1)

# approximate the derivative --> needed in defect constraints
D = poldif(xs, 1)
dζ = D*fτ

#################
# post processing
#################
fp=plot(0,leg=:false)
plot!(x,y,label=string(γ),w=6)
scatter!(xs,fτ,label=string("ftau with N = ",N),markershape = :hexagon, markersize=10)

dp=plot(0,leg=:false)
plot!(x,dy,label="derivative",w=6)
plot!(xs,dζ,label="approx. derivative",line=(4,:dash))

ip=plot(0,leg=:false)
plot!(x,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=6)
plot!(xs,ζ,label=@sprintf("approx. integral = %0.3f",ζ[end]),line=(4,:dash))

plot(fp,ip,dp,layout=(3,1))

#line=(4,:dash,0.6,[:lightgreen :green :darkgreen]),background_color=RGB(0.2,0.2,0.2)
