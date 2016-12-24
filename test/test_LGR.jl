using NLOptControl
using Polynomials
using Plots
using FastGaussQuadrature

pyplot()
default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

# define problem
x₀= Float64(-1); xₙ= Float64(1);  # TODO change x to t and y to x
x = Array(linspace(x₀,xₙ,100));
α₁ =  3; α₂ = -3; α₃ = -8; α₄ =  7;
N = Int64(4); # TODO eventually make it multiple-interval

####################################
# perform analytical calcualtions
####################################
γ = Poly([α₁,α₂,α₁,α₁]); #TODO check on that imported binding warning
# 1/(25x^2+1)
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
τ,ω,I,D = LGR(N);

# scale the problem to the interval
xₛ = scale_tau(τ,x₀,xₙ) # scale the interval
ωₛ = scale_w(ω,x₀,xₙ)   # scale the weights
#TODO check integral matrix
# TODO figure out how to get the last part of the integral!

# approximate the integral
fτ = polyval(γ,xₛ)
ζ =  cumsum(ωₛ.*fτ,1)

# approximate the derivative --> needed in defect constraints
dζ = D*fτ;
dζ = dζ[1:N+1];

#################
# post processing
#################
# actual
plot(x,y,label=string(γ),w=6)
plot!(x,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=6)
plot!(x,dy,label="derivative",w=6)
# approximate
plot!(xₛ,ζ,label=@sprintf("approx. integral = %0.3f",ζ[end]),line=(4,:dash),size=(800,800))
scatter!(xₛ,fτ,label=string("ftau with N = ",N),markershape = :hexagon, markersize=10)
plot!(xₛ,dζ,label="approx. derivative",line=(4,:dash))
