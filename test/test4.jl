using NLOptControl
using Polynomials
using Plots
pyplot()
default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

x₀= Float64(-1);
xₙ= Float64(1);
x = Array(linspace(x₀,xₙ,100));
α₁ =  3;
α₂ = -3;
α₃ = -8;
α₄ =  7;
γ = Poly([α₁,α₂,α₁]);
y = polyval(γ,x)
plot(x,y,label=string(γ),w=4)

# evaluate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,x[end]) - polyval(∫γ,x[1]);
C = Y - polyval(∫γ,x[end]); # find the constant of integration
∫y = polyval(∫γ,x) + C;
plot!(x,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=4)

# construct polynomial approximation
N = 3;  # N = 2 can exactly represent a 3rd order polynomial
if N < 2 error("Please increase N to at least 2!") end
τ = LGL_nodes(N);
ω = LGL_weights(τ);
D = LGL_Dmatrix(τ);

# scale the problem to the interval
xₛ = scale_tau(τ,x₀,xₙ) # scale the interval
ωₛ = scale_w(ω,x₀,xₙ)   # scale the weights
fτ = polyval(γ,xₛ)
ζ =  cumsum(ω.*fτ,1)
plot!(τ,ζ,label=@sprintf("approx. integral = %0.3f",ζ[end]),w=4,size=(800,800))
scatter!(τ,fτ,label=string("tau with N = ",N),markershape = :hexagon, markersize=10)

# approximate the derivative
dγ = polyder(γ); # actual derivative
dy = polyval(dγ,x)
plot!(x,dy,label="derivative",w=5)
D = LGL_Dmatrix(s);
dζ = D*fτ;       # approximate derivative
plot!(τ, dζ, label="approx. derivative",line=(4,:dash))
savefig("test4c.png")
pe = (∫y[end]-ζ[end])/∫y[end]*100;
@sprintf("The percent error is = %0.2f", pe)
