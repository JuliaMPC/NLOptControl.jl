using NLOptControl
using Polynomials
using Plots
pyplot()
default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

x₀= Float64(0);
xₙ= Float64(10);
x = Array(linspace(x₀,xₙ,100));
α₁ =  3;
α₂ = -3;
α₃ = -8;
α₄ =  7;
γ = Poly([α₁,α₂]);
y = polyval(γ,x)
plot(x,y,label=string(γ),w=4)

# evaluate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,x[end]) - polyval(∫γ,x[1]);
C = Y - polyval(∫γ,x[end]); # find the constant of integration
∫y = polyval(∫γ,x) + C;
plot!(x,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=4)

# construct polynomial approximation
N = 4;  # N = 2 can exactly represent a 3rd order polynomial
if N < 2 error("Please increase N to at least 2!") end
τ = LGL_nodes(N);
ω = LGL_weights(τ);
D = LGL_Dmatrix(τ);

# scale the problem to the interval --> [-1,1]
xₛ = scale_X(x,x₀,xₙ,:quad_scale)
s = (xₙ - x₀)/2*τ + (xₙ + x₀)/2
fτ = polyval(γ,s)
ζ = (xₙ - x₀)/2*cumsum(ω.*fτ,1)
plot!(τ,ζ,label=@sprintf("approx. integral = %0.3f",ζ[end]),w=4,size=(800,800))
scatter!(τ,fτ,label=string("tau with N = ",N),markershape = :hexagon, markersize=10)
savefig("test3b.png")
