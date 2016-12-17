using NLOptControl
using Polynomials
using Plots
pyplot()
default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

# N point quadrature gives an exact result for polynomials of 2N - 1

# for the first test, we will approximate a polynomial
x₀=-1; # TODO eventually allow for other intervals
xₙ= 1;
x = linspace(x₀,xₙ,100);
α₁ =  3;
α₂ = -3;
α₃ = -8;
α₄ =  7;
γ = Poly([α₁,α₂,α₃,α₄,α₁,α₂,α₃,α₄]);
y = polyval(γ,x)
plot(x,y,label=string(γ),w=4)

# construct polynomial approximation
N = 4;  # N = 2 can exactly represent a 3rd order polynomial
if N < 2 error("Please increase N to at least 2!") end
τ = LGL_nodes(N);
ω = LGL_weights(τ);
D = LGL_Dmatrix(τ);

#= test the derivative approximation
dγ = polyder(γ);
dy = polyval(dγ,x)
plot!(x,dy,label="actual polynomial derivative")
=#
# approximate the integral
∫γ = polyint(γ);
Y = polyval(∫γ,x[end]) - polyval(∫γ,x[1]);
C = Y - polyval(∫γ,x[end]); # find the constant of integration
∫y = polyval(∫γ,x) + C;
plot!(x,∫y,label=@sprintf("integral = %0.3f",∫y[end]),w=4)
#y_integral_approx = [(x[end]-x[1])/2*(w[idx]*polyval(actual_poly_y,tau[idx])) for idx in 1:N+1];
fτ = polyval(γ,τ)
ζ = cumsum(ω.*fτ,1)
plot!(τ,ζ,label=@sprintf("approx. integral = %0.3f",ζ[end]),w=4,size=(800,800))
scatter!(τ,fτ,label=string("tau with N = ",N),markershape = :hexagon, markersize=10)
savefig("test2c.png")
