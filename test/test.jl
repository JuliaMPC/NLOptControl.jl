using NLOptControl
using Polynomials
using Plots
pyplot()
default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

# define problem
x₀= Float64(-1); xₙ= Float64(1);  # TODO change x to t and y to x
x = Array(linspace(x₀,xₙ,100));
α₁ =  3; α₂ = -3; α₃ = -8; α₄ =  7;
N = Int64(3); # TODO eventually make it multiple-interval

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
# TODO eventually check to see which scheme we are using, LGL etc.
τ = LGL_nodes(N);    # calculate nodes
ω = LGL_weights(τ);  # calculate weights
D = LGL_Dmatrix(τ);  # calculate differention matrix

# scale the problem to the interval
xₛ = scale_tau(τ,x₀,xₙ) # scale the interval
ωₛ = scale_w(ω,x₀,xₙ)   # scale the weights

# approximate the integral --> needed in cost function
fτ = polyval(γ,xₛ)
ζ =  cumsum(ωₛ.*fτ,1)

# approximate the derivative --> needed in defect constraints
#D = LGL_Dmatrix(xₛ);  # calculate differention matrix

# From a modified version of what is described in: 
# "Legendre–Gauss–Lobatto Pseudo–spectral Method for One–Dimensional Advection–Diffusion Equation"
function MYD(τ::Vector{Float64},N::Int64)
  D = zeros(N+1,N+1)
  for i in 1:N+1
    for j in 1:N+1
      if i!=j
        D[i,j] = lepoly(N,τ[i],false)[1]/(lepoly(N,τ[j],false)[1]*(τ[i]-τ[j]));
      elseif (i==j && j==1)
        D[i,j] = -N*(N+1)/4;
      elseif (i==j && j==(N+1)) #TODO added the N+1 --> check this
        D[i,j] = N*(N+1)/4;
      else
        D[i,j] = 0;
      end
    end
  end
  return D
end
#=
temp =2*(xₛ - x₀)/(xₙ - x₀) -1;
D = MYD(temp,N)
temp2 = (xₙ - x₀)*(τ+1)/2+x₀
fτ = polyval(γ,temp2)
dζ = ωₛ.*((2/(xₙ - x₀))*fτ'*D)'; # to get high-order derivatives just raise the power of D and (2/(xₙ - x₀))
=#

#=MINE THAt works FOR [-1,1]
D=MYD(xₛ,N)
dζ =  D*fτ
=#######
D=MYD(τ,N)
fτ = polyval(γ,xₛ)
# dζ = ωₛ.*((2/(xₙ - x₀))*fτ'*D)'; # to get high-order derivatives just raise the power of D and (2/(xₙ - x₀))
dζ = D*fτ;

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

#ylims!((-10,10))

savefig("test4d.png")
@sprintf("The percent error in the integral approximation is = %0.2f", (∫y[end]-ζ[end])/∫y[end]*100)
