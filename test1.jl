using NLOptControl
using Polynomials
using Plots
using FastGaussQuadrature
#using Jacobi
pyplot()
#default(guidefont = font(17), tickfont = font(15), legendfont = font(12), titlefont = font(20))

# define problem
t0 = Float64(0); tf = Float64(10);  # TODO change and y to x
t = Array(linspace(t0,tf,100));
α₁ =  3; α₂ = -3; α₃ = -8; α₄ =  7;
Nc = Int64(5); # number of collocation points in each interval
Ni = Int64(1); # number of intervals
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


"""
L = lagrange_basis_poly(x::Float64,N::Int64)
--------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 12/26/2016
--------------------------------------------------------------------\n
# Input Arguments
* `x::Float64`: point to approximate function value at
* `x_data::Array{Float64}`: x data to used calculate basis polynomials
* `N::Int64`: order of Lagrange interpolating polynomial

# Output Arguments
* `L::Array{Float64}`: Lagrange basis polynomials

A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

"""
function lagrange_basis_poly{T<:Number}(x::Float64,x_data::AbstractArray{T},idx::Int64,N::Int64)
    L = 1;
    for j in 1:N
      if j!=idx
        L = L*(x - x_data[j])/(x_data[idx]-x_data[j]);
      end
    end
  return L
end

# example, interpolate f(x) = x^2 over 1<=x<=3 given
x_data = [1,2,3];
y_data = [1,4,9];
x0 = 1; xf = 3;
ns = 100;  # plotting points
x = Array(linspace(x0,xf,ns));
y = x.^2;
plot(x,y, label= "actual polynomial")

N = 2; # order of Lagrange Polynomial
if N > length(x_data)
  error("Maximum N value = length(x_data)")
elseif N > length(x_data)-1
  warn("Reduce N to = length(x_data) - 1")
end

L = zeros(Float64,N,ns);
for idx in 1:N
  for j in 1:ns
    L[idx,j] = lagrange_basis_poly(x[j],x_data,idx,N)
  end
  plot!(x,L[idx,:]) # plot the Lagrange basis polynomials
end

plot!(x,L[1,:])
plot!(x,L[2,:])
#plot!(x,L[3,:])
