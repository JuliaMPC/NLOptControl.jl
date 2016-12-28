using ForwardDiff

"""
L = lagrange_basis_poly(x,x_data,Nc,j)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 12/28/2016
Citations: This function was influenced by the lagrange() function [located here](https://github.com/pjabardo/Jacobi.jl/blob/master/src/gauss_quad.jl)
--------------------------------------------------------------------------------------\n
# Input Arguments
* `x`: point to approximate function value at
* `x_data`: x data to used calculate basis polynomials
* `Nc`: order of Lagrange interpolating polynomial
* `j`: index of interest

# Output Arguments
* `L`: Lagrange basis polynomials

A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

"""
function lagrange_basis_poly(x,x_data,Nc,j)
    L = 1;
    for idx in 1:Nc+1 # use all of the data
      if idx!=j
        L = L*(x - x_data[idx])/(x_data[j]-x_data[idx]);
      end
    end
  return L
end

function lagrange_basis_poly!{T<:Number}(x::AbstractArray{T},x_data,Nc,L::AbstractArray{T})
   if Nc > length(x_data) -1
      error("Maximum N value = length(x_data)-1")
    end
    ns = length(x);
    L = zeros(eltype(x),Nc+1,ns);
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx)
      end
    end
    return L
end
lagrange_basis_poly{T<:Number}(x::AbstractArray{T},x_data,Nc) = lagrange_basis_poly!(x::AbstractArray{T},x_data,Nc,zeros(x))

function lagrange_basis_poly_vec(x, x_data, Nc)
   if Nc > length(x_data) -1
      error("Maximum N value = length(x_data)-1")
    end
    ns = length(x);
    L = zeros(eltype(x),Nc+1,ns);
    for idx in 1:Nc+1
      for j in 1:ns
        L_temp = Float64(1);
        for idx2 in 1:Nc+1 # use all of the data
          if idx!=idx2
            L_temp = L_temp*(x[j] - x_data[idx2])/(x_data[idx]-x_data[idx2]);
          end
        end
        L[idx,j] = L_temp;
      end
    end
    return L
end

# Derivative of Each Element
Nc = 2; x_data = [1,2,3]; j=1; x = 1;
d = ForwardDiff.derivative(lagrange_basis_poly, x, x_data, Nc, j)

function rosenbrock(x::AbstractVector, a, b)
    result = zero(eltype(x))
    for i in 1:length(x)-1
        result += (a - x[i])^2 + b*(x[i+1] - x[i]^2)^2
    end
    return result
end

x = [1,2];c = length(x);a = one(eltype(x));b = 100 * a;
cfg = ForwardDiff.GradientConfig{c}(x)
d = ForwardDiff.gradient(rosenbrock, x, cfg, a, b)

# Gradient of Lagrange Polynomials
x = [1,2]; Nc = 2; x_data = [1,2,3];
d = ForwardDiff.gradient(lagrange_basis_poly, x, cfg, x_data, Nc)

# Gradient of Lagrange Polynomials
Nc = 2; x_data = [1,2,3];
d = ForwardDiff.gradient(lagrange_basis_poly_vec, x, cfg, x_data, Nc)

#=

function lagrange_basis_poly_vec(x,x_data,Nc,L)
   if Nc > length(x_data) -1
      error("Maximum N value = length(x_data)-1")
    end
    ns = length(x);
    #L = zeros(Nc+1,ns);
    for idx in 1:Nc+1
      for j in 1:ns
        L_temp = Float64(1);
        for idx2 in 1:Nc+1 # use all of the data
          if idx!=idx2
            L_temp = L_temp*(x[j] - x_data[idx2])/(x_data[idx]-x_data[idx2]);
          end
        end
        print(L_temp)
        print(ns)
        print(j)
        print(x)
        print(L)
        L[idx,j] = L_temp;
      end
    end
    return L
end
function lagrange_basis_poly(x,x_data,Nc,j)
    L = 1;
    for idx in 1:Nc+1 # use all of the data
      if idx!=j
        L = L*(x - x_data[idx])/(x_data[j]-x_data[idx]);
      end
    end
  return L
end

=#
