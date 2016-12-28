function scale_tau(τ::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*τ + (xₙ + x₀)/2
end

function scale_w(ω::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*ω
end

function create_intervals(t0::Float64,tf::Float64,Ni::Int64,Nc::Int64,τ::Vector{Float64},ω::Vector{Float64})
  di = (tf - t0)/Ni # interval size
  # create mesh points
  tm = zeros(Float64,Ni+1); tm[1] = t0;
  for idx in 1:Ni
    tm[idx+1] = tm[idx] + di;
  end
  # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
  ts = zeros(Float64,Nc,Ni);  ωₛ = zeros(Float64,Nc,Ni);
  for idx in 1:Ni
    # scale the problem to the interval
    ts[:,idx] = scale_tau(τ,tm[idx],tm[idx+1]); # scale the interval
    ωₛ[:,idx] = scale_w(ω,tm[idx],tm[idx+1]);   # scale the weights
  end
  return di, tm, ts, ωₛ
end

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
    L = zeros(Float64,Nc+1,ns);
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx)
      end
    end
    return L
end
lagrange_basis_poly{T<:Number}(x::AbstractArray{T},x_data,Nc) = lagrange_basis_poly!(x::AbstractArray{T},x_data,Nc,zeros(x))

function interpolate_lagrange{T<:Number}(x::AbstractArray{T},x_data,y_data,Nc)
    if Nc > length(x_data) -1
      error("Maximum N value = length(x_data)-1")
    end
    ns = length(x);
    L = zeros(Float64,Nc+1,ns);
    x = x[:]; x_data = x_data[:]; y_data = y_data[:]; # make sure data is in a column
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx)
      end
    end
    y = y_data'*L;
    return y
end

"""
D = poldif(x, malpha, B...);
--------------------------------------------------------------------------\n
Last modifed for julia on December 23, 2016 by Huckleberry Febbo\n
Original Author: JJ.A.C. Weideman, S.C. Reddy 1998\n
Original Function Name: poldif.m  |  Source: [matlabcentral](https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite)\n
--------------------------------------------------------------------------\n
# Input Arguments
* `x::Vector`: Vector of N distinct nodes.
* `malpha::Int65`: M, the number of derivatives required
    * Note:     0 < M < N-1
*  `B`: Currently this functionality is not tested with this input (only works for constant weights)

# Output Arguments
* `D::Array{Float64,2}`: Differention Matrix
    * DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M

The function DM =  poldif(x, maplha, B) computes the differentiation matrices D1, D2, ..., DM on arbitrary nodes. The function is called with either two or three input arguments. If two input arguments are supplied, the weight function is assumed to be constant.   If three arguments are supplied, the weights should be defined as the second and third arguments. (CURRENTLY NOT TESTED IN julia)

"""

function poldif(x, malpha, B...)
  B_given = Bool(length(B));
         N = length(x);
         x = x[:];                     # Make sure x is a column vector

  if !B_given                       # Check if constant weight function
         M = malpha;                   # is to be assumed.
     alpha = ones(N,1);
         B = zeros(M,N);
  elseif B_given
     alpha = malpha(:);                # Make sure alpha is a column vector
         M = length(B[:,1]);           # First dimension of B is the number
  end                                  # of derivative matrices to be computed

        #  I = eye(N);                  # Identity matrix.
          #L = logical(I);              # Logical identity matrix.
          L=eye(Bool, N, N)
      # XX = x(:,ones(1,N));
         XX = repmat(x,1,N);
         DX = XX-XX';                  # DX contains entries x(k)-x(j).

      DX[L] = ones(N,1);               # Put 1's one the main diagonal.

          c = alpha.*prod(DX,2);       # Quantities c(j).

          #C = c[:,ones(1,N)];
          C = repmat(c,1,N);
          C = C./C';                   # Matrix with entries c(k)/c(j).

          Z = 1./DX;                   # Z contains entries 1/(x(k)-x(j))
       Z[L] = zeros(N,1);              # with zeros on the diagonal.

          X = Z';                      # X is same as Z', but with
       #X[L] = [];                      # diagonal entries removed.
       flag = trues(size(X));
    flag[L] = false;
          X = X[flag];
          X = reshape(X,N-1,N);

          Y = ones(N-1,N);             # Initialize Y and D matrices.
          D = eye(N);                  # Y is matrix of cumulative sums,

  DM = zeros(Float64,N,N,M);                                   # D differentiation matrices.
  for ell = 1:M
          temp = reshape(B[ell,:],1,N)
          Y   = cumsum([temp; ell*Y[1:N-1,:].*X]); # Diagonals
          D   = ell*Z.*(C.*repmat(diag(D),1,N) - D);   # Off-diagonals
       D[L]   = Y[N,:];                                # Correct the diagonal
  DM[:,:,ell] = D;                                     # Store the current D

  end

  return squeeze(DM,3)
end

# the following code uses the function lepoly() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
function lepoly(n::Int64,x,derivative_option::Bool)
# lepoly  Legendre polynomial of degree n
    # y=lepoly(n,x) is the Legendre polynomial
    # The degree should be a nonnegative integer
    # The argument x should be on the closed interval [-1,1];
    # [dy,y]=lepoly(n,x) also returns the values of 1st-order
    #  derivative of the Legendre polynomial stored in dy
# Last modified on August 30, 2011
# Verified with the chart in http://keisan.casio.com/has10/SpecExec.cgi
  if !derivative_option
       if n==0 y=ones(size(x));  return y end
       if n==1 y=x; return y end
       polylst=ones(size(x)); poly=x;   # L_0(x)=1, L_1(x)=x
       polyn=0;
       for k=2:n                      # Three-term recurrence relation:
  	     polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k; # kL_k(x)=(2k-1)xL_{k-1}(x)-(k-1)L_{k-2}(x)
         polylst=poly; poly=polyn;
       end
       y=polyn;
  else
     if n==0 y=ones(size(x)); dy=zeros(size(x)); return dy,y end
     if n==1 y=x; dy=ones(size(x)); return dy, y end

     # TODO MAKE SURE THIS IS CALCULATING SIZE correctly!! i.e. N=4 and size(N+1,N+1) = 1
      polylst=ones(size(x)); pderlst=zeros(size(x));poly=x; pder=ones(size(x));
       # L_0=1, L_0'=0, L_1=x, L_1'=1
      polyn=0;
      pdern=0;
      for k=2:n                          # Three-term recurrence relation:
        polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k; # kL_k(x)=(2k-1)xL_{k-1}(x)-(k-1)L_{k-2}(x)
        pdern=pderlst+(2*k-1)*poly; # L_k'(x)=L_{k-2}'(x)+(2k-1)L_{k-1}(x)
     	  polylst=poly; poly=polyn;
    	  pderlst=pder; pder=pdern;
      end
      y=polyn; dy=pdern;
      return dy, y
  end

end
