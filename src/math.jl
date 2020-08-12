
"""
L = lagrange_basis_poly(x0,x,N,j)
--------------------------------------------------------------------------------------
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 1/25/2016
Citations: This function was influenced by the lagrange() function [located here](https://github.com/pjabardo/Jacobi.jl/blob/master/src/gauss_quad.jl)
--------------------------------------------------------------------------------------
# Input Arguments
* `x0`: point to approximate function value at
* `x`: x data to used calculate basis polynomials
* `N`: order of Lagrange interpolating polynomial
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
"""
y = interpolate_lagrange(x::Vector{T}, x_data::Vector{T}, y_data::Vector{T}) where {T <: Number}
"""
function interpolate_lagrange(x, x_data, y_data)


    Nc = length(x_data) - 1

    if length(x_data)!=length(y_data)
        error(string("\n",
                      "-------------------------------------------------------", "\n",
                      "There is an error with the data vector lengths!!", "\n",
                      "-------------------------------------------------------", "\n",
                      "The following variables should be equal:", "\n",
                      "length(x_data) = ",length(x_data),"\n",
                      "length(y_data) = ",length(y_data),"\n"
                      )
              )
      end
      ns = length(x)
      L = zeros(Float64,Nc+1,ns)
      x = x[:]; x_data = x_data[:]; y_data = y_data[:]; # make sure data is in a column
      for idx in 1:Nc+1
        for j in 1:ns
            L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
        end
      end
      return y_data'*L
end

"""
scale_w(w,ta,tb)
--------------------------------------------------------------------------------------
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 9/18/2017
--------------------------------------------------------------------------------------
"""
function scale_w(w::Vector{T},ta::T,tb::T)::Vector{T} where { T <: Number }
    return w .* (tb - ta) ./ 2
end

"""
scale_tau(tau,ta,tb)
--------------------------------------------------------------------------------------
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 12/06/2019
--------------------------------------------------------------------------------------
"""
function scale_tau(tau::Vector{T}, ta::T, tb::T)::Vector{T} where { T <: Number }
    return ( tau .* (tb - ta) .+ (ta + tb) ) ./ 2
end


"""
L = lagrange_basis_poly!(x,x_data,L)
--------------------------------------------------------------------------------------
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 12/02/2019
Last Modified By: Alexander Buck
Citations: This function was influenced by the lagrange() function [located here](https://github.com/pjabardo/Jacobi.jl/blob/master/src/gauss_quad.jl)
--------------------------------------------------------------------------------------
# Input Arguments
* `x`: point to approximate function value at
* `x_data`: x data to used calculate basis polynomials
* `L`: array of interest

# Output Arguments
* `L`: Lagrange basis polynomials

A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

# TODO get ride of Nc and replace it with Nck[int]
"""
function lagrange_basis_poly!(x::AbstractArray{T},x_data,L::AbstractArray{T}) where {T<:Number}
    Nc = length(x_data) - 1
    ns = length(x);
    L = zeros(Float64,Nc+1,ns);
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
      end
    end
    return L
end
lagrange_basis_poly(x::AbstractArray{T},x_data,Nc) where {T<:Number} = lagrange_basis_poly!(x::AbstractArray{T},x_data,Nc,zeros(x))


"""
D = polyDiff(x);
--------------------------------------------------------------------------
Last modifed for julia on 12/06/2019 by Mattias Fält
Original Author: JJ.A.C. Weideman, S.C. Reddy 1998
Original Function Name: poldif.m  |  Source: [matlabcentral](https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite)
https://pdfs.semanticscholar.org/bae2/1eb9458f194887bc8d7808383f56d7f4dca0.pdf
https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
--------------------------------------------------------------------------
# Input Arguments
* `x::Vector`: Vector of N distinct nodes.

# Output Arguments
* `D::Array{Float64,2}`: Differention Matrix

"""
function polyDiff(x)  #TODO get rid of B stuff

    N = length(x)

    x = x[:]

    M = 1

    alpha = ones(N, 1)

    B = zeros(M, N)

    L = Matrix{Bool}(LinearAlgebra.I, N, N)

    XX = repeat(x, 1, N)

    # DX contains entries x(k)-x(j).
    DX = XX - XX'

    # Put 1's one the main diagonal.
    DX[L] = ones(N, 1)

    # Quantities c(j)
    c = alpha .* prod(DX, dims=2)

    C = repeat(c, 1, N)

    # Matrix with entries c(k)/c(j).
    C = C ./ C'

    # Z contains entries 1/(x(k)-x(j))
    Z = 1 ./ DX

    # with zeros on the diagonal.
    Z[L] = zeros(N, 1)

    # X is same as Z', but with
    X = Z'

    # diagonal entries removed.
    flag = trues(size(X));
    flag[L] .= false;
    X = X[flag];

    X = reshape(X, N-1, N)

    # Initialize Y and D matrices.
    Y = ones(N-1, N)

    # Y is matrix of cumulative sums,
    D = Matrix{Float64}(LinearAlgebra.I, N, N)

    temp = reshape(B[1,:], 1, N)

    # Diagonals
    Y = cumsum([ temp ; Y[1:N-1,:] .* X ], dims=1)

    # Off-diagonals
    D = Z .* (C .* repeat(LinearAlgebra.diag(D),1,N) - D)

    # Correct the diagonal
    D[L] = Y[N, :]

    return D

end


"""
# t must be the time vector
# V is any vector that you are interpolating
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 2/06/2017, Last Modified: 12/06/2019 \n
--------------------------------------------------------------------------------------\n
"""
function linearSpline(t::Vector, V::Vector)
    if !isequal(length(t),length(V))
        error("!isequal(length(t),length(V))")
    end
    # remove any repeating values
    M = Array{Bool}(undef, length(t)); M[1:length(t)] .= false;
    for i in 1:length(t)-1
        if t[i]==t[i+1]
            M[i]=true
        else
            M[i]=false
        end
    end
    rm_idx = findall(M)
    if (length(t)==length(rm_idx))
        error("No time has elapsed and there will be an issue with interpolation. \n
            Cannot simulate the vehicle.")
    end

    # initialize vetors
    t_new = Array{Float64}(undef, length(t)-length(rm_idx))
    V_new = Array{Float64}(undef, length(t)-length(rm_idx))
    q = 1
    for i in 1:length(V) #TODO put an error message here if V and t are different sizes
        if !M[i]
            t_new[q] = t[i]
            V_new[q] = V[i]
            q = q + 1
        end
    end

    # make interpolant using Dierckx.jl
    #Spline1D(t_new,V_new,k=1)    # linear spline

    # make interpolant using Interpolations.jl
    knots = (t_new,)
    it = interpolate(knots,V_new,Gridded(Linear()))
    et = extrapolate(it, Flat())

    return et
end
