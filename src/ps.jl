"""
n = DMatrix!(n)
n = DMatrix!(n, (:mode=>:symbolic))
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/15/2017, Last Modified: 1/27/2017 \n
--------------------------------------------------------------------------------------\n
"""
function DMatrix!(n::NLOpt, kwargs...)         #TODO make IMatrix and option
  kw = Dict(kwargs);

  if !haskey(kw,:mode); mode=:defaut;
  else; mode = get(kw,:mode,0);
  end

  if mode==:defaut
    check_ts = zeros(Float64,n.Ni);     # check input
    for int in 1:n.Ni; check_ts[int]=maximum(n.ts[int]); end
    if maximum(check_ts) < 10*eps()
      error("\n ts is full of zeros: make sure that you call create_intervals() first to calculate ts! \n
                NOTE: This may have occured because  (:finalTimeDV => true) and the final time dv is not working properly!! \n")
    end
    D = [zeros((n.Nck[int]+1),(n.Nck[int]+1)) for int in 1:n.Ni];
    for int in 1:n.Ni
        D[int] = polyDiff(n.ts[int]) # +1 is already appended onto ts
    end
    n.DMatrix = [zeros((n.Nck[int]),(n.Nck[int]+1)) for int in 1:n.Ni];
    DM = [zeros((n.Nck[int]),(n.Nck[int])) for int in 1:n.Ni];
    n.IMatrix = [zeros((n.Nck[int]),(n.Nck[int]+1)) for int in 1:n.Ni];
    for int in 1:n.Ni
        n.DMatrix[int] = D[int][1:end-1,:];    # [Nck]X[Nck+1]
        if n.s.integrationScheme==:lgrImplicit
          DM[int] = n.DMatrix[int][:,2:end];   # [Nck]X[Nck]
          n.IMatrix[int] = inv(DM[int]);       # I = inv(D[:,2:N_k+1])
        end
    end
  elseif mode==:symbolic # for validation only, too slow otherwise
    error(" \n cannot precompile with SymPy \n
             so this fucntion was turned off for typical use!! \n
               -> do a (using SymPy) in NLOptControl.jl then remove this error message and rerun \n")
    Dsym = [Array{Any}(n.Nck[int],n.Nck[int]+1) for int in 1:n.Ni];
    n.DMatrix = [Array{Any}(n.Nck[int],n.Nck[int]+1) for int in 1:n.Ni]
    test = [Array{Any}(n.Nck[int]+1) for int in 1:n.Ni];
    val = 1; # since this is always = 1 this funtion is useful for testing, without scaling the problem from [-1,1] this was useful becuase tf was a design variable
    tf = Sym("tf")
    createIntervals!(n,tf); # gives symbolic expression
    for int in 1:n.Ni
      for i in 1:n.Nck[int]+1
        test[int][i] =  n.ts[int][i](tf=>val)
      end
    end
    for int in 1:n.Ni
        for idx in 1:n.Nck[int]+1
            for j in 1:n.Nck[int]
                f = lagrange_basis_poly(tf, test[int], idx)
                Dsym[int][j,idx] = diff(f,tf) # symbolic differentiation --> slow but useful TODO include this in test functions
                n.DMatrix[int][j,idx] = Dsym[int][j,idx](tf=>test[int][j])
            end
        end
    end
  end
  nothing
end

"""
scale_w(w,ta,tb)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 9/18/2017 \n
--------------------------------------------------------------------------------------\n
"""
function scale_w(w,ta,tb)
  (tb - ta)/2*w
end

"""
n = create_intervals(n)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 9/18/2017 \n
--------------------------------------------------------------------------------------\n
"""
function createIntervals!(n::NLOpt)
    tm = linspace(-1,1,n.Ni+1)     # create mesh points
    di = 2/n.Ni                           # interval size
    # go through each mesh interval creating time intervals; map [tm[i-1],tm[i]] --> [-1,1]
    n.ts = [[scale_tau(n.tau[int],tm[int],tm[int+1]);di*int-1] for int in 1:n.Ni]
    n.ws = [scale_w(n.w[int],tm[int],tm[int+1]) for int in 1:n.Ni]
    nothing
end


# NOTE this function was used for testing, but is currently depreciated. When it is used again figure out why and explain why
# di = (tf + 1).n.Ni
function createIntervals!(n::NLOpt, tf)
    tm = linspace(-1,1,n.Ni+1)       # create mesh points
    di = (tf + 1)/n.Ni                      # interval size
    # go through each mesh interval creating time intervals; map [tm[i-1],tm[i]] --> [-1,1]
    n.ts = [[scale_tau(n.tau[int],tm[int],tm[int+1]);di*int-1] for int in 1:n.Ni]
    n.ws = [scale_w(n.w[int],tm[int],tm[int+1]) for int in 1:n.Ni]
    nothing
end

"""
L = lagrange_basis_poly!(x,x_data,L)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 1/25/2016
Citations: This function was influenced by the lagrange() function [located here](https://github.com/pjabardo/Jacobi.jl/blob/master/src/gauss_quad.jl)
--------------------------------------------------------------------------------------\n
# Input Arguments
* `x`: point to approximate function value at
* `x_data`: x data to used calculate basis polynomials
* `L`: array of interest

# Output Arguments
* `L`: Lagrange basis polynomials

A basic description of Lagrange interpolating polynomials is provided [here](http://127.0.0.1:8000/lagrange_poly.html#lagrange-poly)

# TODO get ride of Nc and replace it with Nck[int]
"""
function lagrange_basis_poly!{T<:Number}(x::AbstractArray{T},x_data,L::AbstractArray{T})
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
lagrange_basis_poly{T<:Number}(x::AbstractArray{T},x_data,Nc) = lagrange_basis_poly!(x::AbstractArray{T},x_data,Nc,zeros(x))

"""
D = polyDiff(x);
--------------------------------------------------------------------------\n
Last modifed for julia on 1/25/2016 by Huckleberry Febbo\n
Original Author: JJ.A.C. Weideman, S.C. Reddy 1998\n
Original Function Name: poldif.m  |  Source: [matlabcentral](https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite)\n
https://pdfs.semanticscholar.org/bae2/1eb9458f194887bc8d7808383f56d7f4dca0.pdf
https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
--------------------------------------------------------------------------\n
# Input Arguments
* `x::Vector`: Vector of N distinct nodes.

# Output Arguments
* `D::Array{Float64,2}`: Differention Matrix

"""
function polyDiff(x)  #TODO get ride of B stuff
  N = length(x);
  x = x[:];                     # Make sure x is a column vector
  M = 1;                        # assume
  alpha = ones(N,1);
  B = zeros(M,N);
  L=eye(Bool, N, N);
  XX = repmat(x,1,N);
  DX = XX-XX';                 # DX contains entries x(k)-x(j).
  DX[L] = ones(N,1);           # Put 1's one the main diagonal.
  c = alpha.*prod(DX,2);       # Quantities c(j).
  C = repmat(c,1,N);
  C = C./C';                   # Matrix with entries c(k)/c(j).
  Z = 1./DX;                   # Z contains entries 1/(x(k)-x(j))
  Z[L] = zeros(N,1);           # with zeros on the diagonal.
  X = Z';                      # X is same as Z', but with
  # diagonal entries removed.
  flag = trues(size(X));
  flag[L] = false;
  X = X[flag];
  X = reshape(X,N-1,N);
  Y = ones(N-1,N);             # Initialize Y and D matrices.
  D = eye(N);                  # Y is matrix of cumulative sums,
  temp = reshape(B[1,:],1,N)
  Y   = cumsum([temp; Y[1:N-1,:].*X]);     # Diagonals
  D   = Z.*(C.*repmat(diag(D),1,N) - D);   # Off-diagonals
  D[L]   = Y[N,:];                         # Correct the diagonal
  return D
end
