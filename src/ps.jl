"""
n = DMatrix!(n)
n = DMatrix!(n, (:mode=>:symbolic))
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/15/2017, Last Modified: 12/06/2019 \n
--------------------------------------------------------------------------------------\n
"""
function DMatrix!(n::NLOpt, kwargs...)         #TODO make IMatrix and option

    kw = Dict(kwargs)

    if !haskey(kw,:mode)
        mode = :defaut
    else
        mode = get(kw,:mode,0)
    end

    if mode == :defaut

        check_ts = maximum.(n.ocp.ts)

        if maximum(check_ts) < 10 * eps()
            error("""
                ts is full of zeros: make sure that you call create_intervals() first to calculate ts!
                NOTE: This may have occured because  (:finalTimeDV => true) and the final time dv is not working properly!!
                """)
        end

        D = [ zeros((n.ocp.Nck[int]+1),(n.ocp.Nck[int]+1)) for int in 1:n.ocp.Ni ]

        for int in 1:n.ocp.Ni
            D[int] = polyDiff(n.ocp.ts[int]) # +1 is already appended onto ts
        end

        n.ocp.DMatrix = [zeros((n.ocp.Nck[int]),(n.ocp.Nck[int]+1)) for int in 1:n.ocp.Ni]
        DM = [zeros((n.ocp.Nck[int]),(n.ocp.Nck[int])) for int in 1:n.ocp.Ni]
        n.ocp.IMatrix = [zeros((n.ocp.Nck[int]),(n.ocp.Nck[int]+1)) for int in 1:n.ocp.Ni]

        for int in 1:n.ocp.Ni
            n.ocp.DMatrix[int] = D[int][1:end-1,:]      # [Nck]X[Nck+1]
            if n.s.ocp.integrationScheme == :lgrImplicit
                DM[int] = n.ocp.DMatrix[int][:,2:end]   # [Nck]X[Nck]
                n.ocp.IMatrix[int] = inv(DM[int])       # I = inv(D[:,2:N_k+1])
            end
        end

    elseif mode == :symbolic # for validation only, too slow otherwise
        error("""
            Cannot precompile with SymPy
            so this fucntion was turned off for typical use!!
            -> do a (using SymPy) in NLOptControl.jl then remove this error message and rerun
            """)
        Dsym = [Array{Any}(undef, n.ocp.Nck[int],n.ocp.Nck[int]+1) for int in 1:n.ocp.Ni];
        n.ocp.DMatrix = [Array{Any}(undef, n.ocp.Nck[int],n.ocp.Nck[int]+1) for int in 1:n.ocp.Ni]
        test = [Array{Any}(undef, n.ocp.Nck[int]+1) for int in 1:n.ocp.Ni];
        val = 1; # since this is always = 1 this funtion is useful for testing, without scaling the problem from [-1,1] this was useful becuase tf was a design variable
        tf = Sym("tf")
        createIntervals!(n, tf); # gives symbolic expression
        for int in 1:n.ocp.Ni
            for i in 1:n.ocp.Nck[int]+1
                test[int][i] =  n.ocp.ts[int][i](tf=>val)
            end
        end
        for int in 1:n.ocp.Ni
            for idx in 1:n.ocp.Nck[int]+1
                for j in 1:n.ocp.Nck[int]
                    f = lagrange_basis_poly(tf, test[int], idx)
                    Dsym[int][j,idx] = diff(f,tf) # symbolic differentiation --> slow but useful # TODO include this in test functions
                    n.ocp.DMatrix[int][j,idx] = Dsym[int][j,idx](tf=>test[int][j])
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
Date Create: 12/23/2017, Last Modified: 12/06/2019 \n
--------------------------------------------------------------------------------------\n
"""
function createIntervals!(n::NLOpt)
    tm = range(-1,1;length=n.ocp.Ni+1)     # create mesh points
    di = 2/n.ocp.Ni                           # interval size
    # go through each mesh interval creating time intervals; map [tm[i-1],tm[i]] --> [-1,1]
    n.ocp.ts = hcat([ [ scale_tau(n.ocp.tau[int],tm[int],tm[int+1]) ; di*int-1 ] for int in 1:n.ocp.Ni ]...)
    n.ocp.ws = hcat([ scale_w(n.ocp.w[int],tm[int],tm[int+1]) for int in 1:n.ocp.Ni ]...)
    return nothing
end


# NOTE this function was used for testing, but is currently depreciated. When it is used again figure out why and explain why
# di = (tf + 1).n.ocp.Ni
function createIntervals!(n::NLOpt, tf)
    tm = range(-1,1; length=n.ocp.Ni+1)       # create mesh points
    di = (tf + 1)/n.ocp.Ni                      # interval size
    # go through each mesh interval creating time intervals; map [tm[i-1],tm[i]] --> [-1,1]
    n.ocp.ts = hcat([ [scale_tau(n.ocp.tau[int],tm[int],tm[int+1]) ; di*int-1 ] for int in 1:n.ocp.Ni]...)
    n.ocp.ws = hcat([ scale_w(n.ocp.w[int],tm[int],tm[int+1]) for int in 1:n.ocp.Ni]...)
    return nothing
end

"""
L = lagrange_basis_poly!(x,x_data,L)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/26/2016, Last Modified: 12/02/2019
Last Modified By: Alexander Buck
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
--------------------------------------------------------------------------\n
Last modifed for julia on 12/06/2019 by Mattias Fält\n
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
    x = vcat(x...);                     # Make sure x is a column vector
  M = 1;                        # assume
  alpha = ones(N,1);
  B = zeros(M,N);
  L= Matrix{Bool}(LinearAlgebra.I, N, N);
  XX = repeat(x,1,N);
  DX = XX-XX';                 # DX contains entries x(k)-x(j).
  DX[L] = ones(N,1);           # Put 1's one the main diagonal.
  c = alpha.*prod(DX,dims=2);       # Quantities c(j).
  C = repeat(c,1,N);
  C = C./C';                   # Matrix with entries c(k)/c(j).
  Z = 1 ./ DX;                   # Z contains entries 1/(x(k)-x(j))
  Z[L] = zeros(N,1);           # with zeros on the diagonal.
  X = Z';                      # X is same as Z', but with
  # diagonal entries removed.
  flag = trues(size(X));
  flag[L] .= false;
  X = X[flag];
  X = reshape(X,N-1,N);
  Y = ones(N-1,N);             # Initialize Y and D matrices.
  D = Matrix{Float64}(LinearAlgebra.I,N,N);                  # Y is matrix of cumulative sums,
  temp = reshape(B[1,:],1,N)
  Y   = cumsum([temp; Y[1:N-1,:].*X], dims=1);     # Diagonals
  D   = Z.*(C.*repeat(LinearAlgebra.diag(D),1,N) - D);   # Off-diagonals
  D[L]   = Y[N,:];                         # Correct the diagonal
  return D
end
