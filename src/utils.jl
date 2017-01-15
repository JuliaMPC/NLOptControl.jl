"""
F_matrix(nlp,ps,int)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/7/2017, Last Modified: 1/12/2017 \n
--------------------------------------------------------------------------------------\n
"""
function F_matrix(nlp::NLP_data, ps::PS_data)
    @unpack FMatrix  = ps
    @unpack stateEquations = nlp

    # TODO if stateEquations() is not defined warn the user
    FMatrix = stateEquations(nlp,ps)

    print(FMatrix)
    @pack ps = FMatrix
end

"""

integral = integrate(mdl,ps,u;(:mode=>:control))

# integrates everything
stInt, stIntVal, ctrInt, ctrIntVal = integrate(ps,nlp;(:mode=>:LGRIM))
stInt, stIntVal, ctrInt, ctrIntVal = integrate(ps,nlp)
# integrates specific states and or controls
  * TODO have an option to integrate only certain states and controls
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/2/2017, Last Modified: 1/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function integrate(ps::PS_data,nlp::NLP_data; kwargs...)
    @unpack Nck, Ni, stateMatrix, controlMatrix, ωₛ, IMatrix = ps
    @unpack numStates, numControls = nlp

    kw = Dict(kwargs);
    # if there was nothing passed -> set default
    if !haskey(kw,:mode); kw = Dict(:mode => :default) end
    mode = get(kw, :mode, 0);

    # state integral
    stInt = [zeros(Float64, numStates, Nck[int],)for int in 1:Ni]; idx = 1;
    stTemp = zeros(Float64, numStates, Ni);
    for st in 1:numStates
        for int in 1:Ni
            if mode == :default
                stInt[int][st,1:Nck[int]] =  cumsum(ωₛ[int].*stateMatrix[int][1:end-1,st]);
            elseif mode == :LGRIM     # Legendre-Gauss-Radau integration matrix (LGRIM)
                stInt[int][st,1:Nck[int]] =  IMatrix[int]*stateMatrix[int][1:end-1,st];
            else
                error("Pick a mode or leave argument blank and it will use Gaussian quadrature")
            end
            # for the entire integral value; only add up what was accumulated in this interval
            stTemp[st,int] = stInt[int][st,end];

            # for the actual points on the graph, we need to offset by the initial value
            if int > 1 # add on intitial value
              stInt[int][st,1:Nck[int]] = stInt[int][st,1:Nck[int]] + stInt[int-1][st,end];
            end
            idx=idx+1;
        end
    end
    stIntVal = sum(stTemp,2);

    # control integral
    ctrInt = [zeros(Float64, numControls, Nck[int],)for int in 1:Ni]; idx = 1;
    ctrTemp = zeros(Float64, numControls, Ni);
    for ctr in 1:numControls
        for int in 1:Ni
            if mode == :default
                ctrInt[int][ctr,1:Nck[int]] =  cumsum(ωₛ[int].*controlMatrix[int][:,ctr]);
            elseif mode == :LGRIM     # Legendre-Gauss-Radau integration matrix (LGRIM)
                ctrInt[int][ctr,1:Nck[int]] =  IMatrix[int]*controlMatrix[int][:,ctr];
            else
                error("Pick a mode or leave argument blank and it will use Gaussian quadrature")
            end
            # for the entire integral value; only add up what was accumulated in this interval
            ctrTemp[ctr,int] = ctrInt[int][ctr,end];

            # for the actual points on the graph, we need to offset by the initial value
            if int > 1 # add on intitial value
              ctrInt[int][ctr,1:Nck[int]] = ctrInt[int][ctr,1:Nck[int]] + ctrInt[int-1][ctr,end];
            end
            idx=idx+1;
        end
    end
    ctrIntVal = sum(ctrTemp,2);
    return  stInt, stIntVal, ctrInt, ctrIntVal
end

#Expr = integrate(mdl,ps,u,idx=1;C=0.5,(:variable=>:control),(:integrand=>:squared))

function integrate(mdl::JuMP.Model,ps::PS_data,V::Array{JuMP.Variable,1}; C::Float64=1.0, kwargs...)
  @unpack ωₛ, Ni, Nck = ps

  kw = Dict(kwargs);
  if !haskey(kw,:mode); kw_ = Dict(:mode => :quadrature); mode = get(kw_,:mode,0);
  else; mode  = get(kw,:mode,0);
  end

  if !haskey(kw,:integrand); kw_ = Dict(:integrand => :default); integrand = get(kw_,:integrand,0);
  else; integrand = get(kw,:integrand,0);
  end

  variable = get(kw,:variable,0);
  if variable == :state; Nck_cum  = [0;cumsum(Nck+1)];
  elseif variable == :control; Nck_cum = [0;cumsum(Nck)];
  else; error("\n Set the variable to either (:variable => :state) or (:variable => :control). \n")
  end

  if mode == :quadrature
    if integrand == :default      # integrate V
      @NLexpression(mdl, temp[int=1:Ni], sum((ωₛ[int])[j] * (V[Nck_cum[int] + 1:Nck_cum[int + 1]])[j] for j = 1:Nck[int]));
      Expr =  @NLexpression(mdl, C*sum(temp[int] for int = 1:Ni));
    elseif integrand == :squared # integrate V^2
      @NLexpression(mdl, temp[int=1:Ni],C*sum((ωₛ[int])[j] * (V[Nck_cum[int] + 1:Nck_cum[int + 1]])[j] * (V[Nck_cum[int] + 1:Nck_cum[int + 1]])[j] for j = 1:Nck[int]));
      Expr =  @NLexpression(mdl, sum(temp[int] for int = 1:Ni));
    else
      error("It should not get here...")
    end
  elseif mode == :LGRIM# TODO add in option to allow for integration using IMatrix
      error("\n Not implemented yet!! \n")
  end
  return Expr
end
"""
dx = differentiate_state(ps,nlp;(:mode=>:something...no modes yet!))
dx = differentiate_state(ps,nlp)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/4/2017, Last Modified: 1/4/2017 \n
--------------------------------------------------------------------------------------\n
"""
function differentiate_state(ps::PS_data,nlp::NLP_data; kwargs...)
    @unpack Nck, Ni, stateMatrix, ωₛ, DMatrix = ps
    @unpack numStates = nlp

    kw = Dict(kwargs);
    if !haskey(kw,:mode); kw = Dict(:mode => :default) end
    mode = get(kw, :mode, 0);

    dx = [zeros(Float64, numStates, Nck[int],)for int in 1:Ni]; idx = 1;
    for st in 1:numStates
        for int in 1:Ni
            if mode == :default
                dx[int][st,1:Nck[int]] = DMatrix[int]*stateMatrix[int][:,st];
            else
                error("Pick a mode or leave argument blank default")
            end
            idx=idx+1;
        end
    end
    return  dx
end

function scale_tau(τ::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*τ + (xₙ + x₀)/2;
end

function scale_w(ω::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*ω;
end

"""
di, tm, ts, ωₛ=create_intervals(t0,tf,Ni,Nck,τ,ω);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/23/2017, Last Modified: 1/2/2017 \n
--------------------------------------------------------------------------\n
"""
function create_intervals(t0::Float64,tf::Float64,Ni::Int64,Nck::Array{Int64,1},τ::Array{Vector{Float64},1},ω::Array{Vector{Float64},1})
  di = (tf - t0)/Ni; # interval size
  # create mesh points
  tm = zeros(Float64,Ni+1); tm[1] = t0;
  for idx in 1:Ni
    tm[idx+1] = tm[idx] + di;
  end
    # go through each mesh interval creating time intervals; [t(i-1),t(i)] --> [-1,1]
    ts=[[scale_tau(τ[int],tm[int],tm[int+1]);di*int] for int in 1:Ni];
    ωₛ=[scale_w(ω[int],tm[int],tm[int+1]) for int in 1:Ni];
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
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
      end
    end
    return L
end
lagrange_basis_poly{T<:Number}(x::AbstractArray{T},x_data,Nc) = lagrange_basis_poly!(x::AbstractArray{T},x_data,Nc,zeros(x))

"""
y=interpolate_lagrange(ts[int],ts[int],stateMatrix[int][:,st],Nck[int])
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/2/2017, Last Modified: 1/9/2017 \n
--------------------------------------------------------------------------------------\n
"""
function interpolate_lagrange{T<:Number}(x::AbstractArray{T},x_data,y_data,Nc)
    if Nc > length(x_data) -1
      error("\n Maximum Nc value = length(x_data)-1 \n")
    end
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

    ns = length(x);
    L = zeros(Float64,Nc+1,ns);
    x = x[:]; x_data = x_data[:]; y_data = y_data[:]; # make sure data is in a column
    for idx in 1:Nc+1
      for j in 1:ns
        L[idx,j] = lagrange_basis_poly(x[j],x_data,Nc,idx);
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
          L=eye(Bool, N, N);
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
