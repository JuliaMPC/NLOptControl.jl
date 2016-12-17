
# the following code uses the function lepoly() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
function lepoly(n::Int64,x::Vector{Float64},derivative_option::Bool)

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


# the following code uses the function LGL_nodes() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
#--------------------------------------------------------------------------
# LGL_nodes.m
# determines Lagrange-Gauss-Lobatto (LGL) nodes
#--------------------------------------------------------------------------
# tau = LGL_nodes(N)
#   N: number of nodes minus 1, should be an integer greater than 0
# tau: LGL nodes
#--------------------------------------------------------------------------
# Examples:
# tau = LGL_nodes(1)
# -1     1
# tau = LGL_nodes(2)
# -1     0     1
# tau = LGL_nodes(3)
# -1  -0.44721   0.44721   1
#--------------------------------------------------------------------------
# Author: Daniel R. Herber, Graduate Student, University of Illinois at
# Urbana-Champaign
# Date: 06/04/2015
#--------------------------------------------------------------------------
function LGL_nodes(N::Int)
    # See Page 99 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
    # Algorithms, Analysis and Applications, Springer Series in Compuational
    # Mathematics, 41, Springer, 2011.
    # Uses the function: lepoly()
    # Original function: [varargout] = legslb(n) located at
    # http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslb.m

    # Compute the initial guess of the interior LGL points
    thetak = [(4*(idx)-1)*pi/(4*N+2) for idx in 1:N];
    sigmak = [-(1-(N-1)/(8*N^3)-(39-28/sin(thetak[idx])^2)/(384*N^4))*cos(thetak[idx]) for idx in 1:N];
    ze = (sigmak[1:N-1]+sigmak[2:N])/2;
    ep = eps()*10;                           # error tolerance for stopping iteration
    ze1 = ze+ep+1;

    if size(ze)[1]!=0
      while maximum(abs(ze1-ze))>=ep            # Newton's iteration procedure
        ze1 = ze;
        dy,y = lepoly(N,ze,true);
        ze = ze-(1-ze.*ze).*dy./(2*ze.*dy-N*(N+1)*y);  # see Page 99 of the book
      end                                   # around 6 iterations are required for n=100
    end

    left=Float64[-1]; right=Float64[1];
    tau=append!(left,ze); tau=append!(tau,right);

    return tau
end

# the following code uses the function LGL_weights() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
#--------------------------------------------------------------------------
# LGL_weights.m
# determines Gaussian quadrature weights using Lagrange-Gauss-Lobatto (LGL)
# nodes
#--------------------------------------------------------------------------
# w = LGL_weights(tau)
# tau: LGL nodes
#   w: Gaussian quadrature weights
#--------------------------------------------------------------------------
# Author: Daniel R. Herber, Graduate Student, University of Illinois at
# Urbana-Champaign
# Date: 06/04/2015
#--------------------------------------------------------------------------
function LGL_weights(tau::Vector{Float64})
    # number of nodes
    N = length(tau)-1;

    # turn tau into a Vector
    tau_V = Vector{Float64}(1);
    tau_V=tau;
    # See Page 99 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
    # Algorithms, Analysis and Applications, Springer Series in Compuational
    # Mathematics, 41, Springer, 2011.
    # Uses the function: lepoly()
    # Original function: [varargout] = legslb(n) located at
    # http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslb.m
    y = lepoly(N,tau_V[2:end-1],false);
    # Use the weight expression (3.188) to compute the weights
    w = [2/(N*(N+1));2./(N*(N+1)*y.^2);2/(N*(N+1))];

    return w
end

# the following code uses the function LGL_Dmatrix() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
#--------------------------------------------------------------------------
# LGL_Dmatrix.m
# determines approximate differentiation matrix for Legendre-based method
# with LGL nodes
#--------------------------------------------------------------------------
# D = LGL_Dmatrix(tau)
# tau: LGL nodes
#   D: differentiation matrix
#--------------------------------------------------------------------------
# Author: Daniel R. Herber, Graduate Student, University of Illinois at
# Urbana-Champaign
# Date: 06/04/2015
#--------------------------------------------------------------------------
function LGL_Dmatrix(tau::Vector{Float64})
    # number of nodes
    N = length(tau)-1;

    # See Page 110 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
    # Algorithms, Analysis and Applications, Springer Series in Compuational
    # Mathematics, 41, Springer, 2011.
    # Uses the function: lepoly()
    # Original function: D = legslbdiff(n,x) located at
    # http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslbdiff.m
    n = N + 1;
    if n==0 D = []; return D end   # null differentiation matrix
    xx = tau; y = lepoly(n-1,xx,false);
    D = (xx./y)*y'-(1./y)*(xx.*y)'; # compute L_{n-1}(x_j) (x_k-x_j)/L_{n-1}(x_k);
                                    # 1/d_{kj} for k not= j (see (3.203))
    D = D + eye(n);                 # add the identity matrix so that 1./D can be operated
    D = 1./D;
    D = D - eye(n);
    D[1,1] = -n*(n-1)/4;
    D[n,n] = -D[1,1];  # update the diagonal entries

    return D
end

# the following code uses the function LGL_Dmatrix() located at:
# https://www.mathworks.com/matlabcentral/fileexchange/51104-basic-implementation-of-multiple-interval-pseudospectral-methods-to-solve-optimal-control-problems
# Last modifed on December 15, 2016 by Huckleberry Febbo
#--------------------------------------------------------------------------
# Fmatrix.m
# determines differential matrix for multiple-interval psuedospectral
# method
#--------------------------------------------------------------------------
# F = Fmatrix(tau,X,U,t0,tf,p)
# tau: nodes
#   X: states
#   U: control
#  t0: initial time
#  tf: final time
#   p: parameter structure
#   F: differential matrix
#--------------------------------------------------------------------------
# Author: Daniel R. Herber, Graduate Student, University of Illinois at
# Urbana-Champaign
# Date: 06/04/2015
#--------------------------------------------------------------------------
#=
function Fmatrix(tau,X,U,t0,tf,p)
    # number of nodes
    N = length(tau)-1;
    # initialize differential matrix
    f = zeros(N+1,size(X,2));
    # determine differential matrix
    for i = 1:N+1
        f[i,:] = eval([p.deriv,'(tau(i),X(i,:),U(i,:),t0,tf,p)']);
    end
    F = (tf-t0)/2*f; # scale
end
=#

function scale_tau(τ::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*τ + (xₙ + x₀)/2
end

function scale_w(ω::Array{Float64,1},x₀::Float64,xₙ::Float64)
  (xₙ - x₀)/2*ω
end

#=
if mode == :quad_scale       # scale to range of; [-1,1]
  (x - x₀)./(xₙ - x₀)
elseif mode == :prob_scale # scale to range of problem; [xₙ,x₀]
   x₀ + x./(xₙ - x₀)
end
=#
