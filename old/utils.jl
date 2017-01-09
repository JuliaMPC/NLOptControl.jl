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
  # Exact for polynomial of order <= 2N - 3
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

#--------------------------------------------------------------------------
# Author: Huckleberry Febbo, Graduate Student, University of Michigan
# Date: 12/23/2015
# Last modifed on December 23, 2016 by Huckleberry Febbo
#--------------------------------------------------------------------------
# From a modified version of what is described in the paper (*) titled:
# "Legendre–Gauss–Lobatto Pseudo–spectral Method for One–Dimensional Advection–Diffusion Equation"
function MYD(τ::Vector{Float64},N::Int64)
  D = zeros(N+1,N+1)
  for i in 1:N+1
    for j in 1:N+1
      if i!=j
        D[i,j] = lepoly(N,τ[i],false)[1]/(lepoly(N,τ[j],false)[1]*(τ[i]-τ[j]));
      elseif (i==j && j==1)
        D[i,j] = -N*(N+1)/4;
      elseif (i==j && j==(N+1)) #TODO added the N+1 --> check this (what I modified from the paper(*))
        D[i,j] = N*(N+1)/4;
      else
        D[i,j] = 0;
      end
    end
  end
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

# https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Julia
# (This code is a simplified version of the Base.gauss subroutine in the Julia standard library.)
#=
julia> x,w = gauss(-3,3,5)
([-2.718539537815986,-1.6154079303170459,1.3322676295501878e-15,1.6154079303170512,2.7185395378159916],[0.7107806551685667,1.4358860114981036,1.706666666666666,1.435886011498098,0.7107806551685655])
julia> sum(exp(x) .* w)
20.03557771838554
=#
#= page 99, of
@Book{Shen2011,
  title     = {Spectral methods: algorithms, analysis and applications},
  publisher = {Springer Science \& Business Media},
  year      = {2011},
  author    = {Shen, Jie and Tang, Tao and Wang, Li-Lian},
  volume    = {41},
  file      = {:bok#3A978-3-540-71041-7.pdf:PDF},
  groups    = {Spectral Methods},
}
"The eigenvalue method is well-suited for the Gauss-quadratures of low or mod-
erate order. However, for high-order quadratures, the eigenvalue method may suffer
from round-off errors, so it is advisable to use a root-finding iterative approach. To
fix the idea, we restrict our attention to the commonly used Legendre-Gauss-Lobatto
case and compute the zeros of L`N (x). In this case, the Newton method"
=#
function gauss(a,b, N)
    λ, Q = eig(SymTridiagonal(zeros(N), [ n / sqrt(4n^2 - 1) for n = 1:N-1 ]))
    return (λ + 1) * (b-a)/2 + a, [ 2*Q[1,i]^2 for i = 1:N ] * (b-a)/2
end

"""
τ,ω,I,D = LGR(10);
--------------------------------------------------------------------------\n
Last modifed on December 23, 2016 by Huckleberry Febbo\n
Original Author: Jiechao Liu, University of Michigan\n
Original Function Name: LGR_Nodes.m  |  Source: OCOA_150312\n
--------------------------------------------------------------------------\n
# Input Arguments
* `N::Integer`:  number of colocation points
# Output Arguments
* `D::Array{Float64,2}`: Radau Psueudospectral Differention Matrix
     * N X (N + 1) non-square matrix
     * Has one more column than row because the of the Lagrange Polynomial associated with the non-collocated point at τ₀ = -1
     * non-singular
     * integration is exact for polynomials of degree 2N - 2 [more info here](http://users.clas.ufl.edu/hager/papers/Control/unified.pdf)

LGR points:

 * roots of: ``P{N-1}(τ)+P_{N}(τ)``
 * exact for polynomials with: ``degree <= 2N-2``
"""
function LGR(N::Int64)
  # The Legendre Vandermonde Matrix
  P	= zeros(N,N+1);
  # row i: x(i)
  # col j: P_{j-1}(x(i))

  # initial guess
  xn	= - cos(2*pi*(0:(N-1))/(2*(N-1)+1))'; # new x

  # any number larger than 1 to initialize the while() loop
  xo	= 2; # old x

  # Newton-Raphson method
  while maximum(abs(xn - xo)) > eps()*10
      xo = xn;
      Pvec = [0:N];       # initialize the P matrix
      for idx = 1:N+1
        P[1,idx]    	= (-1)^Pvec[1][idx];
      end
      P[2:N,1] 	= 1;
      P[2:N,2]  	= xn[2:N];
      # use Bonnet�s recursion formula to complete the P matrix
      for i = 2:N
          P[2:N,i+1] = ((2*i-1)*xn[2:N].*P[2:N, i] - (i-1)*P[2:N, i-1])/ i;
      end
      FCN         =    P[2:N, N+1] + P[2:N, N];
      DER         = N*(P[2:N, N+1] - P[2:N, N])./(xo[2:N] - 1);
      xn[2:N]  	= xo[2:N] - FCN ./ DER;
  end
  TAU = xn;

  # The LGR Weights
  WEIGHT = Array(Float64,1, N)
  WEIGHT[1]    	= 2/N^2;
  #WEIGHT[2:N]     = (1-TAU[2:N])./(N*P[2:N,N]).^2; #TODO seems like an error here --> report it
  WEIGHT[2:N] = 1./((1-TAU[2:N]).*(lepoly(N,TAU[2:N],true)[1]).^2);
  # The Barycentric weights used to calculate the differentiation matrix
  temp = append!(TAU[1:N],1.0);
  M               = length(temp);
  Y               = repmat(temp,1,M);
  YDIFF           = Y - Y' + eye(M);
  BW              = repmat(1./prod(YDIFF,1),M,1);
  # TODO check above with this http://jiao.ams.sunysb.edu/teaching/ams527_spring15/lectures/BarycentricLagrange.pdf



  # The LGR differentiation matrix N x (N + 1)
  DMAT           	= BW./(BW'.*YDIFF');
  DMAT[1:M+1:M*M]	= sum(1./YDIFF, 1) - 1;
  DMAT           	= DMAT[1:N,:];

  # The LGR integration matrix
  IMAT         	= inv(DMAT[:,2:N+1]);

  τ = [TAU[1:N]; 1];  # append +1 on the end
  ω = WEIGHT[1:N];
  I = IMAT;
  D = DMAT;
  return τ,ω,I,D
end
#=
https://math.la.asu.edu/~bdw/MAT535/Sp06/

function D = diffmat(x)
x = x(:);
[x,p] = sort(x);
m = length(x);
D = zeros(m,m);
for j = 1:m
    i = [1:j-1 j+1:m]';
    D(i,j) = 1./(x(i)-x(j));
end
S = zeros(m,m);
a = sum(log(abs(x([2:m])-x(1)))); S(1,1) = 1;
for i = 2:m
    S(i,i) = exp(-(sum(log(abs(x([1:i-1,i+1:m])-x(i))))-a));
end
for i = 2:2:m
    S(i,i) = -S(i,i);
end
D = S\(D*S);
for i = 1:m
    jp = find(D(i,:)>0);
    [z,pp] = sort(D(i,jp));
    jm = find(D(i,:)<0);
    [z,pm] = sort(D(i,jm));
    D(i,i) = -sum(D(i,jm(pm)))-sum(D(i,jp(pp)));
end
=#

#=
function [tau,w,D] = lgrPS(s,N);

#--------------------------------------------------------#
# This function computes the points, weights, and        #
# differentiation matrix for use with a multiple-segment #
# Radau pseudospectral method.                           #
# INPUTS                                                 #
#  s:  a vector of length K+1 that contains mesh points  #
#      on the closed interval [-1,+1]. s must be         #
#      monotonically increasing on [-1,+1].              #
#  N:  a vector of length K that contains the degree of  #
#      the polynomial within each mesh interval          #
# OUTPUTS                                                #
#  tau: a vector of LGR points on [-1,+1]                #
#  w:   a vector of LGR weights                          #
#  D:   the LGR differentiation matrix                   #
#--------------------------------------------------------#


# https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/RadauStuff/lgrPS.m
# s = [-1; s];
Ntotal = sum(N);
tau = zeros(Ntotal+1,1);
w   = zeros(Ntotal,1);
D   = sparse(Ntotal,Ntotal+1);
rowshift = 0;
colshift = 0;
for i=1:length(N)
  [xcurr,wcurr,Pcurr] = lgrnodes(N(i)-1);
  rows = rowshift+1:rowshift+N(i);
  cols = colshift+1:colshift+N(i)+1;
  scurr = (s(i+1)-s(i))*xcurr/2+(s(i+1)+s(i))/2;
  tau(rows) = scurr;
  w(rows) = (s(i+1)-s(i))*wcurr/2;
  scurr = [scurr; s(i+1)];
  Dcurr = collocD(scurr);
  Dcurr = Dcurr(1:end-1,:);
  D(rows,cols) = Dcurr;
  rowshift = rowshift+N(i);
  colshift = colshift+N(i);
end;
tau(end)=+1;

=#


#=
function [x,w,P]=lgrnodes(N)

################################################################################
#
# lgrnodes.m
#
  # Computes the Legendre-Gauss-Radau nodes, weights and the LGR Vandermonde
  # matrix. The LGR nodes are the zeros of P_N(x)+P_{N+1}(x).
#
  # References on LGR nodes and weights:
#   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
#   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
#
							   #   F. B. Hildebrand , "Introduction to Numerical Analysis," Section 8.11
#   Dover 1987
#
# Written by Greg von Winckel - 05/02/2004
							   # Contact: gregvw@chtm.unm.edu
#
################################################################################

# Truncation + 1
							   N1=N+1;

# Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
x=-cos(2*pi*(0:N)/(2*N+1))';

# The Legendre Vandermonde Matrix
P=zeros(N1,N1+1);

# Compute P_(N) using the recursion relation
# Compute its first and second derivatives and
# update x using the Newton-Raphson method.

xold=2;

# Free abscissae
free=2:N1;

while max(abs(x-xold))>eps

    xold=x;

    P(1,:)=(-1).^(0:N1);

    P(free,1)=1;    P(free,2)=x(free);

    for k=2:N1
        P(free,k+1)=( (2*k-1)*x(free).*P(free,k)-(k-1)*P(free,k-1) )/k;
    end

    x(free)=xold(free)-((1-xold(free))/N1).*(P(free,N1)+P(free,N1+1))...
                ./(P(free,N1)-P(free,N1+1));
        end

# The Legendre-Gauss-Radau Vandermonde
P=P(1:N1,1:N1);

# Compute the weights
w=zeros(N1,1);
w(1)=2/N1^2;
w(free)=(1-x(free))./(N1*P(free,N1)).^2;

=#


function collocD(x);

################################################################################
#
# collocD.m
#
# Computes the pseudospectral/collocation differentiation matrix for the
# arbitrary nodes stored in the vector x. Uses the lagrange polynomial
# formulation.
#
# Reference:
#
# Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange Interpolation"
  # http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
#   https://www.mathworks.com/matlabcentral/fileexchange/5511-2d-barycentric-lagrange-interpolation/content/barylag2d.m
#https://www.mathworks.com/matlabcentral/fileexchange/?term=authorid%3A11897
  # Written by: Greg von Winckel       07/18/04
  # Contact:    gregvw@chtm.unm.edu
#
################################################################################

# Make x a column vector if it isn't already and order it
# and get the number of nodes
x=sort(x[:]);                       N=length(x); N1=N+1; N2=N*N;

# Compute the barycentric weights
X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
W=repmat(1./prod(Xdiff,2),1,N);

D=W./(W'.*Xdiff);
D[1:N1:N2]=1-sum(D);                D=-D';
end

#https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite
function poldif(x, malpha, B...)
  B_given = Bool(length(B));
  #  The function DM =  poldif(x, maplha, B) computes the
  #  differentiation matrices D1, D2, ..., DM on arbitrary nodes.
  #
  #  The function is called with either two or three input arguments.
  #  If two input arguments are supplied, the weight function is assumed
  #  to be constant.   If three arguments are supplied, the weights should
  #  be defined as the second and third arguments.
  #
  #  Input (constant weight):
  #
  #  x:        Vector of N distinct nodes.
  #  malpha:   M, the number of derivatives required (integer).
  #  B:        Omitted.
  #
  #  Note:     0 < M < N-1.
  #
  #  Input (non-constant weight):
  #
  #  x:        Vector of N distinct nodes.
  #  malpha:   Vector of weight values alpha(x), evaluated at x = x(k).
  #  B:        Matrix of size M x N,  where M is the highest
  #            derivative required.  It should contain the quantities
  #            B(ell,j) = beta(ell,j) = (ell-th derivative
  #            of alpha(x))/alpha(x),   evaluated at x = x(j).
  #
  #  Output:
  #  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.

  #  J.A.C. Weideman, S.C. Reddy 1998

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
