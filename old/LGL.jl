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
