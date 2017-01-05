"""
LGR_matrices(ps,nlp)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/3/2017, Last Modified: 1/4/2017 \n
--------------------------------------------------------------------------\n
"""
function LGR_matrices(ps::PS_data,nlp::NLP_data)
    @unpack Nck, Ni, ts = ps
    @unpack DMatrix, IMatrix = ps

    # check input
    check_ts = zeros(Float64,Ni);
    for int in 1:Ni
        check_ts[int]=maximum(ts[int]);
    end
    if maximum(check_ts) < 10*eps()
      error("\n ts is full of zeros: make sure that you call create_intervals() first to calculate ts! \n ")
    end

    D = [zeros((Nck[int]+1),(Nck[int]+1)) for int in 1:Ni];
    for int in 1:Ni
        D[int] = poldif(ts[int], 1) # +1 is already appended onto ts
    end

    DM = [zeros((Nck[int]),(Nck[int])) for int in 1:Ni];
    for int in 1:Ni
        DMatrix[int] = D[int][1:end-1,:];   # [Nck]X[Nck+1]
        DM[int] = D[int][1:end-1,2:end];    # [Nck]X[Nck]
        IMatrix[int] = inv(DM[int]);        # I = inv[D_{2:N_k+1}]
    end
    @pack ps = DMatrix, IMatrix
end

"""
xd=lgr_diff(Nc,Ni,x,t_data)
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 12/30/2016, Last Modified: 12/30/2016
--------------------------------------------------------------------------\n
# Input Arguments
* `Nc::Integer`: number of collocation points
* `Ni::Integer`: number of intervals
* `x`: variable data with size = [Nc+1]x[Ni]
* `t_data`: time data with size = [Nc+1]x[Ni]
# Output Arguments
* `xd`: variable derivative

LGR points:

 * roots of: ``P{Nc-1}(τ)+P_{Nc}(τ)``
 * exact for polynomials with: ``degree <= 2Nc-2``
"""
function lgr_diff(Nc::Integer,Ni::Integer,x,t_data)
  if ndims(x)==1
    x=x[:]; # make sure that x is in a column vector
    row = length(x);
    if row != Nc+1
        error("Make sure that `x`: variable data has size of [Nc+1]x[Ni==1]")
    end
  else
    row, column = size(x);
    if row != Nc+1  || column != Ni
        error("Make sure that `x`: variable data has size of [Nc+1]x[Ni]")
    end
  end

  D = zeros(Float64,Nc+1,Nc+1,Ni);
  for int in 1:Ni #TODO see if I can get ride of this, cashe result and only calculate once (if number of nodes are same in each interval)
    D[:,:,int] = poldif(t_data[:,int], 1);
  end
  D=D[1:end-1,:,:];

  # approximate derivative
  xd = zeros(Float64,Nc,Ni);
  for int in 1:Ni
      xd[:,int] = D[:,:,int]*x[:,int]
  end

  return xd
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
