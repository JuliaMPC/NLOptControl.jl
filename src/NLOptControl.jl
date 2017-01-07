module NLOptControl

using Media, DifferentialEquations, Dierckx, Parameters, Interpolations,FastGaussQuadrature, Polynomials

# To copy a particular piece of code (or function) in some location
macro def(name, definition)
  return quote
    macro $name()
      esc($(Expr(:quote,definition)))
    end
  end
end


"""
ps = PS_data(Nc=Nc,Ni=Ni,τ=τ,ω=ω,t0=t0,tf=tf);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 1/3/2017 \n
--------------------------------------------------------------------------\n
"""
@with_kw type PS_data
  Nck::Array{Int64,1} # number of collocation points per interval
  Ni::Int64 # number of intervals
  #TODO --> for now default is gaussradau, eventually add other PS methods
  τ::Array{Array{Float64,1},1}  # Node points ---> Nc increasing and distinct numbers ∈ [-1,1] #TODO might be able to calculate here during problem initialization
  ts::Array{Array{Float64,1},1} # time scaled based off of τ
  ω::Array{Array{Float64,1},1} # weights
  ωₛ::Array{Array{Float64,1},1} # scaled weights
  t0::Float64 # initial time
  tf::Float64 # final time
  DMatrix::Array{Array{Float64,2},1}  # differention matrix
  IMatrix::Array{Array{Float64,2},1}  # integration matrix
  FMatrix::Array{Array{Float64,2},1}
  stateMatrix::Array{Array{Float64,2},1}
  controlMatrix::Array{Array{Float64,2},1}
end

"""
nlp = NLP_data(numStates=numStates, numControls=numControls);
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 1/1/2017, Last Modified: 1/3/2017 \n
Citations: \n
----------\n
Original Author: S. Hughes.  steven.p.hughes@nasa.gov
Source: DecisionVector.m [located here](https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/LowThrust/@DecisionVector/)
--------------------------------------------------------------------------\n
"""

@with_kw type NLP_data

  # general properties
  numStates::Int64      # number of states
  numControls::Int64    # number of controls

  numStatePoints::Array{Int64,1} # number of state discretization within each interval
  numControlPoints::Array{Int64,1} # number of control discretization within each interval

  lengthStateVector::Int64  # total number of state variables
  lengthControlVector::Int64 # total number of control parameters

  # properties for entire decision vector
  lengthDecVector::Int64

  # vector indeces
  timeStartIdx::Int64
  timeStopIdx::Int64
  stateIdx::Array{Tuple{Int64,Int64},1}
  controlIdx::Array{Tuple{Int64,Int64},1}
  stateIdx_all::Array{Tuple{Int64,Int64},1}
  controlIdx_all::Array{Tuple{Int64,Int64},1}
  stateIdx_st::Array{Tuple{Int64,Int64},1}
  controlIdx_ctr::Array{Tuple{Int64,Int64},1}
  # problem data
  decisionVector::Vector{Float64}

  # state equations
  stateEquations::Array{Function,1}
end

# scripts -> need to be called after the data!
include("utils.jl");
include("test_data.jl")
include("LGL.jl");
include("LGR.jl");

  # Objects

  # Functions
export initialize_NLP, nlp2ocp,
       LGL_nodes, LGL_weights, LGL_Dmatrix,
       LGR, lgr_diff, LGR_matrices, F_matrix,
       scale_tau, scale_w, create_intervals,
       lagrange_basis_poly, interpolate_lagrange,
       integrate_state,
       differentiate_state,
       lepoly, poldif,

       generate_Fake_data,

       # objects
       NLP_data, PS_data,

      # Macros and Support Functions
      @unpack_NLP_data,
      @pack_NLP_data,
      @unpack_PS_data,
      @pack_PS_data

  # MAKE SURE YOU REMOVE THE FINAL COMMA!!

# ---------NOTES--------------
# ns: number of state variables
# nc: number of control variables
# τ: interval that collocation points are defined on:
  # LG: τ ∈ [-1, 1)
  # LGR: τ ∈ [-1, 1)
  # backward LGR: τ ∈ (-1, 1]
  # LGL: τ ∈ [-1, 1]
# N: number of collocation points
# P_{N}(τ): Nth degree Legendre polynomial
  # LG points:
    # are the roots of: P_{N}(τ)
    # exact for polynomials with degree <= 2N-1
  # LGR points:
    # are the roots of: P{N-1}(τ)+P_{N}(τ)
    # exact for polynomials with degree <= 2N-2
  # backward LGR points:
    # are the negative roots of: P{N-1}(τ)+P_{N}(τ)
    # exact for polynomials with degree <= 2N-2
  # LGL points:
    # are the roots of: \dot{P}{N-1}(τ) with -1 and _1
    # exact for polynomials with degree <= 2N-3

#= states:
  # ξ₀: intial state vectors

  # all states at τᵢ
    # ξᵢ(τ) = [ξ₁(τᵢ),....,ξns(τᵢ)], 1<=i<=N
    # vectors are row vectors

  # full state matrix--> size (N+1) x n ; includes ξ₀
      Ξ(τ) = [ξ₀;
              ξᵢ;
              .
              .
              .
              ξN]

   # partial state matrix--> size (N) x n
      Ξp(τ) = ξᵢ;
              .
              .
              .
              ξN]

=#
#-----------------------------

#@with_kw immutable CollocationPoint @deftype Float64

  # Define Parameters

end
