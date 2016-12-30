module NLOptControl


using Media, DifferentialEquations, Dierckx, Parameters, Interpolations

# To copy a particular piece of code (or function) in some location
macro def(name, definition)
  return quote
    macro $name()
      esc($(Expr(:quote,definition)))
    end
  end
end

include("utils.jl")
include("LGL.jl")
include("LGR.jl")

immutable NodeData #TODO --> use this
  #Nₜ :: Int  # where Nₜ + 1 is the total number of time steps
  #τ :: Vector{Float64} # Node points ---> Nₜ increasing and distinct numbers ∈ [-1,1]
  # these might have to be design variables
end

  # Objects

  # Functions
export LGL_nodes,LGL_weights,LGL_Dmatrix,
       LGR, lgr_diff,
       scale_tau,scale_w,create_intervals,
       lagrange_basis_poly,interpolate_lagrange,
       lepoly, poldif

  # Macros and Support Functions

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
