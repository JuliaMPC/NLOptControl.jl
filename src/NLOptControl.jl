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

immutable NodeData #TODO --> use this
  Nₜ :: Int  # where Nₜ + 1 is the total number of time steps
  τ :: Vector{Float64} # Node points ---> Nₜ increasing and distinct numbers ∈ [-1,1]
  # these might have to be design variables
end

  # Objects

  # Functions
export LGL_nodes,LGL_weights,LGL_Dmatrix,scale_X

  # Macros and Support Functions

  # MAKE SURE YOU REMOVE THE FINAL COMMA!!

#@with_kw immutable CollocationPoint @deftype Float64

  # Define Parameters

end
