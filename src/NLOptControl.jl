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


export
  # Objects

  # Functions

  # Macros and Support Functions

  # MAKE SURE YOU REMOVE THE FINAL COMMA!!

@with_kw immutable Params @deftype Float64
  # Define Parameters

end
