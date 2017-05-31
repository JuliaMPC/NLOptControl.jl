
"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 5/31/2017, Last Modified: 5/31/2017 \n
-------------------------------------------------------------------------------------\n
"""
macro DiffEq(name,eq)
  return quote
    function $(esc(name)){T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
      if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      dx=Array(Any,L,n.numStates);
        for st in 1:n.numStates
          dx[:,st]=@NLexpression(n.mdl,[j=1:L], $(eq)[st] )
        end
      return dx
    end
  end
end
