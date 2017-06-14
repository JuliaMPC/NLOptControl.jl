"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 5/31/2017, Last Modified: 6/14/2017 \n
-------------------------------------------------------------------------------------\n
"""
macro DiffEq(n,expr)
  n=esc(n)
  return quote
      $n.stateEquations=$(Expr(:quote,expr))
    end
end


function dx_EQ{T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
  dx_expr=parse_DiffEq(n,x,u)
  for st in 1:n.numStates
    for j in 1:n.numStatePoints
      dx[j,st]=create_DiffEq(n,dx_expr[j,st])
    end
  end
  return dx
end

"""
#TODO consider appending a zero so it is not just a jump variable
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/11/2017, Last Modified: 6/11/2017 \n
-------------------------------------------------------------------------------------\n
"""

function parse_DiffEq{T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
  if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end

  dx=Array(Expr,L,n.numStates);
  for qq in 1:L
    dx[qq,:]=deepcopy(n.DXexpr)
  end

  #TODO use the statenames here! -> currently running too often!
  state=Array(Symbol,n.numStates,1)
  for st in 1:n.numStates
    state[st]=Symbol(string("x",:($(st))))
  end
  control=Array(Symbol,n.numControls,1)
  for ctr in 1:n.numControls
    control[ctr]=Symbol(string("u",:($(ctr))))
  end

  for st_high in 1:n.numStates
    for qq in 1:L #TODO this will not start at 1 for the :ps method...unless you go bac to x and u
      # states
      for st in 1:n.numStates
        dx[qq,st_high]=rePlace(dx[qq,st_high],state[st],:($(n.r.x[qq,st])))
      end

      # controls
      for ctr in 1:n.numControls
        dx[qq,st_high]=rePlace(dx[qq,st_high],control[ctr],:($(n.r.u[qq,ctr])))
      end
    end
  end
  return dx
end
#dx=parse_DiffEq(n,n.r.x,n.r.u)
