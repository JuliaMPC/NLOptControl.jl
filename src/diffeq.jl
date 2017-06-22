
"""
#TODO add other expressions
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/11/2017, Last Modified: 6/13/2017 \n
-------------------------------------------------------------------------------------\n
"""
function DiffEq(n::NLOpt,x,u,L::Int64,st::Int64)

  DXexpr=n.DXexpr[st]
  code=quote
    # rename state variables
    state=Array{Expr}($L,$n.numStates);
    for st in 1:$n.numStates
      state[st]=Expr(:(=),Symbol("x",st),$x[:,st])
      eval(state[st])
    end

    # rename control variables
    control=Array{Expr}($L,$n.numControls);
    for ctr in 1:$n.numControls
      control[ctr]=Expr(:(=),Symbol("u",ctr),$u[:,ctr])
      eval(control[ctr])
    end

    @NLexpression($n.mdl,[j=1:$L],$DXexpr)
  end
  return eval(code)
end
