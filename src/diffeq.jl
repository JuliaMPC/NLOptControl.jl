"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified: 6/29/2017 \n
-------------------------------------------------------------------------------------\n
"""
function dynamics!(n::NLOpt,dx::Array{Expr,1})
  if length(dx)!=n.numStates
    error("\n the number of differential equations must equal numStates \n")
  end
  n.DXexpr=dx;
  n.stateEquations=DiffEq; 
  return nothing
end
"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/11/2017, Last Modified2: 6/29/2017 \n
-------------------------------------------------------------------------------------\n
"""
function DiffEq(n::NLOpt,x::Array{JuMP.Variable,2},u::Array{JuMP.Variable,2},L::Int64,st::Int64)
  expr=n.DXexpr[st]
  return NLExpr(n,expr,x,u,L)
end
"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified2: 6/29/2017 \n
-------------------------------------------------------------------------------------\n
"""
function NLExpr(n::NLOpt,expr::Expr,x::Array{JuMP.Variable,2},u::Array{JuMP.Variable,2},L::Int64)
  code=quote
    # rename state variables
    state=Array{Expr}($L,$n.numStates);
    for st in 1:$n.numStates
      state[st]=Expr(:(=),$n.state.name[st],$x[:,st])
      eval(state[st])
    end

    # rename control variables
    control=Array{Expr}($L,$n.numControls);
    for ctr in 1:$n.numControls
      control[ctr]=Expr(:(=),$n.control.name[ctr],$u[:,ctr])
      eval(control[ctr])
    end

    @NLexpression($n.mdl,[j=1:$L],$expr)
  end
  return eval(code)
end
