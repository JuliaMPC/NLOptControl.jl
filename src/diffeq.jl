"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified: 6/30/2017 \n
-------------------------------------------------------------------------------------\n
"""
function dynamics!(n::NLOpt,dx::Array{Expr,1})
  if length(dx)!=n.numStates
    error("\n The number of differential equations must equal numStates. \n")
  end
  n.DXexpr=dx;
  return nothing
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/11/2017, Last Modified2: 6/30/2017 \n
-------------------------------------------------------------------------------------\n
"""
function DiffEq(n::NLOpt,x::Array{JuMP.Variable,2},u::Array{JuMP.Variable,2},L::Int64,st::Int64)
  expr=n.DXexpr[st]
  return NLExpr(n,expr,x,u,L)
end

"""
# returns an array of @NLexpression()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified2: 7/1/2017 \n
-------------------------------------------------------------------------------------\n
"""
function NLExpr(n::NLOpt,expr::Expr,args...)
  if length(args)==3
    xxx=args[1]
    uuu=args[2]
    L=args[3]
  elseif length(args)==0
    xxx=n.r.x
    uuu=n.r.u
    L=n.numStatePoints
  else
    error("\n The length of the args... must be either 3 or 0 \n")
  end

  code=quote
    # rename state variables
    state=Array{Expr}($n.numStates);
    for st in 1:$n.numStates
      state[st]=Expr(:(=),$n.state.name[st],$xxx[:,st])
      eval(state[st])
    end

    # rename control variables
    control=Array{Expr}($n.numControls);
    for ctr in 1:$n.numControls
      control[ctr]=Expr(:(=),$n.control.name[ctr],$uuu[:,ctr])
      eval(control[ctr])
    end
    
    @NLexpression($n.mdl,[j=1:$L],$expr)
  end
  return eval(code)
end
