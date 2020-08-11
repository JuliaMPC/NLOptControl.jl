"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified: 4/13/2018 \n
-------------------------------------------------------------------------------------\n
"""
function dynamics!(n::NLOpt, dx::Vector{Expr})
    if length(dx)!=n.ocp.state.num
        error("\n The number of differential equations must equal ocp.state.num. \n")
    end
    n.ocp.DXexpr = dx
    return nothing
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 7/04/2017, Last Modified: 4/13/2018 \n
-------------------------------------------------------------------------------------\n
"""
function constraints!(n::NLOpt, con::Vector{Expr})
    n.ocp.NLcon = con
    return nothing
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/11/2017, Last Modified: 4/13/2018 \n
-------------------------------------------------------------------------------------\n
"""
#function DiffEq(n::NLOpt,x::Matrix{JuMP.JuMPTypes},u::Matrix{JuMP.JuMPTypes},L::Int,st::Int)::Vector{JuMP.NonlinearExpression}
#function DiffEq(n::NLOpt,x,u,L::Int64,st::Int64)
function DiffEq(n::NLOpt,x::Array{Any,2},u::Array{Any,2},L::Int64,st::Int64)
    expr = n.ocp.DXexpr[st]
    return NLExpr(n,expr,x,u,L)
end

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 7/04/2017, Last Modified: 4/13/2018 \n
-------------------------------------------------------------------------------------\n
"""
#function addCon(n::NLOpt,x::Matrix{JuMP.JuMPTypes},u::Matrix{JuMP.JuMPTypes},L::Int,num::Int)
function addCon(n::NLOpt,x::Array{Any,2},u::Array{Any,2},L::Int64,num::Int64)
    expr = n.ocp.NLcon[num]
    return NLCon(n,expr,x,u,L)
end

"""
# returns an array of @NLexpression()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/29/2017, Last Modified: 12/06/2019 \n
-------------------------------------------------------------------------------------\n
"""
function NLExpr(n::NLOpt, expr::Expr, args...)::Vector{JuMP.NonlinearExpression}
#function NLExpr(n::NLOpt,expr::Expr,args...)

    if length(args) == 3
        xxx = args[1]
        uuu = args[2]
        L = args[3]
    elseif length(args) == 0
        xxx = n.r.ocp.x
        uuu = n.r.ocp.u
        L = n.ocp.state.pts
    else
        error("The length of the args... ($(length(args))) must be either 3 or 0")
    end

    code = quote
        # rename state variables
        state = Array{Expr}(undef, $n.ocp.state.num)
        for st in 1:$n.ocp.state.num
            state[st] = Expr(:(=), $n.ocp.state.name[st], $xxx[:,st])
            eval(state[st])
        end

        # rename control variables
        control = Array{Expr}(undef, $n.ocp.control.num)
        for ctr in 1:$n.ocp.control.num
            control[ctr] = Expr(:(=), $n.ocp.control.name[ctr], $uuu[:,ctr])
            eval(control[ctr])
        end

        @NLexpression($n.ocp.mdl, [j=1:$L], $expr)
    end

    return NLOptControl.eval(code)

end

"""
# returns an array of @NLconstraint()
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 7/04/2017, Last Modified: 12/06/2019 \n
-------------------------------------------------------------------------------------\n
"""
function NLCon(n::NLOpt,expr::Expr,args...)::Vector{JuMP.NonlinearConstraint}
  if length(args) == 3
    xxx = args[1]
    uuu = args[2]
    L = args[3]
  elseif length(args) == 0
    xxx = n.r.ocp.x
    uuu = n.r.ocp.u
    L = n.ocp.state.pts
  else
    error("\n The length of the args... must be either 3 or 0 \n")
  end

  code = quote
    # rename state variables
    state = Array{Expr}(undef, $n.ocp.state.num);
    for st in 1:$n.ocp.state.num
      state[st] = Expr(:(=),$n.ocp.state.name[st],$xxx[:,st])
      eval(state[st])
    end

    # rename control variables
    control = Array{Expr}(undef, $n.ocp.control.num);
    for ctr in 1:$n.ocp.control.num
      control[ctr] = Expr(:(=),$n.ocp.control.name[ctr],$uuu[:,ctr])
      eval(control[ctr])
    end

    @NLconstraint($n.ocp.mdl,[j=1:$L],$expr)
  end
  return eval(code)
end
