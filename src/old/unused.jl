"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/28/2017, Last Modified: 6/28/2017 \n
-------------------------------------------------------------------------------------\n
"""
macro controlVariable(n,name,num)
  n=esc(n)
  return quote
      $name=$n.r.ocp.u[:,num]
    end
end




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
  for st in 1:n.ocp.state.num
    for j in 1:n.ocp.state.pts
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
  if n.s.ocp.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end

  dx=Array{Expr}(L,n.ocp.state.num);
  for qq in 1:L
    dx[qq,:]=deepcopy(n.ocp.DXexpr)
  end

  #TODO use the statenames here! -> currently running too often!
  state=Array{Symbol}(n.ocp.state.num,1)
  for st in 1:n.ocp.state.num
    state[st]=Symbol(string("x",:($(st))))
  end
  control=Array{Symbol}(n.ocp.control.num,1)
  for ctr in 1:n.ocp.control.num
    control[ctr]=Symbol(string("u",:($(ctr))))
  end

  for st_high in 1:n.ocp.state.num
    for qq in 1:L #TODO this will not start at 1 for the :ps method...unless you go bac to x and u
      # states
      for st in 1:n.ocp.state.num
        dx[qq,st_high]=rePlace(dx[qq,st_high],state[st],:($(n.r.ocp.x[qq,st])))
      end

      # controls
      for ctr in 1:n.ocp.control.num
        dx[qq,st_high]=rePlace(dx[qq,st_high],control[ctr],:($(n.r.ocp.u[qq,ctr])))
      end
    end
  end
  return dx
end
#dx=parse_DiffEq(n,n.r.ocp.x,n.r.ocp.u)

"""
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 6/12/2017, Last Modified: 6/12/2017 \n
Citations: \n
----------\n
Based off roughly of:
Source: [located here](https://github.com/JuliaLang/julia/blob/master/base/cartesian.jl#L278-L353)
-------------------------------------------------------------------------------------\n
"""
type Replace
    sym::Symbol
    val::JuMP.Variable
end

function rePlaceEX(ex::Expr,r::Replace)
  for i in 1:length(ex.args)
    if ex.args[i]==r.sym
      ex.args[i]=r.val
    end
  end
  return ex
end

function dig(ex)
  branch=falses(length(ex.args));
  for i in 1:length(ex.args)
    if typeof(ex.args[i])==Expr;
      branch[i]=true;
    end
  end
  return branch
end

function rePlace(ex::Expr,sym::Symbol,val::JuMP.Variable)
  r=Replace(sym,val)
  ex_tmp=[deepcopy(ex)]
  unwrap=true
  level=1
  while unwrap
    expr_arr=[];startCT=true;
    for j in 1:length(ex_tmp)
      ex_tmp[j]=rePlaceEX(ex_tmp[j],r);   # replace any on this branch
      branch=dig(ex_tmp[j]);
      if any(branch) && startCT
        expr_arr=ex_tmp[j].args[branch];
        tree=[any(branch)];
        startCT=false;
      elseif any(branch)
        expr_arr=[expr_arr,ex_tmp[j].args[branch]];
        tree=[tree,any(branch)];
      end
    end

    if level==1; ex=ex_tmp; end

    if isempty(expr_arr)
      unwrap=false
    else
      ex_tmp=expr_arr
      level+=1;
      if level>100
        error("\n Are the variables in the expression nested > 100 levels deep? If not report a bug, if so, consider increasing this limit \n ")
      end
    end
  end
  return ex[1]
end

#= TODO make test function
using JuMP
m=Model()
@variable(m,x[1:2,1:4])
expr=:($:x3*3*sin($:x3)*cos(cos($:x3)))
expr=rePlace(expr,:x3,:($(x[2,3])))
=#
