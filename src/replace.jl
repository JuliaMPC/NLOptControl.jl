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
