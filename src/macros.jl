"""

--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 5/31/2017, Last Modified: 6/11/2017 \n
-------------------------------------------------------------------------------------\n
"""
macro DiffEq(name,expr)
  return quote
   macro $name()
      $(esc(Expr(:quote,expr)))
    end
  end
end


 ########### NOTE this is an attemt to simplfy the macro
#=
macro DiffEq(name,definition)
  return quote
    #function $(esc(name)){T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
    function $(esc(name))(n::NLOpt,x,u)
      if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      dx=Array(Any,L,n.numStates);
      @show x
      esc($(Expr(:quote,definition)))
    #  @show u
    #  return dx
    end
  end
end

macro DiffEq(name,eq)
  return quote
    function $(esc(name)){T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
      if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      dx=Array(Any,L,n.numStates);
      for st in 1:n.numStates
        dx[:,st]=JuMP.@NLexpression(n.mdl,[j=1:L], (eq)[st] )
      end
      return dx
    end
  end
end
=#
################

########
#= NOTE this is a contained example
using JuMP
type R
  eq
  mdl::JuMP.Model
end
function R()
  R(Any,JuMP.Model())
end

macro defODE(r,eq)
  r=esc(r);
  quote
      r.eq=$(Meta.quot(eq))
      function $(esc(:DE))(r,x)
         @NLexpression($r.mdl,$r.eq)
      end
  end
end

r=R();
@defODE(r,2+sin(x));
@variable(r.mdl,x);
expr=DE(r,x)
=#
#########

###############################
#= NOTE this is a contained example #2

using JuMP
type NLOpt
  eq
  mdl::JuMP.Model
  numStates
  numStatePoints
end
function NLOpt()
  NLOpt(Any,JuMP.Model(),1,4)
end

function create_DiffEq(n::NLOpt,expr::Expr)
  JuMP.NonlinearExpression(n.mdl,JuMP.NonlinearExprData(n.mdl,expr))
end

function create_DiffEq_array(n::NLOpt,dx_expr::Array{Expr,2})
  code=quote
    function dx(n::NLOpt,x)
      #if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      L=size(x)[1]; # temp
      dx=Array(Any,L,n.numStates);
    #  X1=x[:,1];  #TODO think about expanding this
      for st in 1:n.numStates
        for j in 1:n.numStatePoints
          dx[j,st]=create_DiffEq(n,dx_expr[j,st])
        end
      end
      return dx
    end
  end

return eval(code)
end

# Example
n=NLOpt();
@variable(n.mdl,x[1:5,1:2]);
JuMP.initNLP(n.mdl)

#create_DiffEq_array(n,[:(:x1*5)])

# working on parsing
diffEqs=[:($(:x1)*5);:($(:x2)*$(:x1))]

numStates=2;L=5;
dx=Array(Any,L,numStates);
state=Array(Symbol,numStates,1)
for st in 1:numStates
  state[st]=Symbol(string("x",:($(st))))
end

for st_high in 1:numStates
  for st in 1:numStates
    idx=find(diffEqs[st_high].args.==state[st])
    if isempty(idx)
      println(string("empty"))
    else
      println(string("the idx=",idx[1]))
      diffEqs[st_high].args[idx[1]]=:($x)
    end
  end
end
=#
#TODO replace :X1 $x[1,1] etc.

############################ NOTE the pervious one was the latest it started to get very complicated, will continue later...
#=
using JuMP
type NLOpt
  eq
  mdl::JuMP.Model
  numStates
  numStatePoints
end
function NLOpt()
  NLOpt(Any,JuMP.Model(),1,4)
end

function create_DiffEq(n::NLOpt,expr::Expr)
  JuMP.NonlinearExpression(n.mdl,JuMP.NonlinearExprData(n.mdl,expr))
end

function create_DiffEq_array(n::NLOpt,dx_expr::Array{Expr,2})
  code=quote
    function dx(n::NLOpt,x)
      #if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      L=size(x)[1]; # temp
      dx=Array(Any,L,n.numStates);
    #  X1=x[:,1];  #TODO think about expanding this
      for st in 1:n.numStates
        for j in 1:n.numStatePoints
          dx[j,st]=create_DiffEq(n,dx_expr[j,st])
        end
      end
      return dx
    end
  end

return eval(code)
end

# Example
n=NLOpt();
@variable(n.mdl,x[1:4,1:1]);
JuMP.initNLP(n.mdl)

# Variable objects should be spliced into the expression.
dx_expr=Array(Expr,n.numStatePoints,n.numStates)
for j in 1:n.numStatePoints #TODO needs to only go to L
  dx_expr[j,1]=:($(x[j,1])*5);
end

create_DiffEq_array(n,dx_expr)
=#
#########################NOTE the above is the latest working example

#=
# scalar ex
macro test(name,ex)
  return quote
    macro $name(a)
      c=length(a)
      $ex
    end
  end
end
@test t1 begin
  a+c
end
@t1(2)

# symbol ex
macro test(name,ex)
  return quote
    macro $name()
      $(esc(Expr(:quote,ex)))
    end
  end
end

a=
c=length(a.args);
@test t1 begin
  a+c
end

@t1([:x1,:x2])



############################
using JuMP
type NLOpt
  eq
  mdl::JuMP.Model
  numStates
  numStatePoints
end
function NLOpt()
  NLOpt(Any,JuMP.Model(),1,4)
end

function create_DiffEq(n::NLOpt,expr::Expr)
  JuMP.NonlinearExpression(n.mdl,JuMP.NonlinearExprData(n.mdl,expr))
end

n=NLOpt();
@variable(n.mdl,x[1:4,1:1]);
JuMP.initNLP(n.mdl)

# Variable objects should be spliced into the expression.
function create_DiffEq_array(n::NLOpt,expr::Array{Expr,2})

  dx=Array(Any,n.numStatePoints,n.numStates)
  for st in 1:n.numStates
    for j in 1:n.numStatePoints
      dx[j,st]=create_DiffEq(n,expr[j,st])
    end
  end
  return dx

end

# Example
n=NLOpt();
@variable(n.mdl,x[1:4,1:1]);
JuMP.initNLP(n.mdl)
=#
