
"""
#TODO pass n
--------------------------------------------------------------------------------------\n
Author: Huckleberry Febbo, Graduate Student, University of Michigan
Date Create: 5/31/2017, Last Modified: 5/31/2017 \n
-------------------------------------------------------------------------------------\n
"""
macro DiffEq(name,eq)
#  name = esc(name)
println("test")

  return quote
    function $(esc(name)){T<:Any}(n::NLOpt,x::Array{T,2},u::Array{T,2})
      println("test1")
      @show n
      if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      dx=Array(Any,L,n.numStates);
      for st in 1:n.numStates
      dx[:,st]=JuMP.@NLexpression(n.mdl,[j=1:L], (eq)[st] ) 
      #=
      for j in 1:L
      dx[j,st]=JuMP.NonlinearExpression($n.mdl,JuMP.@processNLExpr($n.mdl, $(esc(eq))[j,st] ))
      end
      =#
      end
      return dx
    end
  end
end

#=
macro DiffEq(name,eq)
  @show eq
  return quote
    function $(esc(name)){T<:JuMP.Variable}(n::NLOpt,x::Array{T,2},u::Array{T,2})
      if n.s.integrationMethod==:tm; L=size(x)[1]; else; L=size(x)[1]-1; end
      dx=Array(Any,L,n.numStates);
        for st in 1:n.numStates
          println("test \n")
          @show $(eq.args)[st]
          @show $(esc(eq.args))[st]
          j=1
          #@show $(eq)[st]
          println("test2")
          #println($(eq)[st])
          @show x
        #  dx[:,st]=@NLexpression(n.mdl,[j=1:L], $(eq)[st] )

          for j in 1:L
          #  println($(eq)[j,st])
          println("test3")
          @show $(eq)[st]
          println("test4")

            dx[j,st]=JuMP.NonlinearExpression(n.mdl,JuMP.@processNLExpr(n.mdl, $(eq)[j,st] ))
          end

        end
      return dx
    end
  end
end
=#
