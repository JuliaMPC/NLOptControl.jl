#TODO
# 1) why is it failing when we run it again?
# 4) try implicit method to see if that helps
# 6) try other linear solvers
# 9) allow for functions of state constraints of these so we can calculate them on the fly!
# 11) let the user define the objective function above
# 12) allow user to select from using the IMatrix or quadrature
# 13) make a bool to tell the user to restart julia if they change the model significantly -> or do this for them
# 14) make functionality to easily compare all intergration schemes
# 15) implement method with variable time steps for :tm methods, that is how it was before!
try to register functions with JuMP

#mdl = Model(solver=IpoptSolver(linear_solver = "mumps")) #linear_solver = "ma57"
#d = JuMP.NLPEvaluator(mdl)
#MathProgBase.initialize(d, [:Grad,:Hess, :Jac, :ExprGraph])

#MathProgBase.constr_expr(d,i) #6
#MathProgBase.hesslag_structure(d)
#MathProgBase.jac_structure(d)
#MathProgBase.obj_expr(d)


WARNING: The addition operator has been used on JuMP expressions a large number of times.
This warning is safe to ignore but may indicate that model generation is slower than necessary.
For performance reasons, you should not add expressions in a loop. Instead of x += y, use append!(x,y) to modify x in place.
If y is a single variable, you may also use push!(x, coef, y) in place of x += coef*y.
